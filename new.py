import psycopg2
import pandas as pd
import urllib.request
import re
import traceback
import os
from typing import Dict, List, Optional, Tuple
from psycopg2.extras import Json
import subprocess

# --- Ρυθμίσεις ---
CLINVAR_VARIANT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
DB_CONFIG = {
    "dbname": "clinvar_db",
    "user": "ilianam",
    "password": "genekor123!",
    "host": "localhost",
    "port": 5432
}

# --- Βοηθητικές Συναρτήσεις (ΑΚΡΙΒΩΣ όπως στο pipeline) ---
def create_tables(conn: psycopg2.extensions.connection) -> None:
    """Δημιουργία πινάκων βάσης (ίδια με pipeline)"""
    with conn.cursor() as cur:
        cur.execute("""
        CREATE TABLE IF NOT EXISTS gene_variants (
            variation_id BIGINT PRIMARY KEY,
            gene_symbol TEXT NOT NULL,
            transcript_id TEXT,
            hgvs_c TEXT,
            hgvs_p TEXT,
            molecular_consequence TEXT,
            clinical_significance TEXT,
            review_status TEXT,
            phenotype_list TEXT,
            assembly TEXT NOT NULL,
            chromosome TEXT,
            start_pos INTEGER,
            end_pos INTEGER,
            reference_allele TEXT,
            alternate_allele TEXT,
            acmg_criteria JSONB,
            conflicting_interpretations JSONB,
            rcvaccession TEXT[],
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            last_evaluated DATE,
            protein_pos BIGINT
        );
        """)
        conn.commit()

def extract_HGVS(name: str) -> Dict[str, Optional[str]]:
    """Εξαγωγή HGVS (ίδια με pipeline)"""
    result = {'hgvs_c': None, 'hgvs_p': None}
    if not pd.isna(name) and isinstance(name, str):
        dna_pattern = re.compile(r'(c\.[^*\s]+)')
        dna_star_pattern = re.compile(r'(c\.\*[\d_]+[^\s)]*)')
        protein_pattern = re.compile(r'(p\.[^\s)]+)')
        
        dna_match = dna_pattern.search(name)
        dna_star_match = dna_star_pattern.search(name)
        protein_match = protein_pattern.search(name)
        
        if dna_match:
            result['hgvs_c'] = dna_match.group(1)
        elif dna_star_match:
            result['hgvs_c'] = dna_star_match.group(1)
        if protein_match:
            result['hgvs_p'] = protein_match.group(1)
    return result

def extract_protein_pos(hgvs_p: str) -> Optional[int]:
    """Εξαγωγή θέσης πρωτεΐνης (ίδια με pipeline)"""
    if not hgvs_p:
        return None
    match = re.search(r'[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}', hgvs_p)
    return int(match.group(1)) if match else None

def process_clinvar_data(variant_gz_path: str) -> pd.DataFrame:
    """Επεξεργασία ClinVar με zcat + grep (ίδια με pipeline)"""
    print("Φιλτράρισμα δεδομένων με zcat + grep...")
    grep_cmd = (
        f"zcat {variant_gz_path} | head -n 1 > filtered_variants.tsv && "
        f"zcat {variant_gz_path} | grep -E 'KLHL10.*GRCh38' >> filtered_variants.tsv"
    )
    subprocess.run(grep_cmd, shell=True, check=True)

    print("Φόρτωση δεδομένων...")
    df = pd.read_csv("filtered_variants.tsv", sep='\t', low_memory=False)

    # Normalize column names
    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(r"[^\w]+", "_", regex=True)
        .str.strip("_")
    )

    # Rename clinical significance column to your preferred format
    if 'clinicalsignificance' in df.columns:
        df = df.rename(columns={'clinicalsignificance': 'clinical_significance'})

    print("Στήλες διαθέσιμες στο αρχείο:")
    print(df.columns.tolist())

    # Clean up temporary file
    os.remove("filtered_variants.tsv")

    # Εξαγωγή HGVS
    hgvs_data = df['name'].apply(extract_HGVS).apply(pd.Series)

    # Normalize HGVS subcolumn names
    hgvs_data.columns = (
        hgvs_data.columns
        .str.strip()
        .str.lower()
        .str.replace(r"[^\w]+", "_", regex=True)
        .str.strip("_")
    )

    # Συνένωση HGVS δεδομένων με το αρχικό df
    df = pd.concat([df, hgvs_data], axis=1)

    # Προσθήκη protein_pos (αν υπάρχει hgvs_p)
    if 'hgvs_p' in df.columns:
        df['protein_pos'] = df['hgvs_p'].apply(extract_protein_pos)
    else:
        df['protein_pos'] = None

    # Number of submitters
    if 'numbersubmitters' in df.columns:
        df['numbersubmitters'] = pd.to_numeric(df['numbersubmitters'], errors='coerce').fillna(0).astype(int)
    else:
        df['numbersubmitters'] = 0

    # Φιλτράρισμα: κρατάμε μόνο εγγραφές με έγκυρο protein_pos
    df = df[df['protein_pos'].notna()]
    df['protein_pos'] = df['protein_pos'].astype(int)

    return df


def apply_acmg_criteria(df: pd.DataFrame) -> pd.DataFrame:
    """Εφαρμογή ACMG κριτηρίων (ΑΚΡΙΒΩΣ όπως στο pipeline)"""
    known_pathogenic = {
        'p.Arg504Gly': {'dna': 'c.1510A>T', 'significance': 'Pathogenic'},
        'p.Trp41*': {'dna': 'c.123G>A', 'significance': 'Pathogenic'}
    }
    trusted_submitters = {'ClinVar', 'ExpertLab'}
    pathogenic_positions = {41, 504}

    def _apply_criteria(row):
        criteria = []
        
        # PS1
        if row['hgvs_p'] in known_pathogenic:
            if row['hgvs_c'] != known_pathogenic[row['hgvs_p']]['dna']:
                criteria.append("PS1")

        # PM5
        protein_pos = row['protein_pos']
        if protein_pos in pathogenic_positions and row['hgvs_p'] not in known_pathogenic:
            criteria.append('PM5')

        # PP5
        if (row['clinical_significance'] == 'Pathogenic' 
                and row.get('NumberSubmitters', 0) > 1 
                and row.get('conflictinginterpretations', '.') == '.'):
            criteria.append('PP5')

        # BP6
        if (row['clinical_significance'] == 'Benign' 
                and row.get('NumberSubmitters', 0) > 1 
                and row.get('conflictinginterpretations', '.') == '.'):
            criteria.append('BP6')

        return criteria 

    # Εφαρμογή σε όλο το dataframe
    df['ACMG_criteria'] = df.apply(_apply_criteria, axis=1)
    return df


def variant_assortments(df, ref_gene, ref_c, ref_p=None, ref_pos=None):
    """Εύρεση παρόμοιων μεταλλάξεων"""
    # Make sure we're working with a copy that has all columns
    same_c = df[(df['genesymbol'] == ref_gene) & (df['hgvs_c'] == ref_c)].copy()
    same_p = df[(df['genesymbol'] == ref_gene) & (df['hgvs_p'] == ref_p) & (df['hgvs_c'] != ref_c)].copy() if ref_p else pd.DataFrame()
    same_pos = df[(df['genesymbol'] == ref_gene) & (df['protein_pos'] == ref_pos) & (df['hgvs_p'] != ref_p)].copy() if ref_pos else pd.DataFrame()
    return same_c, same_p, same_pos

def split_by_significance(df):
    """Ομαδοποίηση κατά κλινική σημασία"""
    if df.empty or 'clinical_significance' not in df.columns:
        print("Warning: Empty DataFrame or missing clinical_significance column")
        print("Available columns:", df.columns.tolist())
        return {
            'Pathogenic': pd.DataFrame(),
            'Benign': pd.DataFrame(),
            'VUS': pd.DataFrame()
        }
    
    return {
        'Pathogenic': df[df['clinical_significance'].str.contains('Pathogenic', na=False)],
        'Benign': df[df['clinical_significance'].str.contains('Benign', na=False)],
        'VUS': df[df['clinical_significance'].str.contains('Uncertain significance', na=False)]
    }

def build_acmg_support_tables(same_c_groups, same_p_groups, same_pos_groups):
    """Κατασκευή πινάκων υποστήριξης ACMG (με έλεγχο ύπαρξης δεδομένων και στηλών)"""

    def safe_select(df, columns):
        if df is None or df.empty:
            return pd.DataFrame(columns=columns)
        # Έλεγχος αν όλες οι στήλες υπάρχουν
        missing_cols = [col for col in columns if col not in df.columns]
        if missing_cols:
            return pd.DataFrame(columns=columns)
        return df[columns]

    pp5_df = safe_select(same_c_groups.get('Pathogenic', pd.DataFrame()), ['hgvs_c'])
    bp6_df = safe_select(same_c_groups.get('Benign', pd.DataFrame()), ['hgvs_c'])
    ps1_df = safe_select(same_p_groups.get('Pathogenic', pd.DataFrame()), ['hgvs_p', 'hgvs_c'])
    pm5_df = safe_select(same_pos_groups.get('Pathogenic', pd.DataFrame()), ['hgvs_p', 'hgvs_c', 'protein_pos'])

    return {
        'PP5': pp5_df,
        'BP6': bp6_df,
        'PS1': ps1_df,
        'PM5': pm5_df
    }


def mark_acmg_criteria(df, support):
    """Σήμανση variants με βάση τα groups (ίδια με pipeline)"""
    def _determine_criteria(row):
        crit = []
        if row['hgvs_c'] in support['PP5']['hgvs_c'].values:
            crit.append('PP5')
        if row['hgvs_c'] in support['BP6']['hgvs_c'].values:
            crit.append('BP6')
        if pd.notna(row['hgvs_p']) and pd.notna(row['hgvs_c']):
            ps1_matches = support['PS1'][
                (support['PS1']['hgvs_p'] == row['hgvs_p']) & 
                (support['PS1']['hgvs_c'] != row['hgvs_c'])]
            if not ps1_matches.empty:
                crit.append('PS1')
        if pd.notna(row['protein_pos']) and pd.notna(row['hgvs_p']):
            pm5_matches = support['PM5'][
                (support['PM5']['protein_pos'] == row['protein_pos']) & 
                (support['PM5']['hgvs_p'] != row['hgvs_p'])]
            if not pm5_matches.empty:
                crit.append('PM5')
        return "; ".join(crit) if crit else ""
    
    df['acmg_from_grouping'] = df.apply(_determine_criteria, axis=1)
    return df

def insert_to_database(conn, df):
    """Εισαγωγή στη βάση (ίδια με pipeline)"""
    with conn.cursor() as cur:
        for _, row in df.iterrows():
            cur.execute("""
            INSERT INTO gene_variants (
                variation_id, gene_symbol, hgvs_c, hgvs_p, molecular_consequence,
                clinical_significance, review_status, phenotype_list, assembly,
                chromosome, start_pos, end_pos, reference_allele, alternate_allele,
                acmg_criteria, conflicting_interpretations, rcvaccession, protein_pos
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (variation_id) DO UPDATE SET
                gene_symbol = EXCLUDED.gene_symbol,
                hgvs_c = EXCLUDED.hgvs_c,
                hgvs_p = EXCLUDED.hgvs_p,
                molecular_consequence = EXCLUDED.molecular_consequence,
                clinical_significance = EXCLUDED.clinical_significance,
                review_status = EXCLUDED.review_status,
                phenotype_list = EXCLUDED.phenotype_list,
                assembly = EXCLUDED.assembly,
                chromosome = EXCLUDED.chromosome,
                start_pos = EXCLUDED.start_pos,
                end_pos = EXCLUDED.end_pos,
                reference_allele = EXCLUDED.reference_allele,
                alternate_allele = EXCLUDED.alternate_allele,
                acmg_criteria = EXCLUDED.acmg_criteria,
                conflicting_interpretations = EXCLUDED.conflicting_interpretations,
                rcvaccession = EXCLUDED.rcvaccession,
                protein_pos = EXCLUDED.protein_pos,
                last_updated = CURRENT_TIMESTAMP;
            """, (
                row['variationid'], row['genesymbol'], row['hgvs_c'], row['hgvs_p'],
                row['molecular_consequence'], row['clinical_significance'], row['reviewstatus'],
                row['phenotypelist'], row['assembly'], row['chromosome'], row['start'], row['stop'],
                row['referenceallele'], row['alternateallele'], Json(row['acmg_criteria']),
                Json(row['conflictinginterpretations']), row['rcvaccession'], row['protein_pos']
            ))
        conn.commit()


def group_based_acmg(row: pd.Series, df: pd.DataFrame) -> str:
    """
    Επιστρέφει grouping-based ACMG κριτήρια (PS1, PM5) για το δεδομένο variant.
    Βασίζεται σε: ίδιο protein position, ίδιο protein HGVS ή ίδιο cDNA HGVS με άλλα κλινικά σημαντικά variants.
    """
    criteria = []

    # PS1: Ίδιο πρωτεϊνικό αποτέλεσμα (HGVS_p) αλλά διαφορετικό HGVS_c, με Pathogenic απόφαση
    same_p = df[(df['hgvs_p'] == row['hgvs_p']) & (df['hgvs_c'] != row['hgvs_c'])]
    if not same_p.empty and any(same_p['clinical_significance'].str.contains('Pathogenic')):
        criteria.append('PS1')

    # PM5: Ίδια θέση (protein_pos), διαφορετικό HGVS_p, αλλά κάποια απόφαση Pathogenic
    same_pos = df[(df['protein_pos'] == row['protein_pos']) & (df['hgvs_p'] != row['hgvs_p'])]
    if not same_pos.empty and any(same_pos['clinical_significance'].str.contains('Pathogenic')):
        criteria.append('PM5')

    return "; ".join(criteria)

def apply_acmg_criteria_to_row(row):
    """Επιστρέφει τα ACMG criteria (PS1, PM5, PP5, BP6) για μία γραμμή"""
    known_pathogenic = {
        'p.Arg504Gly': {'dna': 'c.1510A>T', 'significance': 'Pathogenic'},
        'p.Trp41*': {'dna': 'c.123G>A', 'significance': 'Pathogenic'}
    }
    pathogenic_positions = {41, 504}

    criteria = []

    # PS1
    if row['hgvs_p'] in known_pathogenic:
        if row['hgvs_c'] != known_pathogenic[row['hgvs_p']]['dna']:
            criteria.append("PS1")

    # PM5
    protein_pos = row.get('protein_pos')
    if protein_pos in pathogenic_positions and row['hgvs_p'] not in known_pathogenic:
        criteria.append('PM5')

    # PP5
    if row.get('clinical_significance') == 'Pathogenic' and row.get('NumberSubmitters', 0) > 1 and row.get('conflictinginterpretations') == '.':
        criteria.append('PP5')

    # BP6
    if row.get('clinical_significance') == 'Benign' and row.get('NumberSubmitters', 0) > 1 and row.get('conflictinginterpretations') == '.':
        criteria.append('BP6')

    return criteria

def categorize_variant_name(name: str) -> dict:
    """
    Κατηγοριοποιεί το όνομα μετάλλαξης από το πεδίο 'Name' του ClinVar
    σε μεταλλάξεις DNA (c.), πρωτεΐνης (p.) και άλλους τύπους
    """
    result = {
        'molecular_consequence': None,  # 'DNA', 'Protein', 'Other'
        'DNA_variant': None,
        'Protein_variant': None,
        'Other_variant': None
    }
    
    if not pd.isna(name) and isinstance(name, str):
        # Κανονικές εκφράσεις για αναγνώριση τύπων μεταλλάξεων
        dna_pattern = re.compile(r'(c\.[^*\s]+)')  # DNA μεταλλάξεις (c.) που δεν περιέχουν *
        dna_star_pattern = re.compile(r'(c\.\*[\d_]+[^\s)]*)')  # DNA μεταλλάξεις με * (π.χ. c.*103_*106del)
        protein_pattern = re.compile(r'(p\.[^\s)]+)')  # Πρωτεϊνικές μεταλλάξεις (p.)
        
        # Έλεγχος για DNA μεταλλάξεις (κανονικές και με *)
        dna_match = dna_pattern.search(name)
        dna_star_match = dna_star_pattern.search(name)
        
        if dna_match:
            result['molecular_consequence'] = 'DNA'
            result['DNA_variant'] = dna_match.group(1)
        elif dna_star_match:
            result['molecular_consequence'] = 'DNA'
            result['DNA_variant'] = dna_star_match.group(1)
        
        # Έλεγχος για πρωτεϊνικές μεταλλάξεις
        protein_match = protein_pattern.search(name)
        if protein_match and not result['DNA_variant']:  # Προτεραιότητα στις DNA μεταλλάξεις
            result['molecular_consequence'] = 'Protein'
            result['Protein_variant'] = protein_match.group(1)
        
        # Αν δεν βρέθηκε τίποτα από τα παραπάνω
        if not result['molecular_consequence']:
            result['molecular_consequence'] = 'Other'
            result['Other_variant'] = name
    
    return result

def split_by_significance(df):
    """Ομαδοποίηση κατά κλινική σημασία (ίδια με pipeline)"""
    if df.empty:
        # Επιστρέφουμε άδειες DataFrames για κάθε κατηγορία
        return {
            'Pathogenic': df,
            'Benign': df,
            'VUS': df
        }

    return {
        'Pathogenic': df[df['clinical_significance'].str.contains('Pathogenic', na=False)],
        'Benign': df[df['clinical_significance'].str.contains('Benign', na=False)],
        'VUS': df[df['clinical_significance'].str.contains('Uncertain significance', na=False)]
    }


def simplify_clinical_significance(df):
    # Μετατροπή σε μικρά γράμματα για ευκολότερο mapping
    def map_significance(cs):
        if cs is None or cs == '' or cs.lower() == 'not provided':
            return 'not provided'
        cs = cs.lower()
        if 'pathogenic' in cs and 'likely' not in cs:
            return 'pathogenic'
        elif 'likely pathogenic' in cs:
            return 'likely pathogenic'
        elif 'benign' in cs and 'likely' not in cs:
            return 'benign'
        elif 'likely benign' in cs:
            return 'likely benign'
        elif 'uncertain significance' in cs or 'vus' in cs:
            return 'vus'
        elif 'conflicting interpretations' in cs:
            return 'conflicting'
        else:
            return 'other'
    df['clinsigsimple'] = df['clinical_significance'].apply(map_significance)
    return df

def compute_conflictinginterpretations(df):
    conflict_dict = {}
    for vid, group in df.groupby('variationid'):
        unique_sigs = set(group['clinsigsimple'].dropna())
        # Αν έχει πάνω από μία διαφορετική κατηγορία, σημαίνει σύγκρουση
        conflict_dict[vid] = len(unique_sigs) > 1
    df['conflictinginterpretations'] = df['variationid'].map(conflict_dict)
    return df

def convert_pipe_string_to_list(s):
    if pd.isna(s) or s == 'na' or s == '':
        return []
    return s.split('|')


def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)
    
    try:
        print("Ξεκίνημα script...")
        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)
        
        df_final = process_clinvar_data(variant_gz)
        df_final = apply_acmg_criteria(df_final)
        # Εδώ προσθέτουμε τις νέες στήλες clinsigsimple & conflictinginterpretations
        df_final = simplify_clinical_significance(df_final)
        df_final = compute_conflictinginterpretations(df_final)

        # Add molecular consequence information
        mol_consequence = df_final['name'].apply(categorize_variant_name).apply(pd.Series)
        df_final = pd.concat([df_final, mol_consequence], axis=1)

        # Δημιουργία υποστηρικτικών ομάδων grouping
        df_final["acmg_from_grouping"] = df_final.apply(
            lambda row: group_based_acmg(row, df_final), axis=1
        )

        # Create support groups for the entire dataframe
        same_c, same_p, same_pos = variant_assortments(
            df_final,
            ref_gene="KLHL10",
            ref_c=None,
            ref_p=None,
            ref_pos=None
        )

          # Verify groups before processing
        print(f"Group sizes - same_c: {len(same_c)}, same_p: {len(same_p)}, same_pos: {len(same_pos)}")
        
        # Process groups safely
        if same_c.empty:
            print("same_c είναι κενό DataFrame, παραλείπουμε split_by_significance για same_c")
            same_c_groups = {'Pathogenic': same_c, 'Benign': same_c, 'VUS': same_c}
        else:
            same_c_groups = split_by_significance(same_c)

        if same_p.empty:
            print("same_p είναι κενό DataFrame, παραλείπουμε split_by_significance για same_p")
            same_p_groups = {'Pathogenic': same_p, 'Benign': same_p, 'VUS': same_p}
        else:
            same_p_groups = split_by_significance(same_p)

        if same_pos.empty:
            print("same_pos είναι κενό DataFrame, παραλείπουμε split_by_significance για same_pos")
            same_pos_groups = {'Pathogenic': same_pos, 'Benign': same_pos, 'VUS': same_pos}
        else:
            same_pos_groups = split_by_significance(same_pos)

        support_tables = build_acmg_support_tables(same_c_groups, same_p_groups, same_pos_groups)

        # Mark and combine criteria
        print("Marking ACMG criteria...")
        df_final = mark_acmg_criteria(df_final, support_tables)
        df_final["acmg_criteria"] = df_final.apply(lambda row: apply_acmg_criteria_to_row(row), axis=1)

        # Combine all ACMG criteria
        print("Combining ACMG criteria...")
        df_final["acmg_combined_criteria"] = df_final.apply(
            lambda row: "; ".join(sorted(
                set(row['acmg_criteria'] + [
                    x.strip() for x in row['acmg_from_grouping'].split(";") if x
                ])
            )),
            axis=1
        )

        if 'rcvaccession' in df_final.columns:
            df_final['rcvaccession'] = df_final['rcvaccession'].apply(convert_pipe_string_to_list)
        
        # Εισαγωγή στη βάση
        print("Inserting to database...")
        insert_to_database(conn, df_final)
        print("Ολοκληρώθηκε η επεξεργασία!")


    except Exception as e:
        print(f"Error occurred: {str(e)}")
        traceback.print_exc()

    finally:
        conn.close()
        if os.path.exists(variant_gz):
            os.remove(variant_gz)

if __name__ == "__main__":
    main()