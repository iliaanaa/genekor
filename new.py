import psycopg2
import pandas as pd
import urllib.request
import re
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
            RCVaccession TEXT[],
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            last_evaluated DATE,
            protein_pos BIGINT
        );
        """)
        conn.commit()

def extract_HGVS(name: str) -> Dict[str, Optional[str]]:
    """Εξαγωγή HGVS (ίδια με pipeline)"""
    result = {'HGVS_c': None, 'HGVS_p': None}
    if not pd.isna(name) and isinstance(name, str):
        dna_pattern = re.compile(r'(c\.[^*\s]+)')
        dna_star_pattern = re.compile(r'(c\.\*[\d_]+[^\s)]*)')
        protein_pattern = re.compile(r'(p\.[^\s)]+)')
        
        dna_match = dna_pattern.search(name)
        dna_star_match = dna_star_pattern.search(name)
        protein_match = protein_pattern.search(name)
        
        if dna_match:
            result['HGVS_c'] = dna_match.group(1)
        elif dna_star_match:
            result['HGVS_c'] = dna_star_match.group(1)
        if protein_match:
            result['HGVS_p'] = protein_match.group(1)
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
    print("Στήλες διαθέσιμες στο αρχείο:")
    print(df.columns.tolist())
    os.remove("filtered_variants.tsv")

    # Εξαγωγή HGVS
    hgvs_data = df['Name'].apply(extract_HGVS).apply(pd.Series)
    df = pd.concat([df, hgvs_data], axis=1)

    # Προσθήκη protein_pos
    df['protein_pos'] = df['HGVS_p'].apply(extract_protein_pos)

    # Μετατροπή NumberSubmitters σε ακέραιο (αν υπάρχει), αλλιώς μηδέν
    if 'NumberSubmitters' in df.columns:
        df['NumberSubmitters'] = pd.to_numeric(df['NumberSubmitters'], errors='coerce').fillna(0).astype(int)
    else:
        df['NumberSubmitters'] = 0

    # Φιλτράρισμα: να κρατήσουμε μόνο εγγραφές με έγκυρο protein_pos
    df = df[df['protein_pos'].notna()]
    df['protein_pos'] = df['protein_pos'].astype(int)

    # In your process_clinvar_data function, add:
    df = df.rename(columns={
        'ClinicalSignificance': 'ClinicalSignificance',
        # Add other columns if needed
    })

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
        if row['HGVS_p'] in known_pathogenic:
            if row['HGVS_c'] != known_pathogenic[row['HGVS_p']]['dna']:
                criteria.append("PS1")

        # PM5
        protein_pos = row['protein_pos']
        if protein_pos in pathogenic_positions and row['HGVS_p'] not in known_pathogenic:
            criteria.append('PM5')

        # PP5
        if (row['ClinicalSignificance'] == 'Pathogenic' 
                and row.get('NumberSubmitters', 0) > 1 
                and row.get('ConflictingInterpretations', '.') == '.'):
            criteria.append('PP5')

        # BP6
        if (row['ClinicalSignificance'] == 'Benign' 
                and row.get('NumberSubmitters', 0) > 1 
                and row.get('ConflictingInterpretations', '.') == '.'):
            criteria.append('BP6')

        return criteria 

    # Εφαρμογή σε όλο το dataframe
    df['ACMG_criteria'] = df.apply(_apply_criteria, axis=1)
    return df


def variant_assortments(df, ref_gene, ref_c, ref_p=None, ref_pos=None):
    """Εύρεση παρόμοιων μεταλλάξεων (ίδια με pipeline)"""
    same_c = df[(df['GeneSymbol'] == ref_gene) & (df['HGVS_c'] == ref_c)]
    same_p = df[(df['GeneSymbol'] == ref_gene) & (df['HGVS_p'] == ref_p) & (df['HGVS_c'] != ref_c)] if ref_p else pd.DataFrame()
    same_pos = df[(df['GeneSymbol'] == ref_gene) & (df['protein_pos'] == ref_pos) & (df['HGVS_p'] != ref_p)] if ref_pos else pd.DataFrame()
    return same_c, same_p, same_pos

def split_by_significance(df):
    """Ομαδοποίηση κατά κλινική σημασία (ίδια με pipeline)"""
    return {
        'Pathogenic': df[df['ClinicalSignificance'].str.contains('Pathogenic', na=False)],
        'Benign': df[df['ClinicalSignificance'].str.contains('Benign', na=False)],
        'VUS': df[df['ClinicalSignificance'].str.contains('Uncertain significance', na=False)]
    }

def build_acmg_support_tables(same_c_groups, same_p_groups, same_pos_groups):
    """Κατασκευή πινάκων υποστήριξης ACMG (ίδια με pipeline)"""
    return {
        'PP5': same_c_groups['Pathogenic'][['HGVS_c']],
        'BP6': same_c_groups['Benign'][['HGVS_c']],
        'PS1': same_p_groups['Pathogenic'][['HGVS_p', 'HGVS_c']],
        'PM5': same_pos_groups['Pathogenic'][['HGVS_p', 'HGVS_c', 'protein_pos']]
    }

def mark_acmg_criteria(df, support):
    """Σήμανση variants με βάση τα groups (ίδια με pipeline)"""
    def _determine_criteria(row):
        crit = []
        if row['HGVS_c'] in support['PP5']['HGVS_c'].values:
            crit.append('PP5')
        if row['HGVS_c'] in support['BP6']['HGVS_c'].values:
            crit.append('BP6')
        if pd.notna(row['HGVS_p']) and pd.notna(row['HGVS_c']):
            ps1_matches = support['PS1'][
                (support['PS1']['HGVS_p'] == row['HGVS_p']) & 
                (support['PS1']['HGVS_c'] != row['HGVS_c'])]
            if not ps1_matches.empty:
                crit.append('PS1')
        if pd.notna(row['protein_pos']) and pd.notna(row['HGVS_p']):
            pm5_matches = support['PM5'][
                (support['PM5']['protein_pos'] == row['protein_pos']) & 
                (support['PM5']['HGVS_p'] != row['HGVS_p'])]
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
                acmg_criteria, conflicting_interpretations, RCVaccession, protein_pos
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
                RCVaccession = EXCLUDED.RCVaccession,
                protein_pos = EXCLUDED.protein_pos,
                last_updated = CURRENT_TIMESTAMP;
            """, (
                row['VariationID'], row['GeneSymbol'], row['HGVS_c'], row['HGVS_p'],
                row['molecular_consequence'], row['ClinicalSignificance'], row['ReviewStatus'],
                row['PhenotypeList'], row['Assembly'], row['Chromosome'], row['Start'], row['Stop'],
                row['ReferenceAllele'], row['AlternateAllele'], Json(row['acmg_criteria']),
                Json(row['ConflictingInterpretations']), row['RCVaccession'], row['protein_pos']
            ))
        conn.commit()
def group_based_acmg(row: pd.Series, df: pd.DataFrame) -> str:
    """
    Επιστρέφει grouping-based ACMG κριτήρια (PS1, PM5) για το δεδομένο variant.
    Βασίζεται σε: ίδιο protein position, ίδιο protein HGVS ή ίδιο cDNA HGVS με άλλα κλινικά σημαντικά variants.
    """
    criteria = []

    # PS1: Ίδιο πρωτεϊνικό αποτέλεσμα (HGVS_p) αλλά διαφορετικό HGVS_c, με Pathogenic απόφαση
    same_p = df[(df['HGVS_p'] == row['HGVS_p']) & (df['HGVS_c'] != row['HGVS_c'])]
    if not same_p.empty and any(same_p['ClinicalSignificance'].str.contains('Pathogenic')):
        criteria.append('PS1')

    # PM5: Ίδια θέση (protein_pos), διαφορετικό HGVS_p, αλλά κάποια απόφαση Pathogenic
    same_pos = df[(df['protein_pos'] == row['protein_pos']) & (df['HGVS_p'] != row['HGVS_p'])]
    if not same_pos.empty and any(same_pos['ClinicalSignificance'].str.contains('Pathogenic')):
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
    if row['HGVS_p'] in known_pathogenic:
        if row['HGVS_c'] != known_pathogenic[row['HGVS_p']]['dna']:
            criteria.append("PS1")

    # PM5
    protein_pos = row.get('protein_pos')
    if protein_pos in pathogenic_positions and row['HGVS_p'] not in known_pathogenic:
        criteria.append('PM5')

    # PP5
    if row.get('ClinicalSignificance') == 'Pathogenic' and row.get('NumberSubmitters', 0) > 1 and row.get('ConflictingInterpretations') == '.':
        criteria.append('PP5')

    # BP6
    if row.get('ClinicalSignificance') == 'Benign' and row.get('NumberSubmitters', 0) > 1 and row.get('ConflictingInterpretations') == '.':
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
    return {
        'Pathogenic': df[df['ClinicalSignificance'].str.contains('Pathogenic', na=False)],
        'Benign': df[df['ClinicalSignificance'].str.contains('Benign', na=False)],
        'VUS': df[df['ClinicalSignificance'].str.contains('Uncertain significance', na=False)]
    }

def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)
    
    try:
        print("Ξεκίνημα script...")
        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)
        
        df_final = process_clinvar_data(variant_gz)
        df_final = apply_acmg_criteria(df_final)

        # Add molecular consequence information
        mol_consequence = df_final['Name'].apply(categorize_variant_name).apply(pd.Series)
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

        # Process groups
        same_c_groups = split_by_significance(same_c)
        same_p_groups = split_by_significance(same_p)
        same_pos_groups = split_by_significance(same_pos)
        support_tables = build_acmg_support_tables(same_c_groups, same_p_groups, same_pos_groups)

        # Mark and combine criteria
        df_final = mark_acmg_criteria(df_final, support_tables)
        df_final["acmg_criteria"] = df_final.apply(lambda row: apply_acmg_criteria_to_row(row), axis=1)

        # Combine all ACMG criteria
        df_final["acmg_combined_criteria"] = df_final.apply(
            lambda row: "; ".join(sorted(
                set(row['acmg_criteria'] + [
                    x.strip() for x in row['acmg_from_grouping'].split(";") if x
                ])
            )),
            axis=1
        )

        # Εισαγωγή στη βάση
        insert_to_database(conn, df_final)
        print("Ολοκληρώθηκε η επεξεργασία!")

    finally:
        conn.close()

if __name__ == "__main__":
    main()