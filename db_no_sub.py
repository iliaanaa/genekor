import pandas as pd
import psycopg2
import os
import re
import gzip
import subprocess
import shutil
import urllib.request
import json
from datetime import datetime
from typing import Dict, List, Optional
import re
import numpy as np
from psycopg2.extras import Json


# --- Ρυθμίσεις ---
CLINVAR_VARIANT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_SUBMISSION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz" #μορφή gzip
GENE_FILTER = "KLHL10"
ASSEMBLY_FILTER = "GRCh38"  #Επιλέγει την έκδοση GRCh38 του γονιδιώματος

#Ρυθμίσεις Βάσης Δεδομένων
DB_CONFIG = {
    "dbname": "clinvar_db",
    "user": "ilianam",
    "password": "genekor123!",
    "host": "localhost",
    "port": 5432
}

#Βοηθητικές Συναρτήσεις
def download_file(url: str, output_path: str) -> None:
    """Κατέβασμα και αποσυμπίεση αρχείου"""
    print(f"Κατέβασμα {url}...")
    urllib.request.urlretrieve(url, output_path)
    
    if output_path.endswith(".gz"):
        print(f"Αποσυμπίεση {output_path}...")
        with gzip.open(output_path, 'rb') as f_in:
            with open(output_path.replace(".gz", ""), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
#grep h zcat gia na mh skaei
def create_tables(conn: psycopg2.extensions.connection) -> None:
    """Δημιουργία πινάκων στη βάση"""
    with conn.cursor() as cur:
        # Κύριος πίνακας μεταλλάξεων
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

#ACMG Criteria
def apply_acmg_criteria(row: pd.Series) -> List[str]:
    criteria = []

    known_pathogenic = {
        'p.Arg504Gly': {'dna': 'c.1510A>T', 'significance': 'Pathogenic'},
        'p.Trp41*': {'dna': 'c.123G>A', 'significance': 'Pathogenic'}
    }

    trusted_submitters = {'ClinVar', 'ExpertLab'}
    pathogenic_positions = {41, 504}

    # PS1
    if row['HGVS_p'] in known_pathogenic:
        if row['HGVS_c'] != known_pathogenic[row['HGVS_p']]['dna']:
            criteria.append("PS1")

    # PM5
    #protein_pos = int(''.join(filter(str.isdigit, row['Protein_variant']))) if pd.notna(row['Protein_variant']) else None
    protein_pos = extract_protein_pos(row['Protein_variant'])  if pd.notna(row['Protein_variant']) else None
    if protein_pos in pathogenic_positions and row['HGVS_p'] not in known_pathogenic:
        criteria.append('PM5')

    # PP5 / BP6
    if 'Submitter' in row and pd.notna(row['Submitter']) and row['Submitter'] in trusted_submitters:
        if row['ClinicalSignificance'] == 'Pathogenic':
            criteria.append('PP5')
        elif row['ClinicalSignificance'] == 'Benign':
            criteria.append('BP6')

    return {
        'acmg_criteria': criteria,
        'protein_pos': protein_pos
    }





           # if row['Submitter'] in trusted_submitters:
    #if row['ClinicalSignificance'] == 'Pathogenic':
       # criteria.append('PP5')
    #elif row['ClinicalSignificance'] == 'Benign':
     #   criteria.append('BP6')


    
    # Υπολογισμός του βαθμού σύγκλισης/διαφωνίας
    def calculate_conflict_score(conflict_data):
        if not conflict_data or len(conflict_data['interpretations']) <= 1:
            return 0
        total = sum(conflict_data['interpretations'].values())
        max_agree = max(conflict_data['interpretations'].values())
        return 1 - (max_agree / total)
    
    df_brca['conflict_score'] = (
        df_brca['conflicting_interpretations']
        .apply(calculate_conflict_score)
    )
    
    # Προσθήκη στήλης που δείχνει αν υπάρχουν conflicting interpretations
    df_brca['has_conflicts'] = (
        df_brca['conflict_score'] > 0.2
    )  # Ορίζουμε threshold 20% διαφωνία
    
    # Υπολογισμός ACMG criteria
    df_brca['acmg_criteria'] = df_brca.apply(apply_acmg_criteria, axis=1)
    
    # Δημιουργία λίστας RCV accessions
    df_brca['RCV_accession'] = df_brca['RCVaccession'].str.split('|')
    
    return df_brca

def parse_conflict(val):
    """Παίρνει την τιμή από τη στήλη conflicting_interpretations και επιστρέφει καθαρή λίστα"""
    if isinstance(val, bool):
        return [str(val)]
    if pd.isna(val):
        return []
    if isinstance(val, str):
        try:
            # Αν είναι ήδη JSON list σε μορφή string
            return json.loads(val)
        except json.JSONDecodeError:
            return [val]
    if isinstance(val, list):
        return val
    return [str(val)]  # τελευταία γραμμή άμυνας


def parse_rcv(value):
    """Ασφαλής μετατροπή RCVaccession σε λίστα από string"""
    if pd.isna(value):
        return None
    try:
        if isinstance(value, str):
            loaded = json.loads(value)
        else:
            loaded = value
        return [str(x) for x in loaded] if isinstance(loaded, list) else [str(loaded)]
    except Exception:
        return [str(value)]


def safe_jsonb(value):
    """
    Μετατρέπει την τιμή σε κατάλληλο αντικείμενο για jsonb insert.
    Αν είναι boolean, None, list ή string με 'True'/'False', χειρίζεται σωστά.
    """
    if value is None or pd.isna(value):
        return json.dumps(None)

    if isinstance(value, bool):
        return json.dumps(value)

    if isinstance(value, list):
        # Αν έχει μόνο ένα στοιχείο που είναι string "True"/"False", το χειριζόμαστε
        if len(value) == 1 and value[0] in ['True', 'False']:
            return json.dumps(value[0] == 'True')
        return json.dumps(value)

    if isinstance(value, str):
        val = value.strip().lower()
        if val in ['true', 'false']:
            return json.dumps(val == 'true')
        try:
            return json.dumps(json.loads(value))  # π.χ. stringified λίστα
        except:
            return json.dumps([value])  # fallback σε λίστα με το string

    return json.dumps(value)

def insert_to_database(conn: psycopg2.extensions.connection, df: pd.DataFrame) -> None:
    """Εισαγωγή δεδομένων στη βάση με αναλυτικό debugging"""
    with conn.cursor() as cur:
        for index, row in df.iterrows():
            try:
                # Prepare conflicting_interpretations and acmg_criteria
                conflict_value = row.get('conflicting_interpretations')
                conflict_list = parse_conflict(conflict_value)  # always returns a list
                acmg_criteria = row.get('acmg_criteria') or []

                print(f"\nΓραμμή {index}:")
                print(f"VariationID: {row['VariationID']} (τύπος {type(row['VariationID'])})")
                print(f"conflicting_interpretations: {conflict_value} (τύπος {type(conflict_value)})")
                print(f"conflict_list μετά parse_conflict: {conflict_list} (τύπος {type(conflict_list)})")
                print(f"acmg_criteria: {acmg_criteria} (τύπος {type(acmg_criteria)})")
                print(f"protein_pos: {row['protein_pos']} (τύπος {type(row['protein_pos'])})")

                cur.execute("""
                INSERT INTO gene_variants (
                    variation_id, gene_symbol, transcript_id, molecular_consequence,
                    hgvs_c, hgvs_p, clinical_significance, review_status, phenotype_list,
                    assembly, chromosome, start_pos, end_pos,
                    reference_allele, alternate_allele, acmg_criteria,
                    conflicting_interpretations, rcvaccession, protein_pos
                ) VALUES (
                    %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
                )
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
                    transcript_id = EXCLUDED.transcript_id,
                    last_updated = CURRENT_TIMESTAMP,
                    protein_pos = EXCLUDED.protein_pos;
                """, (
                    row['VariationID'],
                    row['GeneSymbol'],
                    row['transcript_id'],
                    row['molecular_consequence'],
                    row['HGVS_c'],
                    row['HGVS_p'],
                    row['ClinicalSignificance'],
                    row['ReviewStatus'],
                    row['PhenotypeList'],
                    row['Assembly'],
                    row['Chromosome'],
                    row['Start'],
                    row['Stop'],
                    row['ReferenceAllele'],
                    row['AlternateAllele'],
                    Json(acmg_criteria),
                    Json(conflict_list),
                    row['RCVaccession'],
                    None if pd.isna(row['protein_pos']) else int(row['protein_pos'])  # ✅ εδώ η διόρθωση
                ))

            except Exception as e:
                print(f"\nΣφάλμα στη γραμμή {index}")
                print(row[['VariationID', 'conflicting_interpretations']])
                print(f"Τύπος conflicting_interpretations: {type(conflict_value)}")
                print(f"Τύπος conflict_list: {type(conflict_list)}")
                print(f"Σφάλμα: {e}")
                raise
        conn.commit()

    
def extract_HGVS(name: str) -> dict:
    """
    Εξάγει HGVS.c και HGVS.p από το πεδίο name χρησιμοποιώντας τα συγκεκριμένα regex patterns
    """
    result = {
        'HGVS_c': None,
        'HGVS_p': None,
    }
    
    if not pd.isna(name) and isinstance(name, str):
        dna_pattern = re.compile(r'(c\.[^*\s]+)')  # DNA μεταλλάξεις (c.) χωρίς *
        dna_star_pattern = re.compile(r'(c\.\*[\d_]+[^\s)]*)')  # DNA μεταλλάξεις με * (π.χ. c.*103_*106del)
        protein_pattern = re.compile(r'(p\.[^\s)]+)')  # Πρωτεϊνικές μεταλλάξεις (p.)
        
        # Αναζήτηση για DNA μεταλλάξεις
        dna_match = dna_pattern.search(name)
        dna_star_match = dna_star_pattern.search(name)
        
        if dna_match:
            result['HGVS_c'] = dna_match.group(1)
        elif dna_star_match:
            result['HGVS_c'] = dna_star_match.group(1)
        
        # Αναζήτηση για πρωτεϊνικές μεταλλάξεις
        protein_match = protein_pattern.search(name)
        if protein_match:
            result['HGVS_p'] = protein_match.group(1)
    
    return result

# Παράδειγμα DataFrame (βάλε το δικό σου)
df = pd.DataFrame({
    'Name': [
        'c.123A>T',
        'p.Arg117His',
        'c.*103_*106del',
        'p.Gly12Asp c.35G>A',
        'unknown_variant',
        None
    ]
})

# Εφαρμογή συνάρτησης και προσθήκη στηλών
df_HGVS = df['Name'].apply(extract_HGVS).apply(pd.Series)
df = pd.concat([df, df_HGVS], axis=1)

print(df)


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



def process_clinvar_data(variant_gz_path: str) -> pd.DataFrame:
    """
    Επεξεργασία δεδομένων ClinVar από .gz με zcat και grep (για χαμηλή χρήση μνήμης).
    Διατηρεί τη γραμμή τίτλων (header).
    """
    import subprocess

    print("Φιλτράρισμα δεδομένων με zcat + grep...")
    
    # Δημιουργεί προσωρινό αρχείο με header + filtered rows
    grep_cmd = (
        f"zcat {variant_gz_path} | head -n 1 > filtered_variants.tsv && "
        f"zcat {variant_gz_path} | grep -E 'KLHL10.*GRCh38' >> filtered_variants.tsv"
    )
    subprocess.run(grep_cmd, shell=True, check=True)

    print("Φόρτωση φιλτραρισμένων δεδομένων...")
    df = pd.read_csv("filtered_variants.tsv", sep='\t', low_memory=False)

    print("Κατηγοριοποίηση μεταλλάξεων...")
    #df['VariantName_analysis'] = df['Name'].apply(extract_HGVS)
     # --- 1. Κατηγοριοποίηση ονόματος μετάλλαξης (molecular_consequence, DNA, Protein, Other) ---
    df_variants = df['Name'].apply(categorize_variant_name).apply(pd.Series)
    df = pd.concat([df, df_variants], axis=1)

     # --- 2. Εξαγωγή HGVS_c και HGVS_p ---
    df_hgvs = df['Name'].apply(extract_HGVS).apply(pd.Series)
    df = pd.concat([df, df_hgvs], axis=1)
       # --- 3. Εξαγωγή transcript ID ---
    df['transcript_id'] = df['Name'].apply(extract_transcript_id)

      # --- 3. Βιολογικός τύπος μετάλλαξης (missense, deletion, κ.λπ.) ---
    df['molecular_consequence'] = df.apply(lambda row: determine_variant_type(row['HGVS_p'], row['HGVS_c']), axis=1)
    # --- 4. Αντιγραφή HGVS στις στήλες DNA_variant / Protein_variant ---
    df['DNA_variant'] = df['HGVS_c']
    df['Protein_variant'] = df['HGVS_p']
    # --- Υπολογισμός ACMG κριτηρίων ---
    df_acmg = df.apply(apply_acmg_criteria, axis=1).apply(pd.Series)
    df = pd.concat([df, df_acmg], axis=1)

    #df['DNA_variant'] = df['VariantName_analysis'].apply(lambda x: x['DNA_variant'])
    #df['Protein_variant'] = df['VariantName_analysis'].apply(lambda x: x['Protein_variant'])
    #df['Other_variant'] = df['VariantName_analysis'].apply(lambda x: x['Other_variant'])
    
    # --- 5. Υπολογισμός Other_variant μόνο αν λείπουν και τα δύο ---
    df['Other_variant'] = df.apply(
    lambda row: row['Name'] if pd.isna(row['HGVS_c']) and pd.isna(row['HGVS_p']) else None,
    axis=1
    )
    
 #    # Επιστρέφει dict -> μετατρέπεται σε στήλες
   # df_acmg = df.apply(apply_acmg_criteria, axis=1).apply(pd.Series)

  #  # Προσθήκη των στηλών στο αρχικό df
 #   df['acmg_criteria'] = df_acmg['acmg_criteria']
#    df['protein_pos'] = df_acmg['protein_pos']


    df['protein_pos'] = df['Protein_variant'].apply(extract_protein_pos)

    df_invalid = df[df['protein_pos'] > 9223372036854775807]
    print(df_invalid[['Name', 'Protein_variant', 'protein_pos']])


    df['conflicting_interpretations'] = df['ClinicalSignificance'].str.contains('conflicting', case=False, na=False)
    #df.drop(columns=['VariantName_analysis'], inplace=True)
    df['RCVaccession'] = df['RCVaccession'].fillna('').apply(lambda x: x.split('|') if x else [])

    # Διαγραφή προσωρινού αρχείου
    os.remove("filtered_variants.tsv")
    assert 'acmg_criteria' in df.columns, "acmg_criteria ΔΕΝ προστέθηκε!"

    return df

def extract_transcript_id(name: str)->str:
    if pd.isna(name) or not isinstance(name,str):
        return None
    
    transcript_pattern=re.compile(r'(?<!\w)(NM_\d{5,}(?:\.\d{1,2})?)(?=[(])')

    match = transcript_pattern.search(name)
    
    return match.group(1) if match else None

    '''
    (?<!\w) - Negative lookbehind: Βεβαιώνεται ότι δεν υπάρχει word character πριν

([NXY]M_\d{5,}(?:\.\d{1,2})?) - Κύρια ομάδα:

(?<!\w) Negative lookbehind: Βεβαιώνεται ότι το NM_ δεν προηγείται από άλλο word character (π.χ. γράμμα, αριθμό ή _).

NM_ - Ταιριάζει ακριβώς το πρόθεμα των RefSeq mRNA transcripts.

\d{5,} - Τουλάχιστον 5 ψηφία (οι πραγματικοί αριθμοί transcript είναι συνήθως 5-6 ψηφία)

(?:\.\d{1,2})? - Προαιρετική έκδοση (1-2 ψηφία)

(?=[(]) - Positive lookahead: Πρέπει να ακολουθείται από (
'''

def determine_variant_type(hgvs_p: str, hgvs_c: str) -> str:
    """
    Καθορίζει τον τύπο της μετάλλαξης βάσει των HGVS προσδιορισμών.
    Επιστρέφει ένα από:
    - frameshift, nonsense, deletion, duplication, insertion,
    - missense, synonymous, protein_other,
    - splice_site_essential, splice_region, 5'UTR, 3'UTR, non_coding
    - unknown
    """
    if pd.notna(hgvs_p) and isinstance(hgvs_p, str):
        hgvs_p = hgvs_p.strip()
        if "fs" in hgvs_p:
            return "frameshift"
        elif "*" in hgvs_p:
            return "nonsense"
        elif "del" in hgvs_p:
            return "deletion"
        elif "dup" in hgvs_p:
            return "duplication"
        elif "ins" in hgvs_p:
            return "insertion"
        elif hgvs_p == "p.=":
            return "synonymous"
        elif re.match(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", hgvs_p):
            return "missense"
        else:
            return "protein_other"
    
    if pd.notna(hgvs_c) and isinstance(hgvs_c, str):
        hgvs_c = hgvs_c.strip()
        if re.search(r"\+\d+|\-\d+", hgvs_c):
            if re.search(r"\+1\+2|\-1\-2", hgvs_c):
                return "splice_site_essential"
            else:
                return "splice_region"
        elif hgvs_c.startswith("c.-"):
            return "5'UTR"
        elif "*" in hgvs_c:
            return "3'UTR"
    
    return "unknown"




def extract_protein_pos(protein_variant):
    """
    Εξάγει τον αριθμό θέσης από HGVS.p string, π.χ. από 'p.Arg167His' -> 167
    """
    try:
        match = re.search(r'[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}', protein_variant)
        if match:
            pos = int(match.group(1))
            # ✅ Φιλτράρουμε παράλογες θέσεις (π.χ. πάνω από 10000)
            if pos < 10000:
                return pos
    except:
        pass
    return np.nan



def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)
    
    try:
        print("Ξεκίνημα script...")  # Για επιβεβαίωση ότι τρέχει
        # Λήψη αρχείων (ΜΟΝΟ τα .gz, χωρίς αποσυμπίεση)
        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)
        
        # Επεξεργασία ΑΜΕΣΑ από το .gz
        df_final = process_clinvar_data(variant_gz)
        print(df.columns)
        print(df.head())
        print("Columns in DataFrame before insert:", df.columns.tolist())

          # Εκτύπωση εγγραφών με πολύ μεγάλο protein_pos (> 2000)
        high_pos_df = df_final[df_final['protein_pos'] > 2000]
        if not high_pos_df.empty:
            print(f"Βρέθηκαν {len(high_pos_df)} εγγραφές με protein_pos > 2000:")
            print(high_pos_df[['Name', 'HGVS_p', 'protein_pos']].to_string(index=False))
        else:
            print("Καμία εγγραφή με protein_pos > 2000.")

        invalid_protein_pos = df_final[df_final['protein_pos'] > 9223372036854775807]
        if not invalid_protein_pos.empty:
            print("Προσοχή: Βρέθηκαν τιμές protein_pos εκτός ορίου και παραλείφθηκαν:")
            print(invalid_protein_pos[['Name', 'Protein_variant', 'protein_pos']])



        # Φιλτράρισμα τιμών εκτός BIGINT εύρους
        df_final = df_final[df_final['protein_pos'].isna() | (df_final['protein_pos'] <= 9223372036854775807)]
        # Επιβεβαίωση
        print("Πλήθος εγγραφών πριν insert:", len(df_final))

        bad_rows = df_final[df_final['protein_pos'] > 9223372036854775807]

        if not bad_rows.empty:
            print("Εγγραφές με υπερβολικά μεγάλο protein_pos:")
            print(bad_rows[['Name', 'HGVS_p', 'Protein_variant', 'protein_pos']])
        else:
            print("Καμία εγγραφή με υπερβολικά μεγάλο protein_pos.")


        print("Τύποι δεδομένων:")
        print(df_final.dtypes)

        # Δες τις max τιμές σε αριθμητικά πεδία
        print("Μέγιστες αριθμητικές τιμές:")
        print(df_final.select_dtypes(include='number').max())
        

        bigint_limit = 9223372036854775807
        too_big = df_final.select_dtypes(include='number') > bigint_limit

        if too_big.any().any():
            print("Βρέθηκαν τιμές που ξεπερνούν το όριο του BIGINT:")
            for col in too_big.columns:
                if too_big[col].any():
                    print(f" - {col}")
                    print(df_final.loc[too_big[col], [col, 'Name']])
        else:
            print(" Όλα τα αριθμητικά πεδία είναι εντός ορίων.")


        # Εισαγωγή στη βάση
        insert_to_database(conn, df_final)
        
    except Exception as e:
        print(f"Σφάλμα: {e}")
    finally:
        if os.path.exists(variant_gz):
            os.remove(variant_gz)
        conn.close()
        
test_vals = [False, True, "['Conflicting interpretations']", None, float('nan')]
for val in test_vals:
    print(f"Είσοδος: {val} -> Έξοδος: {parse_conflict(val)}")

if __name__ == "__main__":
    main()