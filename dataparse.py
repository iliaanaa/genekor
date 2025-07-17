import pandas as pd
import psycopg2
import os
<<<<<<< HEAD
=======
import re
>>>>>>> 373c3d0 (changes)
import gzip
import shutil
import urllib.request
import json
from datetime import datetime
<<<<<<< HEAD

# --- Ρυθμίσεις ---
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GZ_FILE = "variant_summary.txt.gz"
TSV_FILE = "variant_summary.txt"
GENE_FILTER = "BRCA1"
OUTPUT_CSV = "brca_variants.csv"

# --- Ρυθμίσεις Βάσης Δεδομένων ---
DB_NAME = "clinvar_db"
DB_USER = "postgres"
DB_PASS = "your_password"
DB_HOST = "localhost"
DB_PORT = 5432

# --- Σύνδεση στη Βάση ---
def connect_db():
    return psycopg2.connect(
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASS,
        host=DB_HOST,
        port=DB_PORT
    )

# --- Δημιουργία Πινάκων ---
def create_tables(conn):
    with conn.cursor() as cur:
        # Πίνακας για πληροφορίες έκδοσης
        cur.execute("""
        CREATE TABLE IF NOT EXISTS clinvar_release_info (
            id SERIAL PRIMARY KEY,
            release_version TEXT NOT NULL,
            release_date DATE NOT NULL,
            downloaded_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Κύριος πίνακας μεταλλάξεων
        cur.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            variant_id BIGINT PRIMARY KEY,
            gene_symbol TEXT NOT NULL,
            dna_change TEXT,
            protein_change TEXT,
            clinical_significance TEXT,
            review_status TEXT,
            condition TEXT,
            variant_type TEXT,
            phenotypes TEXT,
            submitter TEXT,
            acmg_criteria JSONB,
            conflicting_interpretations JSONB
=======
from typing import Dict, List, Optional

# --- Ρυθμίσεις ---
CLINVAR_VARIANT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_SUBMISSION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz" #μορφή gzip
GENE_FILTER = "KLHL10"
ASSEMBLY_FILTER = "GRCh38"  #Επιλέγει την έκδοση GRCh38 του γονιδιώματος

#Ρυθμίσεις Βάσης Δεδομένων
DB_CONFIG = {
    "dbname": "clinvar_db",
    "user": "ilianam",
    "password": "wasd1029!",
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

def create_tables(conn: psycopg2.extensions.connection) -> None:
    """Δημιουργία πινάκων στη βάση"""
    with conn.cursor() as cur:
        # Κύριος πίνακας μεταλλάξεων
        cur.execute("""
        CREATE TABLE IF NOT EXISTS brca1_variants (
            variation_id BIGINT PRIMARY KEY,
            gene_symbol TEXT NOT NULL,
            transcript_id TEXT,
            morecular_consequence TEXT,
            hgvs_c TEXT,
            hgvs_p TEXT,
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
            rcv_accessions TEXT[],
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            last_evaluated DATE
>>>>>>> 373c3d0 (changes)
        );
        """)
        conn.commit()

<<<<<<< HEAD
# --- Λήψη έκδοσης ClinVar ---
def get_remote_release_version():
    try:
        with urllib.request.urlopen("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README.txt") as response:
            for line in response.read().decode().splitlines():
                if "Release" in line and "ClinVar" in line:
                    return line.strip()
        return None
    except Exception as e:
        print(f"Σφάλμα κατά τη λήψη έκδοσης: {e}")
        return None

# --- Κατέβασμα και αποσυμπίεση δεδομένων ---
def download_clinvar():
    print("Κατέβασμα δεδομένων ClinVar...")
    urllib.request.urlretrieve(CLINVAR_URL, GZ_FILE)
    print("Αποσυμπίεση...")
    with gzip.open(GZ_FILE, 'rb') as f_in, open(TSV_FILE, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# --- Εξαγωγή και μετασχηματισμός δεδομένων BRCA1 ---
def extract_and_transform():
    print("Φόρτωση και μετασχηματισμός δεδομένων...")
    df = pd.read_csv(TSV_FILE, sep='\t', low_memory=False)
    
    # Φιλτράρισμα για BRCA1
    df = df[df["GeneSymbol"] == GENE_FILTER]
    
    # Μετασχηματισμός στη ζητούμενη μορφή
    df = df[[
        'VariationID', 'GeneSymbol', 'HGVS_c', 'HGVS_p',
        'ClinicalSignificance', 'ReviewStatus', 'Name', 'Type',
        'PhenotypeIDS', 'Submitter'
    ]].rename(columns={
        'VariationID': 'variant_id',
        'GeneSymbol': 'gene_symbol',
        'HGVS_c': 'dna_change',
        'HGVS_p': 'protein_change',
        'ClinicalSignificance': 'clinical_significance',
        'ReviewStatus': 'review_status',
        'Name': 'condition',
        'Type': 'variant_type',
        'PhenotypeIDS': 'phenotypes',
        'Submitter': 'submitter'
    })
    
    # Προσθήκη των ACMG κριτηρίων (θα τα υπολογίσουμε αργότερα)
    df['acmg_criteria'] = None
    df['conflicting_interpretations'] = None
    
    return df

# --- Υπολογισμός ACMG κριτηρίων ---
def calculate_acmg_criteria(df):
    print("Υπολογισμός ACMG κριτηρίων...")
    
    # Λεξικό με γνωστές pathogenic μεταλλάξεις
=======
#ACMG Criteria
def apply_acmg_criteria(row: pd.Series) -> List[str]:
    """Υπολογισμός κριτηρίων ACMG για μια μετάλλαξη"""
    criteria = []
    
    # Γνωστές pathogenic μεταλλάξεις (προσαρμόστε ανάλογα)
>>>>>>> 373c3d0 (changes)
    known_pathogenic = {
        'p.Arg504Gly': {'dna': 'c.1510A>T', 'significance': 'Pathogenic'},
        'p.Trp41*': {'dna': 'c.123G>A', 'significance': 'Pathogenic'}
    }
<<<<<<< HEAD
    
    for index, row in df.iterrows():
        criteria = []
        
        # Έλεγχος για PS1 (ίδιο protein change, διαφορετικό DNA change)
        if row['protein_change'] in known_pathogenic:
            if row['dna_change'] != known_pathogenic[row['protein_change']]['dna']:
                criteria.append("PS1")
        
        # Προσθήκη άλλων κριτηρίων (PM5, PP5, BP6) εδώ...
        
        # Αποθήκευση των κριτηρίων
        df.at[index, 'acmg_criteria'] = json.dumps(criteria) if criteria else None
    
    return df

# --- Αποθήκευση σε CSV ---
def save_to_csv(df, filename=OUTPUT_CSV):
    print(f"Αποθήκευση σε {filename}...")
    
    # Μετατροπή των λιστών/λεξικών σε JSON strings για το CSV
    df_csv = df.copy()
    df_csv['acmg_criteria'] = df_csv['acmg_criteria'].apply(lambda x: json.loads(x) if x else [])
    
    df_csv.to_csv(filename, index=False, encoding='utf-8')

# --- Εισαγωγή στη βάση δεδομένων ---
def insert_to_database(conn, df):
    print("Εισαγωγή δεδομένων στη βάση...")
    with conn.cursor() as cur:
        for _, row in df.iterrows():
            cur.execute("""
            INSERT INTO variants (
                variant_id, gene_symbol, dna_change, protein_change,
                clinical_significance, review_status, condition,
                variant_type, phenotypes, submitter, acmg_criteria,
                conflicting_interpretations
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (variant_id) DO UPDATE SET
                gene_symbol = EXCLUDED.gene_symbol,
                dna_change = EXCLUDED.dna_change,
                protein_change = EXCLUDED.protein_change,
                clinical_significance = EXCLUDED.clinical_significance,
                review_status = EXCLUDED.review_status,
                condition = EXCLUDED.condition,
                variant_type = EXCLUDED.variant_type,
                phenotypes = EXCLUDED.phenotypes,
                submitter = EXCLUDED.submitter,
                acmg_criteria = EXCLUDED.acmg_criteria,
                conflicting_interpretations = EXCLUDED.conflicting_interpretations;
            """, (
                row['variant_id'], row['gene_symbol'], row['dna_change'],
                row['protein_change'], row['clinical_significance'],
                row['review_status'], row['condition'], row['variant_type'],
                row['phenotypes'], row['submitter'], row['acmg_criteria'],
                row['conflicting_interpretations']
            ))
        conn.commit()

# --- Κύρια λειτουργία ---
def main():
    # Σύνδεση στη βάση δεδομένων
    conn = connect_db()
    create_tables(conn)
    
    # Κατέβασμα δεδομένων
    download_clinvar()
    
    # Εξαγωγή και μετασχηματισμός
    df = extract_and_transform()
    df = calculate_acmg_criteria(df)
    
    # Αποθήκευση σε CSV
    save_to_csv(df)
    
    # Εισαγωγή στη βάση δεδομένων
    insert_to_database(conn, df)
    
    # Κλείσιμο σύνδεσης
    conn.close()
    
    # Καθαρισμός προσωρινών αρχείων
    for f in [GZ_FILE, TSV_FILE]:
        if os.path.exists(f):
            os.remove(f)
    
    print("Ολοκληρώθηκε η επεξεργασία!")

if __name__ == "__main__":
    main()
=======

    # Αξιόπιστοι υποβάλλοντες
    trusted_submitters = {'ClinVar', 'ExpertLab'}
    
    # Θέσεις με γνωστές παθογονικές μεταλλάξεις (π.χ. 41, 504)
    pathogenic_positions = {41, 504}

    # PS1: Ίδιο protein change, διαφορετικό DNA change
    if row['HGVS_p'] in known_pathogenic:
        if row['HGVS_c'] != known_pathogenic[row['HGVS_p']]['dna']:
            criteria.append("PS1")
    
        # PM5
    protein_pos = int(''.join(filter(str.isdigit, row['ProteinChange']))) if pd.notna(row['ProteinChange']) else None
    if protein_pos in pathogenic_positions and row['ProteinChange'] not in known_pathogenic:
        criteria.append('PM5')
    
    # PP5/BP6
    if row['Submitter'] in trusted_submitters:
        if row['ClinicalSignificance'] == 'Pathogenic':
            criteria.append('PP5')
        elif row['ClinicalSignificance'] == 'Benign':
            criteria.append('BP6')
    
    
    return criteria


    
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
    df_brca['rcv_accessions'] = df_brca['RCVaccession'].str.split('|')
    
    return df_brca

#Εισαγωγή στη Βάση
def insert_to_database(conn: psycopg2.extensions.connection, df: pd.DataFrame) -> None:
    """Εισαγωγή δεδομένων στη βάση"""
    with conn.cursor() as cur:
        for _, row in df.iterrows():
            cur.execute("""
            INSERT INTO brca1_variants (
                variation_id, gene_symbol,transcript_id, hgvs_c, hgvs_p,
                clinical_significance, review_status, phenotype_list,
                assembly, chromosome, start_pos, end_pos,
                reference_allele, alternate_allele, acmg_criteria,
                conflicting_interpretations, rcv_accessions, transcript_id
            ) VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
            )
            ON CONFLICT (variation_id) DO UPDATE SET
                gene_symbol = EXCLUDED.gene_symbol,
                hgvs_c = EXCLUDED.hgvs_c,
                hgvs_p = EXCLUDED.hgvs_p,
                morecular_consequence = EXCLUDED.variant_type,
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
                rcv_accessions = EXCLUDED.rcv_accessions,
                transcript_id = EXCLUDED.transcript_id,
                last_updated = CURRENT_TIMESTAMP;
            """, (
                row['VariationID'], row['GeneSymbol'],row['transcript_id'], row['morecular_consequence'], row['HGVS_c'],
                row['HGVS_p'], row['ClinicalSignificance'], row['ReviewStatus'],
                row['PhenotypeList'], row['Assembly'], row['Chromosome'],
                row['Start'], row['Stop'], row['ReferenceAllele'],
                row['AlternateAllele'], json.dumps(row['acmg_criteria']),
                json.dumps(row['conflicting_interpretations']),
                row['rcv_accessions']
            ))
        conn.commit()
    


def categorize_variant_name(name: str) -> dict:
    """ 
    Κατηγοριοποιεί το όνομα μετάλλαξης από το πεδίο 'Name' του ClinVar
    σε μεταλλάξεις DNA (c.), πρωτεΐνης (p.) και άλλους τύπους
    """
    result = {
        'variant_type': None,  # 'DNA', 'Protein', 'Other'
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
            result['variant_type'] = 'DNA'
            result['DNA_variant'] = dna_match.group(1)
        elif dna_star_match:
            result['variant_type'] = 'DNA'
            result['DNA_variant'] = dna_star_match.group(1)
        
        # Έλεγχος για πρωτεϊνικές μεταλλάξεις
        protein_match = protein_pattern.search(name)
        if protein_match and not result['DNA_variant']:  # Προτεραιότητα στις DNA μεταλλάξεις
            result['variant_type'] = 'Protein'
            result['Protein_variant'] = protein_match.group(1)
        
        # Αν δεν βρέθηκε τίποτα από τα παραπάνω
        if not result['variant_type']:
            result['variant_type'] = 'Other'
            result['Other_variant'] = name
    
    return result

def process_clinvar_data(variant_path: str) -> pd.DataFrame:
    """
    Ολοκληρωμένη συνάρτηση επεξεργασίας δεδομένων ClinVar
    με κατηγοριοποίηση βάσει του πεδίου Name
    """
    # Φόρτωση δεδομένων ClinVar
    df = pd.read_csv(variant_path, sep='\t', low_memory=False)
    
    # Φιλτράρισμα για BRCA1 και GRCh38
    df_brca = df[(df['GeneSymbol'] == 'BRCA1') & (df['Assembly'] == 'GRCh38')].copy()
    
    # Κατηγοριοποίηση μεταλλάξεων βάσει του πεδίου Name
    df_brca['VariantName_analysis'] = df_brca['Name'].apply(categorize_variant_name)
    
    # Δημιουργία νέων στηλών
    df_brca['Variant_type'] = df_brca['VariantName_analysis'].apply(lambda x: x['variant_type'])
    df_brca['transcript_id'] = df_brca['Name'].apply(extract_transcript_id)
    df_brca['variant_type'] = df_brca.apply(lambda row: determine_variant_type(row['HGVS_p'], row['HGVS_c']), axis=1)
    df_brca['DNA_variant'] = df_brca['VariantName_analysis'].apply(lambda x: x['DNA_variant'])
    df_brca['Protein_variant'] = df_brca['VariantName_analysis'].apply(lambda x: x['Protein_variant'])
    df_brca['Other_variant'] = df_brca['VariantName_analysis'].apply(lambda x: x['Other_variant'])
    
    # Αφαίρεση της προσωρινής στήλης ανάλυσης
    df_brca.drop(columns=['VariantName_analysis'], inplace=True)
    
    return df_brca

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

# --- Κύρια Λειτουργία ---
def main():
    # Σύνδεση στη βάση
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)
    
    try:
        # Λήψη αρχείων
        variant_gz = "variant_summary.txt.gz"
        submission_gz = "submission_summary.txt.gz"
        
        download_file(CLINVAR_VARIANT_URL, variant_gz)  
        download_file(CLINVAR_SUBMISSION_URL, submission_gz)
        
        # Επεξεργασία δεδομένων
        variant_file = variant_gz.replace(".gz", "")
        submission_file = submission_gz.replace(".gz", "")
        
        df_final = process_clinvar_data(variant_file)
        
        # Εισαγωγή στη βάση
        insert_to_database(conn, df_final)
        
        print("Επεξεργασία ολοκληρώθηκε επιτυχώς!")
        
    except Exception as e:
        print(f"Σφάλμα: {e}")
    finally:
        # Καθαρισμός και κλείσιμο σύνδεσης
        for f in [variant_gz, submission_gz, variant_file, submission_file]:
            if os.path.exists(f):
                os.remove(f)
        conn.close()

if __name__ == "__main__":
    main()

   
>>>>>>> 373c3d0 (changes)
