import pandas as pd
import psycopg2
import os
import gzip
import shutil
import urllib.request
import json
from datetime import datetime
from typing import Dict, List, Optional

# --- Ρυθμίσεις ---
CLINVAR_VARIANT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_SUBMISSION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz" #μορφή gzip
GENE_FILTER = "BRCA1"
ASSEMBLY_FILTER = "GRCh38"  #Επιλέγει την έκδοση GRCh38 του γονιδιώματος

# --- Ρυθμίσεις Βάσης Δεδομένων ---
DB_CONFIG = {
    "dbname": "clinvar_db",
    "user": "postgres",
    "password": "your_password",
    "host": "localhost",
    "port": 5432
}

# --- Βοηθητικές Συναρτήσεις ---
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
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            last_evaluated DATE
        );
        """)
        conn.commit()

# --- Λογική ACMG Criteria ---
def apply_acmg_criteria(row: pd.Series) -> List[str]:
    """Υπολογισμός κριτηρίων ACMG για μια μετάλλαξη"""
    criteria = []
    
    # Γνωστές pathogenic μεταλλάξεις (προσαρμόστε ανάλογα)
    known_pathogenic = {
        'p.Arg504Gly': {'dna': 'c.1510A>T', 'significance': 'Pathogenic'},
        'p.Trp41*': {'dna': 'c.123G>A', 'significance': 'Pathogenic'}
    }

    # Αξιόπιστοι υποβάλλοντες
    trusted_submitters = {'ClinVar', 'ExpertLab'}
    
    # Θέσεις με γνωστές παθογονικές μεταλλάξεις (π.χ. 41, 504)
    pathogenic_positions = {41, 504}

    # PS1: Ίδιο protein change, διαφορετικό DNA change
    if row['HGVS_p'] in known_pathogenic:
        if row['HGVS_c'] != known_pathogenic[row['HGVS_p']]['dna']:
            criteria.append("PS1")
    
    # Προσθέστε εδώ άλλα κριτήρια (PM5, PP5, κλπ)
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

# --- Επεξεργασία Δεδομένων ---
def process_clinvar_data(variant_path: str, submission_path: str) -> pd.DataFrame:
    #Φορτώνει το variant_summary και φιλτράρει για BRCA1 στο assembly GRCh38.

#Επιλέγει μόνο κάποιες στήλες.

#Φορτώνει το submission_summary.

#Κάνει merge με βάση το VariationID.

#Τυποποιεί το πεδίο κλινικής σημασίας (ClinicalSignificance).

#Υπολογίζει τα ACMG κριτήρια για κάθε μετάλλαξη.

#Ομαδοποιεί τα conflicting interpretations (διαφωνίες κλινικών αξιολογήσεων).

#Δημιουργεί το τελικό DataFrame με όλα τα πεδία που θα αποθηκευτούν.

    """Φόρτωση και επεξεργασία δεδομένων ClinVar"""
    # Φόρτωση variant data
    df_variants = pd.read_csv(variant_path, sep='\t', low_memory=False)
    
    # Φιλτράρισμα για BRCA1 και συγκεκριμένο assembly
    df_brca = df_variants[
        (df_variants['GeneSymbol'] == GENE_FILTER) & 
        (df_variants['Assembly'] == ASSEMBLY_FILTER)
    ].copy()
    
    # Επιλογή στηλών
    columns_to_keep = [
        'VariationID', 'GeneSymbol', 'HGVS_c', 'HGVS_p', 
        'ClinicalSignificance', 'ReviewStatus', 'PhenotypeList',
        'Submitter', 'Assembly', 'Chromosome', 'Start', 'Stop',
        'ReferenceAllele', 'AlternateAllele', 'RCVaccession'
    ]
    df_brca = df_brca[columns_to_keep]
    
    # Φόρτωση submission data
    df_submissions = pd.read_csv(submission_path, sep='\t')
    
    # Συγχώνευση και επεξεργασία
    df_merged = pd.merge(
        df_brca,
        df_submissions,
        on='VariationID',
        how='left'
    )
    
    # Τυποποίηση κλινικής σημασίας
    significance_map = {
        'Pathogenic': 'Pathogenic',
        'Likely pathogenic': 'Likely_pathogenic',
        'Benign': 'Benign',
        'Uncertain significance': 'VUS'
    }
    df_merged['ClinicalSignificance'] = df_merged['ClinicalSignificance_x'].map(significance_map)
    
    # Υπολογισμός ACMG criteria
    df_merged['acmg_criteria'] = df_merged.apply(apply_acmg_criteria, axis=1)
    
    # Ομαδοποίηση conflicting interpretations
    conflicting = df_merged.groupby('VariationID').apply(
        lambda x: x[['Submitter', 'ClinicalSignificance_y']].rename(
            columns={'ClinicalSignificance_y': 'classification'}
        ).to_dict('records')
    )
    
    # Δημιουργία τελικού DataFrame
    df_final = df_brca.drop_duplicates(subset=['VariationID']).set_index('VariationID')
    df_final['acmg_criteria'] = df_merged.groupby('VariationID')['acmg_criteria'].first()
    df_final['conflicting_interpretations'] = conflicting
    df_final['rcv_accessions'] = df_final['RCVaccession'].str.split('|')
    
    return df_final.reset_index()

# --- Εισαγωγή στη Βάση ---
def insert_to_database(conn: psycopg2.extensions.connection, df: pd.DataFrame) -> None:
    """Εισαγωγή δεδομένων στη βάση"""
    with conn.cursor() as cur:
        for _, row in df.iterrows():
            cur.execute("""
            INSERT INTO brca1_variants (
                variation_id, gene_symbol, hgvs_c, hgvs_p,
                clinical_significance, review_status, phenotype_list,
                assembly, chromosome, start_pos, end_pos,
                reference_allele, alternate_allele, acmg_criteria,
                conflicting_interpretations, rcv_accessions
            ) VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
            )
            ON CONFLICT (variation_id) DO UPDATE SET
                gene_symbol = EXCLUDED.gene_symbol,
                hgvs_c = EXCLUDED.hgvs_c,
                hgvs_p = EXCLUDED.hgvs_p,
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
                last_updated = CURRENT_TIMESTAMP;
            """, (
                row['VariationID'], row['GeneSymbol'], row['HGVS_c'],
                row['HGVS_p'], row['ClinicalSignificance'], row['ReviewStatus'],
                row['PhenotypeList'], row['Assembly'], row['Chromosome'],
                row['Start'], row['Stop'], row['ReferenceAllele'],
                row['AlternateAllele'], json.dumps(row['acmg_criteria']),
                json.dumps(row['conflicting_interpretations']),
                row['rcv_accessions']
            ))
        conn.commit()

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
        
        df_final = process_clinvar_data(variant_file, submission_file)
        
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

    #to api sto telos tha apantaei me enq json to eides rhs metallahs kai to kritrio. ayto to json pairnei plhrofoeis mono apo to csv h kai apo thn vash?