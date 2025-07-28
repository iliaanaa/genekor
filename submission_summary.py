import pandas as pd
import psycopg2
import os
import requests
import shutil

from typing import List

# URL για το submission_summary
CLINVAR_SUBMISSION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
TRUSTED_SUBMITTERS = {'ClinVar', 'ExpertLab'}

# --- Κατέβασμα και αποσυμπίεση ---
def download_file(url: str, dest: str):
    """Κατεβάζει gzip αρχείο"""
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(dest, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    else:
        raise Exception(f"Λάθος στο κατέβασμα: {url} (status {response.status_code})")

def extract_gzip(source: str, dest: str):
    """Αποσυμπιέζει gzip αρχείο"""
    import gzip
    with gzip.open(source, 'rb') as f_in, open(dest, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# --- Υπολογισμός ACMG ---
def apply_acmg(row: pd.Series) -> List[str]:
    """Υπολογισμός κριτηρίων ACMG για μια μετάλλαξη"""
    criteria = []
    if row['Submitter'] in TRUSTED_SUBMITTERS:
        if row['ClinicalSignificance'] == 'Pathogenic':
            criteria.append('PP5')
        elif row['ClinicalSignificance'] == 'Benign':
            criteria.append('BP6')
    return criteria

# --- Επεξεργασία Δεδομένων ---
def process_submission_data(filename: str) -> pd.DataFrame:
    """Φόρτωση και επεξεργασία του submission_summary"""
    df = pd.read_csv(filename, sep='\t', low_memory=False)

    # Υπολογισμός ACMG
    df['ACMG_Criteria'] = df.apply(apply_acmg, axis=1)
    return df

# --- Δημιουργία Πίνακα στη Βάση ---
def create_submissions_table(conn):
    """Δημιουργεί τον πίνακα submissions αν δεν υπάρχει"""
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS submissions (
            variation_id INT PRIMARY KEY,
            clinical_significance TEXT,
            datelastevaluated DATE,
            description TEXT,
            submittedphenotypeinfo TEXT,
            reportedphenotypeinfo TEXT,
            reviewstatus TEXT,
            collectionmethod TEXT,
            origincounts TEXT,
            submitter TEXT,
            scv TEXT,
            submittedgenesymbol TEXT,
            explanationofinterpretation TEXT,
            somaticclinicalimpact TEXT,
            oncogenicity TEXT,
            acmg_criteria TEXT
        );
    """)
    conn.commit()
    cur.close()

# --- Εισαγωγή στη Βάση ---
def insert_submission_to_db(conn, df: pd.DataFrame):
    """Εισαγωγή του submission_summary στη βάση"""
    cur = conn.cursor()
    for _, row in df.iterrows():
        cur.execute("""
            INSERT INTO submissions (
                variation_id, clinical_significance, datelastevaluated, description,
                submittedphenotypeinfo, reportedphenotypeinfo, reviewstatus,
                collectionmethod, origincounts, submitter, scv,
                submittedgenesymbol, explanationofinterpretation,
                somaticclinicalimpact, oncogenicity, acmg_criteria
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (variation_id) DO NOTHING
        """, (
            row['VariationID'], row['ClinicalSignificance'], row['DateLastEvaluated'],
            row['Description'], row['SubmittedPhenotypeInfo'], row['ReportedPhenotypeInfo'],
            row['ReviewStatus'], row['CollectionMethod'], row['OriginCounts'],
            row['Submitter'], row['SCV'], row['SubmittedGeneSymbol'],
            row['ExplanationOfInterpretation'], row['SomaticClinicalImpact'],
            row['Oncogenicity'], ','.join(row['ACMG_Criteria'])
        ))
    conn.commit()
    cur.close()

# --- Κύρια Συνάρτηση ---
def main():
    from config import DB_CONFIG  # Βεβαιώσου ότι έχεις DB_CONFIG στο config.py
    conn = psycopg2.connect(**DB_CONFIG)
    create_submissions_table(conn)

    try:
        # Λήψη και αποσυμπίεση
        submission_gz = "submission_summary.txt.gz"
        submission_txt = "submission_summary.txt"
        download_file(CLINVAR_SUBMISSION_URL, submission_gz)
        extract_gzip(submission_gz, submission_txt)

        # Επεξεργασία
        df_submission = process_submission_data(submission_txt)

        # Εισαγωγή DB
        insert_submission_to_db(conn, df_submission)
        print(" Το submission_summary επεξεργάστηκε και φορτώθηκε επιτυχώς!")
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    main()
