import pandas as pd
import psycopg2
import os
import gzip
import shutil
import urllib.request
import json
from datetime import datetime

#--- Ρυθμίσεις ---
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GZ_FILE = "variant_summary.txt.gz"
TSV_FILE = "variant_summary.txt"
GENE_FILTER = "BRCA1"
OUTPUT_CSV = "brca_variants.csv"

#--- Ρυθμίσεις Βάσης Δεδομένων ---
DB_NAME = "clinvar_db"
DB_USER = "postgres"
DB_PASS = "your_password"
DB_HOST = "localhost"
DB_PORT = 5432

#--- Σύνδεση στη Βάση ---
def connect_db():
return psycopg2.connect(
dbname=DB_NAME,
user=DB_USER,
password=DB_PASS,
host=DB_HOST,
port=DB_PORT
)

#--- Δημιουργία Πινάκων ---
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
    );
    """)
    conn.commit()
#--- Λήψη έκδοσης ClinVar ---
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

#--- Κατέβασμα και αποσυμπίεση δεδομένων ---
def download_clinvar():
print("Κατέβασμα δεδομένων ClinVar...")
urllib.request.urlretrieve(CLINVAR_URL, GZ_FILE)
print("Αποσυμπίεση...")
with gzip.open(GZ_FILE, 'rb') as f_in, open(TSV_FILE, 'wb') as f_out:
shutil.copyfileobj(f_in, f_out)

#--- Εξαγωγή και μετασχηματισμός δεδομένων BRCA1 ---
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
#--- Υπολογισμός ACMG κριτηρίων ---
def calculate_acmg_criteria(df):
print("Υπολογισμός ACMG κριτηρίων...")

# Λεξικό με γνωστές pathogenic μεταλλάξεις
known_pathogenic = {
    'p.Arg504Gly': {'dna': 'c.1510A>T', 'significance': 'Pathogenic'},
    'p.Trp41*': {'dna': 'c.123G>A', 'significance': 'Pathogenic'}
}

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
#--- Αποθήκευση σε CSV ---
def save_to_csv(df, filename=OUTPUT_CSV):
print(f"Αποθήκευση σε {filename}...")

# Μετατροπή των λιστών/λεξικών σε JSON strings για το CSV
df_csv = df.copy()
df_csv['acmg_criteria'] = df_csv['acmg_criteria'].apply(lambda x: json.loads(x) if x else [])

df_csv.to_csv(filename, index=False, encoding='utf-8')
#--- Εισαγωγή στη βάση δεδομένων ---
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

#--- Κύρια λειτουργία ---
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
if name == "main":
main()