# ΡΥΘΜΙΣΕΙΣ
CLINVAR_VARIANT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_SUBMISSION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
GENE_FILTER = "BRCA1"
ASSEMBLY_FILTER = "GRCh38"

Ορίζει:

URLs για κατέβασμα δεδομένων ClinVar (variants και submissions)

Γονίδιο στόχος: BRCA1

Έκδοση γονιδιώματος: GRCh38

# ΣΥΝΔΕΣΗ ΜΕ ΒΑΣΗ ΔΕΔΟΜΕΝΩΝ
DB_CONFIG = {
    "dbname": "clinvar_db",
    "user": "postgres",
    "password": "your_password",
    "host": "localhost",
    "port": 5432
}
Ορίζει τα credentials για σύνδεση σε PostgreSQL βάση.

# download_file(url, output_path)
Κατεβάζει αρχείο και το αποσυμπιέζει αν είναι .gz.

# create_tables(conn)
Δημιουργεί τον πίνακα brca1_variants με πεδία όπως:

variation_id, hgvs_c, hgvs_p, clinical_significance, review_status, κλπ

acmg_criteria και conflicting_interpretations σε μορφή JSONB

rcv_accessions ως array

last_updated, last_evaluated για tracking

# apply_acmg_criteria(row)
Υπολογίζει κριτήρια ACMG, όπως:

PS1: ίδιο protein change, διαφορετικό DNA change

PM5: ίδια θέση με γνωστή παθογόνο μετάλλαξη, αλλά διαφορετικό aa

PP5/BP6: αν υπάρχει έγκυρος submitter

# Χρησιμοποιεί μια dictionary known_pathogenic για να εντοπίσει “γνωστές” μεταλλάξεις.

# categorize_variant_name(name)
Αναλύει το πεδίο Name και εξάγει:

DNA_variant (π.χ. c.1510A>T)

Protein_variant (π.χ. p.Arg504Gly)

Εναλλακτικά, τα κρατά ως Other_variant

Χρησιμοποιεί κανονικές εκφράσεις (regex) για c. και p. μεταλλάξεις (συμπεριλαμβανομένων περιπτώσεων με * όπως c.*103_*106del).

# process_clinvar_data(variant_path)
Φορτώνει το .txt (μετά από αποσυμπίεση).

Φιλτράρει για BRCA1 + GRCh38.

Αναλύει τα πεδία Name και δημιουργεί στήλες για τύπο μετάλλαξης (Variant_type, DNA_variant, Protein_variant, κλπ).

Υπολογίζει:

conflict_score: βαθμός διαφωνίας ερμηνειών

has_conflicts: boolean αν υπάρχουν ουσιώδεις διαφωνίες

acmg_criteria: λίστα με ACMG tags

rcv_accessions: array από accession strings

# insert_to_database(conn, df)
Κάνει INSERT ... ON CONFLICT UPDATE ώστε να ενημερώνει εγγραφές αν υπάρχει ήδη το variation_id.

# main()
Ο κορμός:

Συνδέεται με βάση.

Κατεβάζει και αποσυμπιέζει αρχεία.

Εκτελεί process_clinvar_data().

Εισάγει δεδομένα με insert_to_database().

Διαγράφει προσωρινά αρχεία.

Κλείνει σύνδεση.