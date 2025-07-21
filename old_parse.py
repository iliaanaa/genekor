import pandas as pd
import psycopg2
import os
import re
import gzip
import subprocess
import shutil
import urllib.request
import json
import requests
import gzip    
import subprocess
import shutil
from datetime import datetime
from typing import Dict, List, Optional

# --- Ρυθμίσεις ---
CLINVAR_VARIANT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_SUBMISSION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz" #μορφή gzip
GENE_FILTER = "KLHL10"
ASSEMBLY_FILTER = "GRCh38"  #Επιλέγει την έκδοση GRCh38 του γονιδιώματος
CLINVAR_SUBMISSION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
TRUSTED_SUBMITTERS = {'ClinVar', 'ENIGMA'}

#Ρυθμίσεις Βάσης Δεδομένων
DB_CONFIG = {
    "dbname": "clinvar_db",
    "user": "ilianam",
    "password": "wasd1029!",
    "host": "localhost",
    "port": 5432
}

#Βοηθητικές Συναρτήσεις
'''
def download_file(url: str, output_path: str) -> None:
    """Κατέβασμα και αποσυμπίεση αρχείου"""
    print(f"Κατέβασμα {url}...")
    urllib.request.urlretrieve(url, output_path)
    
    if output_path.endswith(".gz"):
        print(f"Αποσυμπίεση {output_path}...")
        with gzip.open(output_path, 'rb') as f_in:
            with open(output_path.replace(".gz", ""), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)'''

'''
def download_and_extract_clinvar_files() -> None:
    """Κατεβάζει και αποσυμπιέζει τα variant_summary και submission_summary της ClinVar"""

    files = {
        "variant_summary.txt.gz": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
        "submission_summary.txt.gz": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
    }

    for gz_filename, url in files.items():
        txt_filename = gz_filename.replace(".gz", "")
        print(f"\nΚατέβασμα: {url}")

        # Κατέβασμα .gz
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(gz_filename, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"✓ Το αρχείο αποθηκεύτηκε ως {gz_filename}")
        else:
            raise Exception(f"Σφάλμα στο κατέβασμα του {url} (status {response.status_code})")

        # Αποσυμπίεση
        print(f"Αποσυμπίεση: {gz_filename} → {txt_filename}")
        with gzip.open(gz_filename, 'rb') as f_in, open(txt_filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"✓ Δημιουργήθηκε: {txt_filename}")
'''
def get_header_from_commented_tsv(gz_path: str) -> list:
    """
    Εντοπίζει τη γραμμή header που ξεκινάει με # και επιστρέφει τα ονόματα στηλών ως λίστα.
    """
    with gzip.open(gz_path, 'rt') as f:
        for line in f:
            if line.startswith('#') and 'VariationID' in line:
                return line.lstrip('#').strip().split('\t')
    raise ValueError("Δεν βρέθηκε γραμμή κεφαλίδας με 'VariationID'")

def merge_variant_submission(variant_path: str, submission_path: str) -> pd.DataFrame:
        """
        Συνδυάζει variant_summary και μόνο τους Submitter από submission_summary με βάση το VariationID.
        """

        print("Φόρτωση variant_summary...")
        df_variant = pd.read_csv(variant_path, sep='\t', low_memory=False)

        # Φιλτράρισμα variant_summary
        df_variant = df_variant[
            (df_variant['GeneSymbol'] == GENE_FILTER) &
            (df_variant['Assembly'] == ASSEMBLY_FILTER)
        ].copy()

        print(f"Μετά το φιλτράρισμα: {len(df_variant)} εγγραφές")

        # Εξαγωγή unique VariationIDs
        variation_ids = set(df_variant['VariationID'].dropna().unique())

        print("Φόρτωση και φιλτράρισμα submission_summary...")
        submitters_dict = {}

        print("Φόρτωση και φιλτράρισμα submission_summary...")

        # ➤ Πάρε το σωστό header από το αρχείο gzip
        header_line = get_header_from_commented_tsv(submission_path)
        print("🔎 Header columns:", header_line)


        # ➤ Ανάγνωση με το σωστό header
        submitters_dict = {}
        for chunk in pd.read_csv(
            submission_path,
            sep='\t',
            comment='#',
            names=header_line,
            compression='gzip',
            usecols=['VariationID', 'Submitter'],
            chunksize=500_000
        ):
            filtered = chunk[chunk['VariationID'].isin(variation_ids)]
            for _, row in filtered.iterrows():
                vid = row['VariationID']
                sub = row['Submitter']
                if pd.notna(sub):
                    submitters_dict.setdefault(vid, set()).add(sub)

        # Μετατροπή σε DataFrame
        df_submitters = pd.DataFrame([
            {'VariationID': vid, 'Submitter': sorted(list(subs))}
            for vid, subs in submitters_dict.items()
        ])

        print(f"✓ Βρέθηκαν submitters για {len(df_submitters)} VariationIDs")

        # Συγχώνευση
        merged_df = pd.merge(
            df_variant,
            df_submitters,
            how='left',
            on='VariationID'
        )

        # Αν θέλεις pipe-separated string αντί για λίστα:
        merged_df['Submitter'] = merged_df['Submitter'].apply(
            lambda x: '|'.join(x) if isinstance(x, list) else None
        )

        print(f"✓ Συγχωνεύθηκαν: {merged_df['Submitter'].notna().sum()}/{len(merged_df)} με Submitter")

        return merged_df



#grep h zcat gia na mh skaei
def create_tables(conn: psycopg2.extensions.connection) -> None:
    """Δημιουργία πινάκων στη βάση"""
    with conn.cursor() as cur:
        # Κύριος πίνακας μεταλλάξεων
        cur.execute("""
        CREATE TABLE IF NOT EXISTS brca1_variants (
            variation_id BIGINT PRIMARY KEY,
            gene_symbol TEXT NOT NULL,
            transcript_id TEXT,
            variant_type TEXT,
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
        );
        """)
        conn.commit()

#ACMG Criteria
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
    
        # PM5
    protein_pos = int(''.join(filter(str.isdigit, row['HGVS_p']))) if pd.notna(row['HGVS_p']) else None
    if protein_pos in pathogenic_positions and row['HGVS_p'] not in known_pathogenic:
        criteria.append('PM5')
    
    # PP5/BP6

    submitter = row.get('Submitter', None)
    if submitter in trusted_submitters:
        if row['ClinicalSignificance'] == 'Pathogenic':
            criteria.append('PP5')
        elif row['ClinicalSignificance'] == 'Benign':
            criteria.append('BP6')

        '''
    if row['Submitter'] in trusted_submitters:
        if row['ClinicalSignificance'] == 'Pathogenic':
            criteria.append('PP5')
        elif row['ClinicalSignificance'] == 'Benign':
            criteria.append('BP6')
'''
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
                variation_id, gene_symbol, transcript_id, hgvs_c, hgvs_p, variant_type,
                clinical_significance, review_status, phenotype_list,
                assembly, chromosome, start_pos, end_pos,
                reference_allele, alternate_allele, acmg_criteria,
                conflicting_interpretations, rcv_accessions
            ) VALUES (
                %s, %s, %s, %s, %s, %s,
                %s, %s, %s,
                %s, %s, %s, %s,
                %s, %s, %s,
                %s, %s
            )
            ON CONFLICT (variation_id) DO UPDATE SET
                gene_symbol = EXCLUDED.gene_symbol,
                transcript_id = EXCLUDED.transcript_id,
                hgvs_c = EXCLUDED.hgvs_c,
                hgvs_p = EXCLUDED.hgvs_p,
                variant_type = EXCLUDED.variant_type,
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
                row['VariationID'], row['GeneSymbol'], row['transcript_id'], row['HGVS_c'], row['HGVS_p'], row['variant_type'],
                row['ClinicalSignificance'], row['ReviewStatus'], row['PhenotypeList'],
                row['Assembly'], row['Chromosome'], row['Start'], row['Stop'],
                row['ReferenceAllele'], row['AlternateAllele'], json.dumps(row['acmg_criteria']),
                json.dumps(row['conflicting_interpretations']), row['rcv_accessions']
            ))
        conn.commit()

        
    '''
def extract_HGVS(name: str) -> dict:
    """
    Εξάγει HGVS.c και HGVS.p από το πεδίο name χρησιμοποιώντας τα συγκεκριμένα regex patterns
    """
    result = {
        'HGVS_c': None,
        'HGVS_p': None,
        'variant_type': None,
        'Other_variant': None
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

        # Αν δεν βρέθηκε τίποτα από τα παραπάνω
        if not result['variant_type']:
            result['variant_type'] = 'Other'
            result['Other_variant'] = name

    return result
'''
def extract_HGVS(name: str) -> dict:
    """
    Εξάγει HGVS.c, HGVS.p και variant_type από το πεδίο name χρησιμοποιώντας regex patterns.
    """
    result = {
        'HGVS_c': None,
        'HGVS_p': None,
        'variant_type': None,
        'Other_variant': None
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
            result['variant_type'] = 'DNA'
        elif dna_star_match:
            result['HGVS_c'] = dna_star_match.group(1)
            result['variant_type'] = 'DNA'
        
        # Αναζήτηση για πρωτεϊνικές μεταλλάξεις
        protein_match = protein_pattern.search(name)
        if protein_match:
            result['HGVS_p'] = protein_match.group(1)
            # Αν δεν έχει ήδη variant_type DNA, βάλε protein
            if result['variant_type'] is None:
                result['variant_type'] = 'Protein'
        
        # Αν δεν βρέθηκε ούτε DNA ούτε Protein
        if result['variant_type'] is None:
            result['variant_type'] = 'Other'
            result['Other_variant'] = name

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

'''

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
'''
'''
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
    '''

#######################################################
'''
def process_clinvar_data(variant_gz_path: str) -> pd.DataFrame:
    """
    Επεξεργασία δεδομένων ClinVar ΑΜΕΣΑ από το .gz αρχείο,
    χωρίς πλήρη αποσυμπίεση. Χρησιμοποιεί zcat και grep.
    """

    print("Φόρτωση όλων των δεδομένων από το gzip...")
    with gzip.open(variant_gz_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', low_memory=False)

    
    # 1. Φιλτράρισμα για γονίδιο KLHL10 και GRCh38 με grep
    print("Φιλτράρισμα δεδομένων με grep...")
    grep_cmd = f"zcat {variant_gz_path} | grep -E 'KLHL10.*GRCh38' > filtered_variants.tsv"
    subprocess.run(grep_cmd, shell=True, check=True)
    
    # 2. Φόρτωση ΜΟΝΟ των φιλτραρισμένων δεδομένων
    print("Φόρτωση φιλτραρισμένων δεδομένων...")
    df = pd.read_csv("filtered_variants.tsv", sep='\t', low_memory=False)
    
    # 3. Κατηγοριοποίηση μεταλλάξεων (όπως πριν)
    df['VariantName_analysis'] = df['Name'].apply(categorize_variant_name)
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
    
    # 4. Διαγραφή προσωρινού αρχείου
    os.remove("filtered_variants.tsv")
    
    return df

###########################################################################
    '''


def process_clinvar_data(variant_gz_path: str) -> pd.DataFrame:
    """
    Επεξεργασία δεδομένων ClinVar από .gz με zcat και grep (για χαμηλή χρήση μνήμης).
    Διατηρεί τη γραμμή τίτλων (header).
    """


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
    '''
    df['VariantName_analysis'] = df['Name'].apply(extract_HGVS)
    df[['HGVS_c', 'HGVS_p']] = df['VariantName_analysis'].apply(pd.Series)
    if 'HGVS_p' not in df.columns or 'HGVS_c' not in df.columns:
        raise ValueError("Λείπουν οι στήλες HGVS_p ή HGVS_c πριν το determine_variant_type.")

    df['transcript_id'] = df['Name'].apply(extract_transcript_id)
    df['variant_type'] = df.apply(lambda row: determine_variant_type(row['HGVS_p'], row['HGVS_c']), axis=1)
    df['HGVS_c'] = df['VariantName_analysis'].apply(lambda x: x['HGVS_c'])
    df['PHGVS_p'] = df['VariantName_analysis'].apply(lambda x: x['HGVS_p'])
    df['Other_variant'] = df['VariantName_analysis'].apply(lambda x: x['Other_variant'])

    df.drop(columns=['VariantName_analysis'], inplace=True)

    # Διαγραφή προσωρινού αρχείου
    os.remove("filtered_variants.tsv")
'''
    # Εξαγωγή HGVS_c, HGVS_p, variant_type, Other_variant σε μία γραμμή
    df = pd.concat([df, df['Name'].apply(extract_HGVS).apply(pd.Series)], axis=1)

    # Εξαγωγή transcript ID
    df['transcript_id'] = df['Name'].apply(extract_transcript_id)

    # Διαγραφή προσωρινού αρχείου
    os.remove("filtered_variants.tsv")
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

'''
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

        '''


def filter_submission_by_variation_ids(submission_gz_path: str, variation_ids: List[str]) -> pd.DataFrame:
    """
    Φιλτράρει το submission_summary.gz ώστε να κρατήσει μόνο τις γραμμές με τα συγκεκριμένα VariationID.
    """
    temp_filename = "filtered_submission.tsv"
    
    # Γράφουμε το header
    with open(temp_filename, "w") as f_out:
        header_cmd = f"zcat {submission_gz_path} | head -n 1"
        subprocess.run(header_cmd, shell=True, stdout=f_out)

    # Δημιουργούμε regex pattern από όλα τα IDs (π.χ. 12345|67890|...)
    id_pattern = "|".join(map(str, variation_ids))
    
    # Προσθέτουμε γραμμές με αυτά τα IDs
    grep_cmd = f"zcat {submission_gz_path} | grep -E '{id_pattern}' >> {temp_filename}"
    subprocess.run(grep_cmd, shell=True, check=True)
    
    # Φορτώνουμε σε DataFrame
    df_submissions = pd.read_csv(temp_filename, sep='\t', low_memory=False)
    
    # Διαγραφή προσωρινού αρχείου
    os.remove(temp_filename)
    
    return df_submissions

   
'''
def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)
    
    try:
        print("Ξεκίνημα script...")  # Για επιβεβαίωση ότι τρέχει
        # Λήψη αρχείων (ΜΟΝΟ τα .gz, χωρίς αποσυμπίεση)
        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)
        submission_gz = "submission_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_SUBMISSION_URL, submission_gz)

        
        # Επεξεργασία ΑΜΕΣΑ από το .gz
        df_final = process_clinvar_data(variant_gz)
        # Εξαγωγή variation IDs για να φιλτράρουμε το submission_summary
        variation_ids = df_final['VariationID'].astype(str).tolist()
        df_sumbissions = filter_submission_by_variation_ids(submission_gz, variation_ids)
        df_final['acmg_criteria'] = df_final.apply(apply_acmg_criteria, axis=1)
        df_final['conflicting_interpretations'] = [{} for _ in range(len(df_final))]
        df_final['rcv_accessions'] = df_final['RCVaccession'].fillna('').apply(lambda x: x.split('|') if x else [])
        df_final['Submitter'] = None

        # Εισαγωγή στη βάση
        insert_to_database(conn, df_final)
        
    except Exception as e:
        print(f"Σφάλμα: {e}")
    finally:
        if os.path.exists(variant_gz):
            os.remove(variant_gz)
        if os.path.exists(submission_gz):
            os.remove(submission_gz)
        conn.close()
'''
def extract_submitters(submission_path: str, variation_ids: set) -> pd.DataFrame:
    """
    Διαβάζει το submission_summary σε chunks και επιστρέφει μοναδικούς Submitters ανά VariationID.
    """
    submitters = {}

   
    for chunk in pd.read_csv(submission_path, sep='\t', usecols=['VariationID', 'Submitter'], chunksize=500_000):
        filtered = chunk[chunk['VariationID'].isin(variation_ids)]
        for _, row in filtered.iterrows():
            vid = row['VariationID']
            sub = row['Submitter']
            if pd.notna(sub):
                submitters.setdefault(vid, set()).add(sub)

    # Μετατροπή σε DataFrame
    df = pd.DataFrame([
        {'VariationID': vid, 'Submitter': sorted(list(subs))}
        for vid, subs in submitters.items()
    ])

    return df





'''
def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)

    try:
        print("Ξεκίνημα script...")

        # Λήψη .gz αρχείων
        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)
        submission_gz = "submission_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_SUBMISSION_URL, submission_gz)

        # Φόρτωση και φιλτράρισμα variant_summary
        df_final = process_clinvar_data(variant_gz)

        # Εξαγωγή μόνο των VariationID
        variation_ids = set(df_final['VariationID'].unique())

        # Εξαγωγή Submitters για κάθε VariationID
        df_submitters = extract_submitters(submission_gz, variation_ids)

        # Συγχώνευση στο τελικό dataframe
        df_final = df_final.merge(df_submitters, how='left', on='VariationID')

        # Αν θες pipe-separated string αντί για λίστα:
        df_final['Submitter'] = df_final['Submitter'].apply(
            lambda x: '|'.join(x) if isinstance(x, list) else None
        )

        # ACMG και επιπλέον στήλες
        df_final['acmg_criteria'] = df_final.apply(apply_acmg_criteria, axis=1)
        df_final['conflicting_interpretations'] = [{} for _ in range(len(df_final))]
        df_final['rcv_accessions'] = df_final['RCVaccession'].fillna('').apply(lambda x: x.split('|') if x else [])

        # Εισαγωγή στη βάση
        insert_to_database(conn, df_final)

    except Exception as e:
        print(f"Σφάλμα: {e}")

    finally:
        # Καθαρισμός αρχείων
        if os.path.exists(variant_gz):
            os.remove(variant_gz)
        if os.path.exists(submission_gz):
            os.remove(submission_gz)
        conn.close()
'''
def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)

    try:
        print("Ξεκίνημα script...")

        # Κατέβασμα αρχείων
        variant_gz = "variant_summary.txt.gz"
        submission_gz = "submission_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)
        urllib.request.urlretrieve(CLINVAR_SUBMISSION_URL, submission_gz)

        # ➤ Συγχώνευση variant + submission με βάση VariationID
        df_final = merge_variant_submission(variant_gz, submission_gz)

        # ➤ ACMG criteria + επιπλέον στήλες
        df_final['acmg_criteria'] = df_final.apply(apply_acmg_criteria, axis=1)
        df_final['conflicting_interpretations'] = [{} for _ in range(len(df_final))]
        df_final['rcv_accessions'] = df_final['RCVaccession'].fillna('').apply(
            lambda x: x.split('|') if x else []
        )

        # ➤ Εισαγωγή στη βάση
        insert_to_database(conn, df_final)

    except Exception as e:
        print(f"Σφάλμα: {e}")

    finally:
        # ➤ Καθαρισμός προσωρινών αρχείων
        for path in [variant_gz, submission_gz]:
            if os.path.exists(path):
                os.remove(path)
        conn.close()


if __name__ == "__main__":
    main()