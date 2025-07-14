
SYMPERASMA
columns_to_keep = [
    'VariationID', 'GeneSymbol', 'HGVS_c', 'HGVS_p', 
    'ClinicalSignificance', 'ReviewStatus', 'PhenotypeList',
    'Submitter', 'Assembly', 'Chromosome', 'Start', 'Stop',
    'ReferenceAllele', 'AlternateAllele'
]
# filtrarisma
df_variants = pd.read_csv('variant_summary.txt', sep='\t', low_memory=False)
df_brca = df_variants[df_variants['GeneSymbol'] == 'BRCA1']

# sygxwneysh me to submission_summary
df_submissions = pd.read_csv('submission_summary.txt', sep='\t')
df_merged = pd.merge(
    df_brca, 
    df_submissions, 
    on='VariationID', 
    how='left'
) #plhrofories apo alla submissions gia plhresterh eikona   


# KATHARISMOS KAI METASXHMATISMOS DATA
# typopoihsh klinikhs shmasias
significance_map = {
    'Pathogenic': 'Pathogenic',
    'Likely pathogenic': 'Likely_pathogenic',
    'Benign': 'Benign',
    'Uncertain significance': 'VUS'
}
df_merged['ClinicalSignificance'] = df_merged['ClinicalSignificance'].map(significance_map)

# ejagwgh acmg krihthria(logikh idia me krhthria.py)
df_merged['acmg_criteria'] = df_merged.apply(apply_acmg_criteria, axis=1)

# sympiesh plhroforiwn ypovolhs
gia metallajeis me polles ypovoles omadopoiw
conflicting = df_merged.groupby('VariationID').apply(
    lambda x: x[['Submitter', 'ClinicalSignificance']].to_dict('records')
)
df_final['conflicting_interpretations'] = df_final['VariationID'].map(conflicting)

# TELIKH DOMH CSV
variant_id,gene,dna_change,protein_change,clinical_significance,review_status,acmg_criteria,condition,submitter,conflicting_interpretations
1831,BRCA1,c.1510A>G,p.Arg504Gly,Pathogenic,Expert panel,"['PS1', 'PP5']","Breast-ovarian cancer","ENIGMA",null
1832,BRCA1,c.123G>A,p.Trp41*,Likely_pathogenic,Conflicting,"['PM5']","Ovarian cancer","LabX","[{'submitter': 'LabA', 'classification': 'Pathogenic'}, {'submitter': 'LabB', 'classification': 'Benign'}]"


Σημαντικές Σημειώσεις
Μετατροπή HGVS:

HGVS_c (DNA change): c.1510A>G

HGVS_p (Protein change): p.Arg504Gly

Διευθέτηση Conflicts:
Αν μια μετάλλαξη έχει πολλές υποβολές με διαφορετικές ταξινομήσεις, το πεδίο conflicting_interpretations θα περιέχει λίστα με όλες τις υποβολές.

Μείωση Όγκου:
Για μεγάλα datasets, εξετάστε το φιλτράρισμα μόνο για germline variants (OriginSimple == 'germline').



PARDEIGMA KWDIKA
import pandas as pd
import json

# load kai filtering dedomenwn
df_variants = pd.read_csv('variant_summary.txt', sep='\t', low_memory=False)
df_brca = df_variants[df_variants['GeneSymbol'] == 'BRCA1']

# merge me submission data
df_submissions = pd.read_csv('submission_summary.txt', sep='\t')
df_merged = pd.merge(df_brca, df_submissions, on='VariationID', how='left')

# efarmogh ACMG criteria (krithria.py)
df_merged['acmg_criteria'] = df_merged.apply(apply_acmg_criteria, axis=1)

# omadopoihsh gia conflicting interpretations
conflicting = df_merged.groupby('VariationID').apply(
    lambda x: x[['Submitter', 'ClinicalSignificance']].to_dict('records')
)

# teliko DataFrame
df_final = df_merged[[
    'VariationID', 'GeneSymbol', 'HGVS_c', 'HGVS_p',
    'ClinicalSignificance', 'ReviewStatus', 'PhenotypeList', 'Submitter'
]].drop_duplicates()

df_final['conflicting_interpretations'] = df_final['VariationID'].map(conflicting)
df_final['conflicting_interpretations'] = df_final['conflicting_interpretations'].apply(
    lambda x: json.dumps(x) if isinstance(x, list) else None
)

# save
df_final.to_csv('brca1_clinvar_processed.csv', index=False)


# RCV ACCESSION
import pandas as pd

# Φόρτωση δεδομένων
df = pd.read_csv('variant_summary.txt', sep='\t', low_memory=False)

# Φιλτράρισμα για BRCA1 και εξαγωγή RCVaccession
df_brca = df[df['GeneSymbol'] == 'BRCA1']
rcv_data = df_brca[['VariationID', 'RCVaccession']].copy()

# Διαχωρισμός πολλαπλών RCVaccession σε λίστα
rcv_data['RCVaccession'] = rcv_data['RCVaccession'].str.split('|')

# Επεκτεταμένο DataFrame με μία γραμμή ανά RCVaccession
rcv_expanded = rcv_data.explode('RCVaccession')


# RCV STO DATAFRAME
df_final = pd.merge(
    df_brca[['VariationID', 'GeneSymbol', 'HGVS_c', 'HGVS_p']],
    rcv_expanded,
    on='VariationID',
    how='left'
)

# Αποθήκευση
df_final.to_csv('brca1_with_rcv.csv', index=False)

# GENIKA
Αν μια μετάλλαξη έχει πολλά RCVs, αποφασίζουμε αν θα

Χωρίσουμε σε ξεχωριστές γραμμές (.explode()).

Αποθηκεύσουμε ως λίστα σε μία γραμμή (π.χ., ["RCV000123456.1", "RCV000789012.3"]).

Σύνδεση με PhenotypeList:
Κάθε RCVaccession αντιστοιχεί σε ένα PhenotypeList (κλινική κατάσταση). Για ανάλυση:
df_brca[['RCVaccession', 'PhenotypeList']].drop_duplicates()

Αξιοπιστία Ερμηνείας:

Τα RCVs με ReviewStatus = "Expert panel" είναι πιο αξιόπιστα.

Φιλτράρουμε με βάση το ReviewStatus αν χρειάζεται.