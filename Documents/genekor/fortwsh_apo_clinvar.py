import pandas as pd

# Φόρτωση TSV αρχείου ClinVar
df = pd.read_csv('variant_summary.txt', sep='\t', low_memory=False)

# Φιλτράρισμα για BRCA1/2
df_brca = df[df['GeneSymbol'].isin(['BRCA1', 'BRCA2'])]