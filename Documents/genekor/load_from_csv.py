df_loaded = pd.read_csv('brca_acmg_variants.csv')
# Μετατροπή λιστών από strings (αν χρειάζεται)
import ast
df_loaded['acmg_criteria'] = df_loaded['acmg_criteria'].apply(ast.literal_eval)