pm5_hits = df[
    (df['gene_symbol'] == user_gene) &
    (df['protein_pos'] == user_pos) &
    (df['p_hgvs'] != user_p_hgvs) &
    (df['molecular_consequence'] == "missense_variant") &
    (df['clinical_significance'].str.lower() == "pathogenic")
]

if not pm5_hits.empty and "missense" in user_consequence:
    crit.append("PM5")
