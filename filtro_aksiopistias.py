# --- PRIN THN MAIN ---

# Synarthsh poy filtrarei aksiopista VariationID ---
def get_reliable_variation_ids_from_variant_summary(path: str) -> set:
    """
    Epistrefei to synolo ton aksiopistwn VariationID:
    - reviewed by expert panel ή
    - >=3 submitters σε SubmitterCategory 2 ή 3 που συμφωνούν
    """
    df['SubmitterCategories'] = pd.to_numeric(df['SubmitterCategories'], errors='coerce')
    df['VariationID'] = pd.to_numeric(df['VariationID'], errors='coerce')
    df = df.dropna(subset=['VariationID', 'SubmitterCategories'])

    # Expert panel
    expert_panel_ids = set(df[df['ReviewStatus'] == 'reviewed by expert panel']['VariationID'])

    # Oxi expert panel
    # Filtraroyme gia SubmitterCategories 2 ή 3 kai xwris conflicting
    df_valid = df[
        (~df['VariationID'].isin(expert_panel_ids)) &
        (df['SubmitterCategories'].isin([2, 3])) &
        (~df['ClinicalSignificance'].str.lower().str.contains("conflict"))
    ]

    # Grouping me bash VariationID kai ClinicalSignificance
    grouped = (
        df_valid
        .groupby(['VariationID', 'ClinicalSignificance'])
        .size()
        .reset_index(name='count')
    )

    # Kratame osa exoyn toyl. 3 submitters me thn idia clinical shmasia
    reliable_non_expert = grouped[grouped['count'] >= 3]

    # Elegxoyme na yparxei mono mia ClinicalSignificance ana VariationID
    counts_per_id = reliable_non_expert.groupby('VariationID')['ClinicalSignificance'].nunique().reset_index()
    consistent_ids = counts_per_id[counts_per_id['ClinicalSignificance'] == 1]['VariationID']

    reliable_ids = expert_panel_ids.union(set(consistent_ids))

    print(f"Aksiopista VariationIDs: {len(reliable_ids)}")
    return reliable_ids

#--- STHN MAIN ---
# Filtro aksiopistias
reliable_ids = get_reliable_variation_ids_from_variant_summary(variant_gz.replace('.gz', ''))
df_final = df_final[df_final['VariationID'].isin(reliable_ids)]




"""
    # Φιλτράρουμε για clinical testing και curation
    def has_trusted_submitters(category_str):
        if pd.isna(category_str):
            return False
        parts = category_str.lower().split('|')
        return sum(cat in ['clinical testing', 'curation'] for cat in parts) >= 3

    df_trusted = df[
        (~df['VariationID'].isin(expert_panel_ids)) &
        (df['ClinicalSignificance'].str.lower().str.contains("conflicting") == False) &
        (df['SubmitterCategories'].apply(has_trusted_submitters))
    ]

    trusted_ids = set(df_trusted['VariationID'])

    reliable_ids = expert_panel_ids.union(trusted_ids)
    print(f"Συνολικά αξιόπιστα VariationIDs: {len(reliable_ids)}")

    return reliable_ids
"""
    """
    df_valid = df[
        (~df['VariationID'].isin(expert_panel_ids)) &
        (df['SubmitterCategory'].isin([2, 3]))
    ]

    # Ομαδοποίηση: μετράμε submitters ανά VariationID + ClinicalSignificance
    grouped = (
        df_valid
        .groupby(['VariationID', 'ClinicalSignificance'])
        .agg(num_submitters=('Submitter', 'nunique'))
        .reset_index()
    )

    # Κρατάμε μόνο όσα έχουν ≥3 συμφωνούντες submitters
    reliable = grouped[grouped['num_submitters'] >= 3]

    # Ελέγχουμε ότι κάθε VariationID έχει μόνο μια μοναδική ClinicalSignificance
    sig_count = reliable.groupby('VariationID')['ClinicalSignificance'].nunique().reset_index()
    sig_count = sig_count[sig_count['ClinicalSignificance'] == 1]

    non_expert_reliable_ids = set(sig_count['VariationID'])

    all_reliable = expert_panel_ids.union(non_expert_reliable_ids)
    print(f"Αξιόπιστα variation_ids: {len(all_reliable)}")

    return all_reliable
"""