'''
META APO df_final = process_clinvar_data
PRIN APO insert_to_database
'''
def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)
    
    try:
        print("Ξεκίνημα script...")  # Για επιβεβαίωση ότι τρέχει
        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)

        df_final = process_clinvar_data(variant_gz)

        ### NEW - Step 1: Get user input
        print("\nΔώσε τα παρακάτω για την παραλλαγή στόχο:")
        user_gene = input("Gene symbol (e.g., BRCA1): ").strip()
        user_c_hgvs = input("c.HGVS (e.g., c.123G>T): ").strip()
        user_p_hgvs = input("p.HGVS (e.g., p.Val12Cys) [Optional]: ").strip()
        user_p_hgvs = user_p_hgvs if user_p_hgvs else None  # None if left blank

        user_pos = None
        if user_p_hgvs:
            match = re.search(r'[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}', user_p_hgvs)
            if match:
                user_pos = int(match.group(1))

        ### NEW - Step 2: Build assortments
        def variant_assortments(df, ref_gene, ref_c, ref_p=None, ref_pos=None):
            same_c = df[(df['GeneSymbol'] == ref_gene) & (df['HGVS_c'] == ref_c)]
            same_p = df[(df['GeneSymbol'] == ref_gene) & (df['HGVS_p'] == ref_p) & (df['HGVS_c'] != ref_c)] if ref_p else pd.DataFrame()
            same_pos = df[(df['GeneSymbol'] == ref_gene) & (df['protein_pos'] == ref_pos) & (df['HGVS_p'] != ref_p)] if ref_pos else pd.DataFrame()
            return same_c, same_p, same_pos

        same_c, same_p, same_pos = variant_assortments(df_final, user_gene, user_c_hgvs, user_p_hgvs, user_pos)

        ### NEW - Step 3: Split by significance
        def split_by_significance(df):
            return {
                'Pathogenic': df[df['ClinicalSignificance'].str.contains('Pathogenic', na=False)],
                'Benign': df[df['ClinicalSignificance'].str.contains('Benign', na=False)],
                'VUS': df[df['ClinicalSignificance'].str.contains('Uncertain significance', na=False)]
            }

        same_c_groups = split_by_significance(same_c)
        same_p_groups = split_by_significance(same_p)
        same_pos_groups = split_by_significance(same_pos)

        ### NEW - Step 4: Build ACMG support candidates
        def build_acmg_support_tables(same_c_groups, same_p_groups, same_pos_groups):
            pp5_candidates = same_c_groups['Pathogenic'][['HGVS_c']]
            bp6_candidates = same_c_groups['Benign'][['HGVS_c']]
            ps1_candidates = same_p_groups['Pathogenic'][['HGVS_p', 'HGVS_c']]
            pm5_candidates = same_pos_groups['Pathogenic'][['HGVS_p', 'HGVS_c', 'protein_pos']]
            return {
                'PP5': pp5_candidates,
                'BP6': bp6_candidates,
                'PS1': ps1_candidates,
                'PM5': pm5_candidates
            }

        support_tables = build_acmg_support_tables(same_c_groups, same_p_groups, same_pos_groups)

        print("\n=== Υποψήφιες για ACMG κριτήρια ===")
        for crit, df_support in support_tables.items():
            print(f"{crit}: {len(df_support)} υποψήφιες")
            print(df_support.head())

        ### NEW - Step 5: Επισήμανση στο df_final (τροποποιημένη για string)
        def mark_acmg_criteria(df, support):
            def determine_criteria(row):
                crit = []
                if row['HGVS_c'] in support['PP5']['HGVS_c'].values:
                    crit.append('PP5')
                if row['HGVS_c'] in support['BP6']['HGVS_c'].values:
                   crit.append('BP6')
                if row['HGVS_p'] in support['PS1']['HGVS_p'].values:
                   crit.append('PS1')
                if row['protein_pos'] in support['PM5']['protein_pos'].values:
                  crit.append('PM5')
                return "; ".join(crit) if crit else ""  # <=== returns string or empty string
            
            df['acmg_from_grouping'] = df.apply(determine_criteria, axis=1)
            return df
        
        df_final = mark_acmg_criteria(df_final, support_tables)

        # --- original script continues here ---
