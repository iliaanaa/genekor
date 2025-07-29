import psycopg2
import pandas as pd
import urllib.request
import re

def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)

    try:
        print("Ξεκίνημα script...")

        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)

        df_final = process_clinvar_data(variant_gz)

     

        user_pos = None
        if user_p_hgvs:
            match = re.search(r'[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}', user_p_hgvs)
            if match:
                user_pos = int(match.group(1))

        # === Βήμα 2: Εύρεση παρόμοιων παραλλαγών ===
        def variant_assortments(df, ref_gene, ref_c, ref_p=None, ref_pos=None):
            same_c = df[(df['GeneSymbol'] == ref_gene) & (df['HGVS_c'] == ref_c)]
            same_p = df[(df['GeneSymbol'] == ref_gene) & (df['HGVS_p'] == ref_p) & (df['HGVS_c'] != ref_c)] if ref_p else pd.DataFrame()
            same_pos = df[(df['GeneSymbol'] == ref_gene) & (df['protein_pos'] == ref_pos) & (df['HGVS_p'] != ref_p)] if ref_pos else pd.DataFrame()
            return same_c, same_p, same_pos

        same_c, same_p, same_pos = variant_assortments(df_final, user_gene, user_c_hgvs, user_p_hgvs, user_pos)

        # === Βήμα 3: Ομαδοποίηση κατά κλινική σημασία ===
        def split_by_significance(df):
            return {
                'Pathogenic': df[df['ClinicalSignificance'].str.contains('Pathogenic', na=False)],
                'Benign': df[df['ClinicalSignificance'].str.contains('Benign', na=False)],
                'VUS': df[df['ClinicalSignificance'].str.contains('Uncertain significance', na=False)]
            }

        same_c_groups = split_by_significance(same_c)
        same_p_groups = split_by_significance(same_p)
        same_pos_groups = split_by_significance(same_pos)

        # === Βήμα 4: Κατασκευή πίνακα υποψήφιων για ACMG ===
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

        # === Βήμα 5: Σήμανση στο df_final ===
        def mark_acmg_criteria(df, support):
            def determine_criteria(row):
                crit = []

                # PP5
                if row['HGVS_c'] in support['PP5']['HGVS_c'].values:
                    crit.append('PP5')

                # BP6
                if row['HGVS_c'] in support['BP6']['HGVS_c'].values:
                    crit.append('BP6')

                # PS1
                if pd.notna(row['HGVS_p']) and pd.notna(row['HGVS_c']):
                    ps1_matches = support['PS1'][
                        (support['PS1']['HGVS_p'] == row['HGVS_p']) &
                        (support['PS1']['HGVS_c'] != row['HGVS_c'])
                    ]
                    if not ps1_matches.empty:
                        crit.append('PS1')

                # PM5
                if pd.notna(row['protein_pos']) and pd.notna(row['HGVS_p']):
                    pm5_matches = support['PM5'][
                        (support['PM5']['protein_pos'] == row['protein_pos']) &
                        (support['PM5']['HGVS_p'] != row['HGVS_p'])
                    ]
                    if not pm5_matches.empty:
                        crit.append('PM5')

                return "; ".join(crit) if crit else ""

            df['acmg_from_grouping'] = df.apply(determine_criteria, axis=1)
            return df

        df_final = mark_acmg_criteria(df_final, support_tables)

        # === Συνέχεια εισαγωγής στη βάση ===
        insert_to_database(conn, df_final)

    finally:
        conn.close()
