import psycopg2
import pandas as pd
import urllib.request
import re

def main():
    conn = psycopg2.connect(**DB_CONFIG)
    create_tables(conn)

    try:
        # === Κατέβασμα και προεπεξεργασία ClinVar ===
        variant_gz = "variant_summary.txt.gz"
        urllib.request.urlretrieve(CLINVAR_VARIANT_URL, variant_gz)
        df_final = process_clinvar_data(variant_gz)

        # === Εφαρμογή known/trusted κριτηρίων (π.χ. PS1, PM5, PP5, BP6) ===
        df_final = apply_acmg_criteria(df_final)

        # === Λήψη στοιχείων από χρήστη ===
        user_gene = input("Gene symbol: ").strip()
        user_c_hgvs = input("c.HGVS: ").strip()
        user_p_hgvs = input("p.HGVS (optional): ").strip() or None
        user_pos = None
        if user_p_hgvs:
            match = re.search(r'[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}', user_p_hgvs)
            if match:
                user_pos = int(match.group(1))

        # === Εξαγωγή παρόμοιων παραλλαγών ===
        same_c, same_p, same_pos = variant_assortments(df_final, user_gene, user_c_hgvs, user_p_hgvs, user_pos)
        same_c_groups = split_by_significance(same_c)
        same_p_groups = split_by_significance(same_p)
        same_pos_groups = split_by_significance(same_pos)
        support_tables = build_acmg_support_tables(same_c_groups, same_p_groups, same_pos_groups)

        # === Σήμανση grouping-based κριτηρίων ===
        df_final = mark_acmg_criteria(df_final, support_tables)

        # === Ενοποίηση κριτηρίων ===
        def combine_criteria(row):
            known = row.get("acmg_criteria", [])
            grouping = row.get("acmg_from_grouping", "")
            combined = set(known)
            if grouping:
                combined.update([x.strip() for x in grouping.split(";") if x])
            return "; ".join(sorted(combined))

        df_final["acmg_combined_criteria"] = df_final.apply(combine_criteria, axis=1)

        # === Εισαγωγή στη βάση ===
        insert_to_database(conn, df_final)

    finally:
        conn.close()
