import pandas as pd
import re
from sqlalchemy import create_engine

# ✅ 1. Συνάρτηση υπολογισμού consequence
def consequence_dna(hgvs_c: str) -> str:
    if not hgvs_c or hgvs_c.startswith('c.?'):
        return 'unknown'
    if re.search(r'[-+]\d+', hgvs_c):
        return 'intronic_variant'
    elif re.search(r'c\.\*\d+', hgvs_c):
        return '3_prime_UTR_variant'
    elif re.search(r'c\.-\d+', hgvs_c):
        return '5_prime_UTR_variant'
    elif '=' in hgvs_c:
        return 'synonymous_variant'
    elif 'del' in hgvs_c and 'ins' in hgvs_c:
        return 'indel'
    elif 'del' in hgvs_c:
        return 'deletion'
    elif 'ins' in hgvs_c:
        return 'insertion'
    elif 'dup' in hgvs_c:
        return 'duplication'
    elif '>' in hgvs_c:
        return 'substitution'
    return 'other'

# ✅ 2. Φόρτωσε ή φτιάξε το df_final
# Αν έχεις CSV:
# df_final = pd.read_csv("variants.csv")

# Ή αν το έχεις ήδη σε μεταβλητή:
# from your_script import df_final

# 👇 Dummy παράδειγμα, αντικατέστησέ το
data = {
    "gene_symbol": ["KLHL10", "KLHL10", "KLHL10"],
    "hgvs_c": ["c.123A>G", "c.456+1G>T", "c.789delT"],
    "clinical_significance": ["Pathogenic", "Benign", "Pathogenic"]
}
df_final = pd.DataFrame(data)

# ✅ 3. Υπολογισμός consequence
df_final["molecular_consequence"] = df_final["hgvs_c"].apply(consequence_dna)

# ✅ 4. Δημιουργία engine
engine = create_engine("postgresql+psycopg2://ilianam:genekor123!@localhost:5432/clinvar_db")

# ✅ 5. Drop και Replace πίνακα
with engine.connect() as conn:
    conn.execute("DROP TABLE IF EXISTS gene_variants")
    print("🗑️ Drop old gene_variants (αν υπήρχε).")

df_final.to_sql("gene_variants", engine, if_exists="replace", index=False)
print(f"✅ Εισαγωγή {len(df_final)} εγγραφών ολοκληρώθηκε.")

# ✅ 6. Τύπωσε μοναδικά consequences για έλεγχο
with engine.connect() as conn:
    results = conn.execute("SELECT DISTINCT molecular_consequence FROM gene_variants ORDER BY molecular_consequence")
    consequences = [r[0] for r in results]
    print("🎯 Unique consequences στη βάση:", consequences)
