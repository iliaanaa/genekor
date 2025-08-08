 PM5 (ίδια θέση, διαφορετικό αμινοξύ):
Πάρε position από p_hgvs (π.χ. από p.Arg248Gln → 248)

Βρες variants με protein_position = 248 αλλά διαφορετικό p_hgvs

python
Copy
Edit
matches = db.query(GeneVariant).filter(
    GeneVariant.protein_position == position,
    GeneVariant.hgvs_p != p_hgvs,
    GeneVariant.significance.in_(["Pathogenic", "Likely pathogenic"])
).all()
if matches:
    criteria.append("PM5")
▶️ PP5 / BP6 (ίδιο c_hgvs):
python
Copy
Edit
same_c = db.query(GeneVariant).filter(
    GeneVariant.hgvs_c == c_hgvs,
    GeneVariant.hgvs_p != p_hgvs
).all()

if same_c:
    sigs = set(v.significance for v in same_c)
    if sigs.issubset({"Pathogenic", "Likely pathogenic"}):
        criteria.append("PP5")
    elif sigs.issubset({"Benign", "Likely benign"}):
        criteria.append("BP6")



✅ Πώς να το εφαρμόσεις (PS1 logic)
python
Copy
Edit
# Αν δεν υπάρχει το variant στη βάση:
variant = db.query(GeneVariant).filter_by(hgvs_c=c_hgvs, hgvs_p=p_hgvs).first()

# Εξετάζουμε το PS1:
ps1_matches = db.query(GeneVariant).filter(
    GeneVariant.hgvs_p == p_hgvs,
    GeneVariant.significance.in_(['Pathogenic', 'Likely pathogenic'])
).all()

if ps1_matches:
    criteria.append("PS1")