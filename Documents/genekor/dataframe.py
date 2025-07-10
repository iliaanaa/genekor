import pandas as pd

# Παράδειγμα DataFrame με κριτήρια ACMG
data = [
    {
        "variant_id": 1,
        "gene": "BRCA1",
        "dna_change": "c.1510A>G",
        "protein_change": "p.Arg504Gly",
        "clinical_significance": "Pathogenic",
        "review_status": "Expert panel",
        "submitter": "ENIGMA",
        "acmg_criteria": ["PS1", "PP5"],  # PS1: Ίδιο protein change με γνωστό pathogenic, PP5: Αξιόπιστη πηγή
        "conflicting_interpretations": None
    },
    {
        "variant_id": 2,
        "gene": "BRCA1",
        "dna_change": "c.1510A>T",
        "protein_change": "p.Arg504Ser",
        "clinical_significance": "Likely pathogenic",
        "review_status": "Conflicting interpretations",
        "submitter": "LabX",
        "acmg_criteria": ["PM5"],  # PM5: Διαφορετικό amino acid change σε γνωστή θέση
        "conflicting_interpretations": [
            {"submitter": "LabA", "classification": "Pathogenic"},
            {"submitter": "LabB", "classification": "Benign"}
        ]
    },
    {
        "variant_id": 3,
        "gene": "BRCA2",
        "dna_change": "c.123G>A",
        "protein_change": "p.Trp41*",
        "clinical_significance": "Benign",
        "review_status": "Expert panel",
        "submitter": "ENIGMA",
        "acmg_criteria": ["BP6"],  # BP6: Benign από αξιόπιστη πηγή
        "conflicting_interpretations": None
    }
]

df = pd.DataFrame(data)