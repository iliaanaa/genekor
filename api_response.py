from fastapi import FastAPI
import pandas as pd
import ast

app = FastAPI()

# Φόρτωμα του CSV
df = pd.read_csv('brca_acmg_variants.csv')

# Αν η στήλη acmg_criteria περιέχει λίστες ως string, τις μετατρέπουμε σε Python λίστες
df['acmg_criteria'] = df['acmg_criteria'].apply(
    lambda x: ast.literal_eval(x) if isinstance(x, str) else []
)

#endpoint: Δίνεις μια μετάλλαξη και επιστρέφει τύπο και κριτήρια
@app.get("/variant_info/{mutation}")
def get_variant_info(mutation: str):
    # Φιλτράρει με βάση την στήλη 'ProteinChange'
    match = df[df['ProteinChange'] == mutation]

    if match.empty:
        return {
            "mutation": mutation,
            "message": "Η μετάλλαξη δεν βρέθηκε στη βάση δεδομένων."
        }

    # Παίρνει την πρώτη καταχώρηση
    row = match.iloc[0]

    return {
        "mutation": mutation,
        "variant_type": row['variant_type'],
        "acmg_criteria": row['acmg_criteria'],
    }
