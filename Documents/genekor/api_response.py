from fastapi import FastAPI
import pandas as pd

app = FastAPI()
df = pd.read_csv('brca_acmg_variants.csv') #Ορίζει ένα GET endpoint με μια παράμετρο criterion (π.χ., PS1, PM5).

@app.get("/variants/{criterion}")
def get_variants(criterion: str):
    results = df[df['acmg_criteria'].apply(lambda x: criterion in x)] #Φιλτράρει το DataFrame για να βρει μεταλλάξεις όπου το criterion υπάρχει στη λίστα acmg_criteria.
    return {
        "criterion": criterion,
        "count": len(results),
        "variants": results.to_dict(orient='records') #Μετατρέπει τα αποτελέσματα σε λίστα από dictionaries (κάθε dictionary αντιστοιχεί σε μια μετάλλαξη)
    }


#aithma GET /variants/PS1