from fastapi import FastAPI
import pandas as pd

app = FastAPI()
df = pd.read_csv('brca_acmg_variants.csv')

@app.get("/variants/{criterion}")
def get_variants(criterion: str):
    results = df[df['acmg_criteria'].apply(lambda x: criterion in x)]
    return {
        "criterion": criterion,
        "count": len(results),
        "variants": results.to_dict(orient='records')
    }


#aithma GET /variants/PS1