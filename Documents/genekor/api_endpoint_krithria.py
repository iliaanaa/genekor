from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import pandas as pd
import re

app = FastAPI()
df = pd.read_csv('brca_acmg_classified.csv')  # Προϋπόθεση: Το CSV έχει στήλη 'acmg_criteria' ως λίστα

class VariantQuery(BaseModel):
    gene: str                  # Π.χ., "BRCA1"
    variant: str               # Π.χ., "c.123G>T" ή "p.Val12Cys"

@app.post("/get_variant_info/")
async def get_variant_info(query: VariantQuery):
    # Αναγνώριση τύπου μετάλλαξης (DNA ή protein)
    if query.variant.startswith("c."):
        filtered = df[(df['gene'] == query.gene) & (df['dna_change'] == query.variant)]
    elif query.variant.startswith("p."):
        filtered = df[(df['gene'] == query.gene) & (df['protein_change'] == query.variant)]
    else:
        raise HTTPException(status_code=400, detail="Invalid variant format. Use 'c.' for DNA or 'p.' for protein changes.")

    if filtered.empty:
        raise HTTPException(status_code=404, detail=f"Variant {query.variant} not found for gene {query.gene}")

    variant_data = filtered.iloc[0].to_dict()

    # Ανάλυση κριτηρίων ACMG
    acmg_criteria = variant_data.get("acmg_criteria", [])
    acmg_explanations = {
        "PS1": "Same amino acid change as a known pathogenic variant (different nucleotide change).",
        "PM5": "Different amino acid change at the same position as a known pathogenic variant.",
        "PP5": "Pathogenic classification from reputable source without published evidence.",
        "BP6": "Benign classification from reputable source without published evidence."
    }

    # Δημιουργία λίστας με κριτήρια και εξηγήσεις
    criteria_details = [
        {"criterion": criterion, "explanation": acmg_explanations.get(criterion, "No description available.")}
        for criterion in acmg_criteria
    ]

    return {
        "gene": query.gene,
        "variant": query.variant,
        "clinical_significance": variant_data.get("clinical_significance"),
        "review_status": variant_data.get("review_status"),
        "acmg_criteria": criteria_details,
        "conflicting_interpretations": variant_data.get("conflicting_interpretations")
    }

