from fastapi import FastAPI, Query
from typing import List
import psycopg2
import pandas as pd
from psycopg2.extras import RealDictCursor
import json

app=FastAPI()



# Ρυθμίσεις σύνδεσης βάσης
DB_CONFIG={
    "dbname":"clinvar_db",
    "user":"ilianam",
    "password":"genekor123!",
    "host":"localhost",
    "port":5432
}

# Ρίζα - απλό health check
@app.get("/")
def health_check():
    return{"status":"clinvar api is running"}


# Ερώτημα: πάρε μεταλλάξεις για ένα γονίδιο
@app.get("/variants")
def get_variants(
    gene: str = Query(..., description="Γονίδιο π.χ. KLHL10"),
    max_results: int = Query(100, ge=1, le=1000)
):
    """
    Επιστρέφει τις μεταλλάξεις για το συγκεκριμένο γονίδιο
    """
    conn = None

    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)
        cur.execute("""
            SELECT variation_id, hgvs_c, hgvs_p, clinical_significance,
                   molecular_consequence, review_status, acmg_criteria
            FROM gene_variants
            WHERE gene_symbol = %s
            LIMIT %s;
        """, (gene, max_results))
        
        results = cur.fetchall()
        return results

    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn:
            conn.close()

            