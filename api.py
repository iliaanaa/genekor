from fastapi import FastAPI, Query
from typing import List, Optional
import psycopg2
import pandas as pd
from psycopg2.extras import RealDictCursor
import json
from collections import Counter
from collections import defaultdict


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


'''
# --- NEW --- #
@app.get("/user_classify_variant")  
def user_classify_variant(
    gene: str = Query(..., description="Gene symbol (e.g., BRCA1)"),
    c_hgvs: str = Query(..., description="c.HGVS notation (e.g., c.123G>T)"),
    p_hgvs: Optional[str] = Query(None, description="Optional p.HGVS notation (e.g., p.Val12Cys)")
):
    user_gene = gene.strip()
    user_c_hgvs = c_hgvs.strip()
    user_p_hgvs = p_hgvs.strip() if p_hgvs else None

    return {
        "message": "Variant input received.",
        "gene": user_gene,
        "c.HGVS": user_c_hgvs,
        "p.HGVS": user_p_hgvs
    }


'''

@app.get("/user_classify_variant")  
def user_classify_variant(
    gene: str = Query(..., description="Gene symbol (e.g., BRCA1)"),
    hgvs_c: str = Query(..., description="c.HGVS notation (e.g., c.123G>T)"),
    hgvs_p: Optional[str] = Query(None, description="Optional p.HGVS notation (e.g., p.Val12Cys)")
):
    conn = None 
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)

        query = """
        SELECT * 
        FROM gene_variants
        WHERE gene_symbol = %s AND hgvs_c = %s
        """
        params = [gene.strip(), hgvs_c.strip()]
        cur.execute(query, params)
        results = cur.fetchall()

        if not results:
            return{"error": f"No variant found for gene '{gene}' and c.HGVS '{hgvs_c}'"}
        
        variant = results[0]
        variant_summary ={
            "variant id": variant.get("variation_id"),
            "gene": variant.get("gene_symbol"),
            "c.HGVS": variant.get("hgvs_c"),
            "p.HGVS": hgvs_p.strip() if hgvs_p else variant.get("hgvs_p"),
            "protein_pos": variant.get("protein_pos"),
            "molecular_consequence": variant.get("molecular_consequence"),
            "clinicalsignificance": variant.get("clinicalsignificance"),
            "review_status": variant.get("review_status"),
             # "dbSNP": variant.get("rsid"),
            "other_fields": {k: v for k, v in variant.items() if k not in [
                "gene_symbol", "transcript_id", "hgvs_c", "hgvs_p", "protein_pos", 
                "molecular_consequence", "clinicalsignificance", 
                "review_status"
            ]}
        }

        return variant_summary
    except Exception as e:
        return{"error": str(e)}
    finally:
        if conn:
            conn.close()


    # ---  NEW --- #
def calculate_pp5_bp6_from_summary(variants):
    """
    Υπολογισμός PP5 / BP6 από variant summary:
    - Δεν κοιτάμε τα ονόματα των submitters
    - Χρησιμοποιούμε μόνο count + consistency
    """
    significance_counts = Counter([v['clinicalsignificance'] for v in variants if v['clinicalsignificance']])

    # Υπολογισμός PP5
    pp5 = significance_counts.get("Pathogenic", 0) >= 2 and len(significance_counts) == 1

    # Υπολογισμός BP6
    bp6 = significance_counts.get("Benign", 0) >= 2 and len(significance_counts) == 1

    return pp5, bp6


'''
def calculate_pp5_bp6_from_summary(variants):
    """
    Υπολογισμός PP5 / BP6 από variant summary:
    - Δεν κοιτάμε τα ονόματα των submitters
    - Χρησιμοποιούμε μόνο count + consistency
    """
    significance_counts = Counter([v['clinicalsignificance'] for v in variants if v['clinicalsignificance']])

    # Υπολογισμός PP5
    pp5 = significance_counts.get("Pathogenic", 0) >= 2 and len(significance_counts) == 1

    # Υπολογισμός BP6
    bp6 = significance_counts.get("Benign", 0) >= 2 and len(significance_counts) == 1

    return pp5, bp6
'''


@app.get("/acmg_criteria_bp6_pp5")
def get_acmg_criteria(gene: str = Query(..., description="Γονίδιο π.χ. KLHL10")):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)

        cur.execute("""
            SELECT clinicalsignificance
            FROM gene_variants
            WHERE gene_symbol = %s;
        """, (gene,))
        variants = cur.fetchall()

        if not variants:
            return {"error": f"Δεν βρέθηκαν μεταλλάξεις για το γονίδιο {gene}"}

        pp5, bp6 = calculate_pp5_bp6_from_summary(variants)
        conflict_score = int(pp5) + int(bp6)

        return {
            "gene": gene,
            "PP5": pp5,
            "BP6": bp6,
            "conflict_score": conflict_score
        }

    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn:
            conn.close()


@app.get("/variants_by_genomic_range")
def variants_by_genomic_range(
    gene: str = Query(..., description="Gene symbol (e.g., TP53)"),
    start: int = Query(..., description="Start genomic position (e.g., 7668402)"),
    end: int = Query(..., description="End genomic position (e.g., 7687550)")
):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)

        query = """
        SELECT gene_symbol, protein_pos, hgvs_c, hgvs_p, clinicalsignificance, molecular_consequence
        FROM gene_variants
        WHERE gene_symbol = %s
          AND start_pos >= %s
          AND end_pos <= %s
        """
        cur.execute(query, (gene, start, end))
        results = cur.fetchall()

        return {"results": results}
    
    except Exception as e:
        return {"error": str(e)}
    
    finally:
        if conn:
            conn.close()


'''
@app.get("/variants_by_protein_pos")
def variants_by_protein_pos(
    gene: str = Query(..., description="Gene symbol (e.g. TP53)"),
    start_pos: int = Query(..., description="Start of protein position range (e.g. 350)"),
     end_pos: Optional[int] = Query(None, description="End of protein position range (optional)")
    ):
   # end_pos: int = Query(..., description="End of protein position range (e.g. 400)")

    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)

        query = """
        SELECT gene_symbol, protein_pos, hgvs_c, hgvs_p
        FROM gene_variants
        WHERE gene_symbol = %s
          AND protein_pos BETWEEN %s AND %s
        """
        cur.execute(query, (gene, start_pos))
        results = cur.fetchall()

        return {"results": results}
    
    except Exception as e:
        return {"error": str(e)}
    
    finally:
        if conn:
            conn.close()

'''


@app.get("/variants_by_protein_pos")
def get_variants_by_protein_pos(
    gene: str = Query(..., description="Gene symbol (e.g. TP53)"),
    start_pos: int = Query(..., description="Start of protein position range (e.g. 100)"),
    end_pos: Optional[int] = Query(None, description="End of protein position range (optional)")
):
    conn = psycopg2.connect(**DB_CONFIG)
    cur = conn.cursor()
    
    if end_pos is not None:
        query = """
        SELECT * FROM gene_variants
        WHERE gene_symbol = %s
        AND protein_pos BETWEEN %s AND %s
        """
        cur.execute(query, (gene, start_pos, end_pos))
    else:
        query = """
        SELECT * FROM gene_variants
        WHERE gene_symbol = %s
        AND protein_pos = %s
        """
        cur.execute(query, (gene, start_pos))

    rows = cur.fetchall()
    cur.close()
    conn.close()
    return rows


@app.get("/variant_counts")
def get_variant_counts(
    gene: str = Query(...,description="Γονίδιο π.χ. KLHL10"),
    consequence: Optional[str] = Query(None, description="Τύπος μετάλλαξης π.χ. missense"),
    significance: Optional[str] = Query(None, description="Παθογένεια π.χ. Pathogenic"),
    protein_start: Optional[int] = Query(None, description="Αρχή εύρους θέσης πρωτεΐνης"),
    protein_end: Optional[int] = Query(None, description="Τέλος εύρους θέσης πρωτεΐνης"),
    exact_position: Optional[int] = Query(None, description="Συγκεκριμένη θέση πρωτεΐνης")
):
    """
    Επιστρέφει πλήθος μεταλλάξεων με βάση φίλτρα: τύπος, παθογένεια, θέση πρωτεΐνης
    """
    conn = None
    try: 
        conn=psycopg2.connect(**DB_CONFIG)
        cur=conn.cursor()

        query="SELECT COUNT(*) FROM gene_variants WHERE gene_symbol = %s"
        params = [gene]

        if consequence:
            query += " AND molecular_consequence ILIKE %s"
            params.append(f"%{consequence}%")

        if significance:
            query += " AND clinicalsignificance ILIKE %s"
            params.append(f"%{significance}%")


        if exact_position is not None:
            query += " AND protein_pos = %s"
            params.append(exact_position)
        elif protein_start is not None and protein_end is not None:
            query += " AND protein_pos BETWEEN %s AND %s"
            params.extend([protein_start, protein_end])

        cur.execute(query, tuple(params))
        count = cur.fetchone()[0]
        return {"count": count}
    
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn:
            conn.close()

# 1. Summary by molecular consequence
@app.get("/summary")
def summary_by_consequence(gene: str = Query(...)):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)
        cur.execute("""
            SELECT molecular_consequence, COUNT(*) as count
            FROM gene_variants
            WHERE gene_symbol = %s
            GROUP BY molecular_consequence
            ORDER BY count DESC;
        """, (gene,))
        results = cur.fetchall()
        summary = {row['molecular_consequence']: row['count'] for row in results}
        return {"gene": gene, "summary": summary}
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn: conn.close()



#Ομαδοποίηση με βάση παθογένεια
@app.get("/significance_summary")
def significance_summary(gene: str = Query(...)):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)
        cur.execute("""
            SELECT clinicalsignificance, COUNT(*) as count
            FROM gene_variants
            WHERE gene_symbol = %s
            GROUP BY clinicalsignificance
            ORDER BY count DESC;
        """, (gene,))
        rows = cur.fetchall()
        return {"gene": gene, "summary": {r['clinicalsignificance']: r['count'] for r in rows}}
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn: conn.close()


# Αναζήτηση με τύπο παραλλαγής + παθογένεια
def variant_counts(gene: str, consequence: Optional[str] = None, significance: Optional[str] = None):
    conn = None    
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor()
        query = "SELECT COUNT(*) FROM gene_variants WHERE gene_symbol = %s"
        params = [gene]
        if consequence:
            query += " AND molecular_consequence ILIKE %s"
            params.append(f"%{consequence}%")
        if significance:
            query += " AND clinicalsignificance ILIKE %s"
            params.append(f"%{significance}%")
        cur.execute(query, tuple(params))
        count = cur.fetchone()[0]
        return {"gene": gene, "count": count}
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn: conn.close()
    
'''
@app.get("/variants_by_position")
def variants_by_position(gene: str, min_pos: Optional[int] = None, max_pos: Optional[int] = None):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)
        query = "SELECT * FROM gene_variants WHERE gene_symbol = %s"
        params = [gene]
        if min_pos is not None:
            query += " AND start_pos >= %s"
            params.append(min_pos)
        if max_pos is not None:
            query += " AND end_pos <= %s"
            params.append(max_pos)
        cur.execute(query, tuple(params))
        return cur.fetchall()
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn: conn.close()
'''

@app.get("/available_genes")
def available_genes():
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor()
        cur.execute("SELECT DISTINCT gene_symbol FROM gene_variants ORDER BY gene_symbol")
        genes = [r[0] for r in cur.fetchall()]
        return {"genes": genes}
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn: conn.close()

# 6. Available consequence types
@app.get("/available_consequences")
def available_consequences():
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor()
        cur.execute("SELECT DISTINCT molecular_consequence FROM gene_variants ORDER BY molecular_consequence")
        types = [r[0] for r in cur.fetchall()]
        return {"consequences": types}
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn: conn.close()

# 7. Search with multiple filters
@app.get("/search_variants")
def search_variants(
    gene: Optional[str] = None,
    consequence: Optional[str] = None,
    significance: Optional[str] = None,
    protein_pos: Optional[int] = None
):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)
        query = "SELECT * FROM gene_variants WHERE 1=1"
        params = []
        if gene:
            query += " AND gene_symbol = %s"
            params.append(gene)
        if consequence:
            query += " AND molecular_consequence ILIKE %s"
            params.append(f"%{consequence}%")
        if significance:
            query += " AND clinicalsignificance ILIKE %s"
            params.append(f"%{significance}%")
        if protein_pos is not None:
            query += " AND protein_pos = %s"
            params.append(protein_pos)
        cur.execute(query, tuple(params))
        return cur.fetchall()
    except Exception as e:
        return {"error": str(e)}
    finally:
        if conn: conn.close()
