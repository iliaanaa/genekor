from fastapi import FastAPI, Query
from typing import List, Optional
import psycopg2
import pandas as pd
from psycopg2.extras import RealDictCursor
import json
from collections import Counter
from collections import defaultdict
from new2 import apply_ps1_pm5_pp5_bp6
import re



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





def normalize(value):
    if value is None or (isinstance(value, str) and value.strip() == ""):
        return None
    return value


def build_variant_summary_with_acmg(result, gene_symbol, hgvs_c, hgvs_p, df_for_acmg):
    if result is None:
        protein_pos = None
        if hgvs_p:
            match = re.search(r"\d+", hgvs_p)
            if match:
                protein_pos = int(match.group())
        fake_row = pd.Series({
            "gene_symbol": gene_symbol,   # ❌ Διορθώθηκε: χρησιμοποιούμε gene_symbol του χρήστη
            "hgvs_c": hgvs_c,
            "hgvs_p": hgvs_p,
            "protein_pos": protein_pos
        })
        acmg_criteria = apply_ps1_pm5_pp5_bp6(fake_row, df_for_acmg)
        acmg_combined = "; ".join(acmg_criteria) if acmg_criteria else None
    else:
        acmg_criteria = result.get("acmg_criteria") or []
        acmg_combined = result.get("acmg_combined_criteria")
    
    return {
        "variant id": result.get("variation_id") if result else None,
        "gene": result.get("gene_symbol") if result else gene_symbol,
        "c.HGVS": hgvs_c,
        "p.HGVS": hgvs_p,
        "protein_pos": result.get("protein_pos") if result else protein_pos,
        "molecular_consequence": result.get("molecular_consequence") if result else None,
        "clinicalsignificance": result.get("clinicalsignificance") if result else "not provided",
        "review_status": result.get("review_status") if result else None,
        "other_fields": {
            "variation_id": result.get("variation_id") if result else None,
            "phenotype_list": result.get("phenotype_list") if result else None,
            "assembly": result.get("assembly") if result else None,
            "chromosome": result.get("chromosome") if result else None,
            "start_pos": result.get("start_pos") if result else None,
            "end_pos": result.get("end_pos") if result else None,
            "reference_allele": result.get("reference_allele") if result else None,
            "alternate_allele": result.get("alternate_allele") if result else None,
            "acmg_criteria": acmg_criteria,
            "conflicting_interpretations": result.get("conflicting_interpretations") if result else None,
            "rcvaccession": result.get("rcvaccession") if result else [],
            "last_updated": result.get("last_updated") if result else None,
            "last_evaluated": result.get("last_evaluated") if result else None,
            "acmg_combined_criteria": acmg_combined,
            "acmg_from_grouping": result.get("acmg_from_grouping") if result else None,
            "clinsigsimple": result.get("clinsigsimple") if result else None
        }
    }

@app.get("/user_classify_variant")
def user_classify_variant(
    gene_symbol: str = Query(..., description="Gene symbol (π.χ. BRCA1)"),
    hgvs_c: str = Query(..., description="c.HGVS notation (π.χ. c.123G>T)"),
    hgvs_p: Optional[str] = Query(None, description="Optional p.HGVS notation")
):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)

        # Αναζήτηση παραλλαγής
        cur.execute("SELECT * FROM gene_variants WHERE gene_symbol=%s AND hgvs_c=%s", (gene_symbol.strip(), hgvs_c.strip()))
        results = cur.fetchall()
        variant = results[0] if results else None

        # Όλες οι μεταλλάξεις για υπολογισμό ACMG
        cur.execute("SELECT gene_symbol, hgvs_c, hgvs_p, protein_pos, clinicalsignificance, clinsigsimple FROM gene_variants WHERE gene_symbol=%s", (gene_symbol,))
        df_for_acmg = pd.DataFrame(cur.fetchall())

        # Αν δεν υπάρχει η μετάλλαξη, φτιάχνουμε fake row
        if variant is None:
            protein_pos = None
            if hgvs_p:
                match = re.search(r"\d+", hgvs_p)
                if match:
                    protein_pos = int(match.group())
            row = pd.Series({
                "gene_symbol": gene_symbol,
                "hgvs_c": hgvs_c,
                "hgvs_p": hgvs_p,
                "protein_pos": protein_pos
            })
        else:
            row = pd.Series({
                "gene_symbol": variant.get("gene_symbol"),
                "hgvs_c": variant.get("hgvs_c"),
                "hgvs_p": variant.get("hgvs_p"),
                "protein_pos": variant.get("protein_pos")
            })

        # Υπολογισμός ACMG criteria (επιστρέφει dict με Yes/No)
        acmg_criteria = apply_ps1_pm5_pp5_bp6(row, df_for_acmg, return_dict=True)

        # Συγκεντρώνουμε το JSON
        return {
            "variant id": variant.get("variation_id") if variant else None,
            "gene": gene_symbol,
            "c.HGVS": hgvs_c,
            "p.HGVS": hgvs_p,
            "protein_pos": row.get("protein_pos"),
            "molecular_consequence": variant.get("molecular_consequence") if variant else None,
            "clinicalsignificance": variant.get("clinicalsignificance") if variant else "not provided",
            "review_status": variant.get("review_status") if variant else None,
            "other_fields": {
                "variation_id": variant.get("variation_id") if variant else None,
                "phenotype_list": variant.get("phenotype_list") if variant else None,
                "assembly": variant.get("assembly") if variant else None,
                "chromosome": variant.get("chromosome") if variant else None,
                "start_pos": variant.get("start_pos") if variant else None,
                "end_pos": variant.get("end_pos") if variant else None,
                "reference_allele": variant.get("reference_allele") if variant else None,
                "alternate_allele": variant.get("alternate_allele") if variant else None,
                "acmg_criteria": acmg_criteria,
                "conflicting_interpretations": variant.get("conflicting_interpretations") if variant else None,
                "rcvaccession": variant.get("rcvaccession") if variant else [],
                "last_updated": variant.get("last_updated") if variant else None,
                "last_evaluated": variant.get("last_evaluated") if variant else None,
                "acmg_combined_criteria": "; ".join([k for k,v in acmg_criteria.items() if v=="Yes"]) or None,
                "acmg_from_grouping": "; ".join([k for k,v in acmg_criteria.items() if v=="Yes"]) or None,
                "clinsigsimple": variant.get("clinsigsimple") if variant else None
            }
        }

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



@app.get("/acmg_criteria_bp6_pp5")
def get_acmg_criteria(
    gene_symbol: str = Query(..., description="Γονίδιο π.χ. KLHL10"),
    hgvs_p: Optional[str] = Query(None, description="Optional p.HGVS for PS1 evaluation"),
    hgvs_c: Optional[str] = Query(None, description="Optional c.HGVS for full evaluation")
):
    conn = None
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cur = conn.cursor(cursor_factory=RealDictCursor)

        cur.execute("""
            SELECT gene_symbol, hgvs_c, hgvs_p, protein_pos, clinicalsignificance, clinsigsimple
            FROM gene_variants
            WHERE gene_symbol = %s;
        """, (gene_symbol,))
        variants = cur.fetchall()

        if not variants:
            return {"error": f"Δεν βρέθηκαν μεταλλάξεις για το γονίδιο {gene_symbol}"}

        df = pd.DataFrame(variants)

        protein_pos = None
        if hgvs_p:
            match = re.search(r"\d+", hgvs_p)
            if match:
                protein_pos = int(match.group())

        row = pd.Series({
            "gene_symbol": gene_symbol,
            "hgvs_c": hgvs_c,
            "hgvs_p": hgvs_p,
            "protein_pos": protein_pos
        })
        criteria = apply_ps1_pm5_pp5_bp6(row, df)

        return {
            "gene": gene_symbol,
            "hgvs_c": hgvs_c,
            "hgvs_p": hgvs_p,
            "criteria": criteria,
            "conflict_score": len(criteria)
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
