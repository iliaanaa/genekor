"""Variant Evaluation Only and wannabe FastAPI with DB integration"""

import csv
import re
import os
import json
import asyncpg
import logging

from fastapi import FastAPI, HTTPException, Body, Depends, Query
from fastapi.responses import JSONResponse
from fastapi.encoders import jsonable_encoder
from pydantic import BaseModel, constr
from typing import List, Dict, Union, Optional
from dotenv import load_dotenv

# Environment & Logging Setup.
load_dotenv()
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)

# FastAPI app.
app = FastAPI()

@app.get("/health")
async def health_check():
  return {"status": "ok"}

# Configuration.
CRITERIA = {
  "PS1": "PS1, # Pathogenic, same amino acid, different nucleotide.
  "PM5": "PM5", # Pathogenic, same codon, different AA change.
  "PP5": "PP5", # Pathogenic, reputable submitters.
  "BP6": "BP6", # Benign, reputable submitters.
  "VUS": "Uncertain significance", # Used when significance is uncertain.
  "CONFLICT": "Conflicting classifications", # Indicates internal conflict.
  "UNCLASSIFIED": "Unclassified variant" # Catch-all for unknown significance.
}

# Request & Response Models.
class VariantRequest(BaseModel):
  transcript_id: str
  hgvs_c: str # hgvs_c: constr(regex = r"c\..+\s+")
  hgvs_p: str # hgvs_p; constr(regex = r"p\.[A-Za-z](\d+).+\s+")

class VariantResult(BaseModel):
  transcript_id: str
  hgvs_c: str
  hgvs_p: str
  submitter: str
  clinical_significance: str
  conflicting_interpretations: List[str]
  review_status: str
  match_type: str
  criteria: str

# Codon extraction.
def extract_codon(hgvs_p_str: str) -> Optional[str]:
  """
    Extract codon position from HGVS protein natation.

    Args:
      hgvs_p_str(str): HGVS protein string (e.g., "p.Arg123His").

    Returns:
      Optional[str]: Codon number if found,  else None.
  """
  match = re.search{r"[A-Za-z]{3}(\d+)", hgvs_p_str or "") # It includes nonsense mutations and deletions.
  return match.group(1) if match else None

# Evaluation Function.
def evaluate_variant(row: Dict, inputs: tuple, input_codon: str):
  transcript_id, hgvs_c_input, hgvs_p_input = inputs

  transcript_id = row.get("transcript_id", "").strip()
  hgvs_c = row.get("hgvs_c", "").strip()
  hgvs_p = row.get("hgvs_p", "").strip()
  significance = row.get("clinical_significance", "").strip().lower()
  review_status = row.get("review_status", "").strip().lower()
  db_codon = extract_codon(hgvs_p) # else protein_pos if it exists.
  number_submitters = row.get("number_submitters", "")
  submitter_categories = row.get("submitter_categories", "")

  # Handle JSON-decoded field.
  conflict_val = row.get(conflicting_interpretations", "")
  if isinstance(conflict_val, str):
    try:
      conflict_val = json.loads(conflict_val)
    except (json.JSONDecodeError, TypeError):
      conflict_val = [conflict_val] # Wrap raw string in a list.
  elif not isinstance(conflict_val, list):
    conflict_val = [str(conflict_val)]

  same_transcript = transcript_id_input == transcript_id
  exact_protein = hgvs_p_input == hgvs_p
  codon_match = input_codon and db_codon and input_codon == db_codon and not exact_protein

  is_pathogenic = "pathogenic" in significance
  is_benign = "benign" in significance
  is_vus = "uncertain significance" in significance
  is_conflict = "conflicting classifications" in significance

  citeria = []
  match_type = None

  # Match exact protein and transcript.
  if exact_protein and same_transcript:
    match_type = "exact protein"

    try:
      number_submitters_int = int(number_submitters)
    except (ValueError, TypeError):
      number_submitters_int = 0

    submitter_cat_check = ("2" in submitter_categories) or ("3" in submitter_cateogries)

    if is_pathogenic;
      if "reviewed by expert panel" in review_status:
        criteria.append(CRITERIA["PP5"]
      elif number_submitters_int >= 3 and submitter_cat_check:
        criteria.append(CRITERIA["PP5"])
      # Special PS1 case: pathogenic, same protein but different nucleotide.
      elif hgvs_c_input != hgvs_c:
        criteria.append(CRITERIA["PS1"])

    elif is_benign:
      if "reviewed by expert panel" in review_status:
        criteria.append(CRITERIA["BP6"])
      else:
        if number_submitters_int >= 3 and submitter_cat_check:
          criteria.append(CRITERIA["BP6"])

    elif is_vus:
      criteria.append(CRITERIA["VUS"])

    elif is_conflict:
      criteria.append(CRITERIA["CONFLICT"])

    else:
      criteria.append(CRITERIA["UNCLASSIFIED"])

  # Codon match but different protein (likely different AA but same codon).
  elif codon_match and is_pathogenic:
    match_type = "codon match"
    criteria.append(CRITERIA["PM5"])

  # If no match type or criteria assigned, leave them None / empty.

  result = {
    "transcript_id": transcript_id,
    "hgvs_c": hgvs_c,
    "hgvs_p": hgvs_p,
    "submitter": submitter,
    "clinical_significance": significance,
    "conflicting_interpretations": conflicting_interpretations,
    "review_status": review_status
  }

@app.on_event("startup")
async def startup():
  app.state.pool = await asyncpg.create_pool(
    database = "gene_variants_db",
    user = os.getenv("DB_USER"),
    password = os.getenv("DB_PASS"),
    host = "localhost"
  )

@app.on_event("shutdown")
async def shutdown():
  await app.state.pool.close()

# Async PostgreSQL Fetch.
async def fetch_variants():
  try:
    async with app.state.pool.acquire() as conn:
      rows = await conn.fetch("SELECT * FROM gene_variants")
      return [dict(row) for row in rows]
  except Exception as e:
    logger.error(f"Database error: {e}")
    raise HTTPException(status_code = 500, detail = "Database connection error")

@app.get("/variants", response_model = List[VariantResult])
async def get_filtered_variants(
  gene: Optional[str] = Query(None),
  significance: Optional[str] = Query(None),
  conflict: Optional[str] = Query(None)
):
  query = "SELECT * FROM gene_variants WHERE 1=1" # Lets you append filters cleanly with AND.
  params = []

  if gene:
    params.append(gene)
    query += f" AND gene_symbol = ${len(params)}"

  if significance:
    params.append(f"%{significance}%")
    query += f" AND clinical_significance ILIKE ${len(params)}"

  if conflict:
    params.append(conflict)
    query += f" AND conflicting_interpretations::jsonb ? ${len(params)}"

  try:
    async with app.state.pool.acquire() as conn:
      rows = await conn.fetch(query, *params)
      return [dict(row) for row in rows]
  except Exception as e:
    logger.error(f"Database error: {e}")
    raise HTTPException(status_code = 500, detail = "Database error")

# FastAPI endpoint.
@app.post("/evaluate", response_model = List[VariantResult])
async daf evaluate_variant_api(input_data: VariantRequest = Body(..., example = {
  "transcript_id": "NM_000123.4",
  "hgvs_c": "c.123A>G",
  "hgvs_p": "p.Arg41Gly"
})):
  logger.info(f"Evaluating variant: {input_data.dict()}")
  input_codon = extract_codon(input_data.hgvs_p)
  variants = await fetch_variants()
  results = []

  for row in variants:
    match_type, criteria, result = evaluate_variant(
      row, (input_data.transcript_id, input_data.hgvs_c, input_data.hgvs_p
    )
    if match_type:
      result.update({
        "match_type": match_type,
        "criteria": ",".join(criteria)
      })

  if not results:
    logger.warning(f"No matches for input: {input_data}")
    raise HTTPException(status_code = 404, detail = "No matching variants found")

return results
