"""Variant Evaluation - ACMG 2015"""
# Applies ACMG criteria: PS1, PM5, PP5, BP6.
# review_status "criteria provided, single submitter" can be used to classify a variant as PP5 or BP6.

import csv
import re

# --- Configuration ---
gene = "brca1"
INPUT_FILE = f"{gene}_variants.tsv"
OUTPUT_FILE = f"{gene}_variants_evaluated.tsv"

# The dictionary has been created (known_pathogenic{...}).
# --- Example: known_pathogenic dictionary ---
# known_pathogenic = {
#   'p.Arg175His' : {
#       'Clinical Significance': 'Pathogenic',
#       'HGVS_c': 'c.524G>A',
#       'HGVS_p': 'p.Arg175His',
#       'codon': '175'
#   },
#   ...
# }

# --- ACMG Evaluation ---
criteria = []

with open(INPUT_FILE, "r", encoding = "utf-8") as infile:
    reader = csv.DictReader(infile, delimiter = "\t")
    fieldnames = reader.fieldnames + ["acmg_evidence"] # Add ACMG column.

    for row in reader:
        acmg_evidence = []
        review_status = row.get("review_status", "").strip().lower()
        clinical_significance = row.get("clinical_significance","").strip().lower()
        hgvs_c = row.get("hgvs_c", "").strip()
        hgvs_p = row.get("hgvs_p", "").strip()
        transcript_id = row.get("transcript_id", "").strip()

        # Skip exact match in known_pathogenic.
        if hgvs_c in known_pathogenic or hgvs_p in known_pathogenic:
            continue

        # --- PS1: Same protein change, different nucleotide change ---
        if hgvs_p and hgvs_p in known_pathogenic:
            if hgvs_c != known_pathogenic[hgvs_p]["HGVS_c"]:
                acmg_evidence.append("PS1")

        # --- PM5: Different missense change at same codon as a known pathogenic ---
        codon = None
        if hgvs_p:
            codon_match = re.search(r"[A-Za-z]{3}(\d+)[A-Za-z]{3}", hgvs_p)
            if codon_match:
                codon = codon_match.group(1)
                for known in known_pathogenic.values():
                    if known.get("codon") == codon and known.get("HGVS_p") != hgvs_p:
                        acmg_evidence.append("PM5")
                        break

        # --- PP5: Single submitter, claims pathogenicity ---
        if review_status == "criteria provided, single submitter" and "pathogenic" in clinical_significance:
            acmg_evidence.append("PP5")

        # --- BP6: Single submitter, claims benign ---
        if review_status == "criteria provided, single submitter" and "benign" in clinical_significance:
            acmg_evidence.append("BP6")

        # Debug info: unmatched variants.
        if not acmg_evidence:
            print(f"[INFO] No ACMG criteria matched for: {transcript_id}\t{hgvs_c}\t{hgvs_p or 'N/A'}\n")

        # Store result.
        row["acmg_evidence"] = ",".join(acmg_evidence) if acmg_evidence else "None"
        criteria.append(row)

# Write to output file.
with open(OUTPUT_FILE, "w", encoding = "utf-8", newline = "") as outfile:
    writer = csv.DictWriter(outfile, fieldnames = fieldnames, delimiter = "\t")
    writer.writeheader()
    writer.writerows(criteria)

print(f"ACMG evaluation complete.\n Results written to: {OUTPUT_FILE}.\n")
