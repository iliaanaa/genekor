"""Enhanced Evalutaion of Variants with ACMG criteria"""
import csv
import re
import os
import json

# --- Configuration ---
gene = "gene_symbol"
db = f"{gene}_variants.tsv"

# Set of reputable submitters.
reputable_submitters = {
    "GeneDx",
    "Invitae",
    "Ambry Genetics",
    "Color Health",
    "Blueprint Genetics",
    "Laboratory for Molecular Medicine",
    "University of Chicago Genetic Services"
}

# --- Input Collection ---
def collect_input():
    transcript_id_input = input("Transcript ID: ").strip()
    hgvs_c_input = input("coding variant: ").strip()
    hgvs_p_input = input("protein variant: ").strip()
    return transcript_id_input, hgvs_c_input, hgvs_p_input

# --- Codon Extraction ---
def extract_codon(hgvs_p_str):
    match = re.search(r"[A-Za-z]{3}(\d+)[A-Za-z]{3}", hgvs_p_str)
    return match.group(1) if match else None

# --- Criteria Evaluation ---
def evaluat_variant(row, inputs, input_codon):
    criteria = []
    transcript_id_input, hgvs_c_input, hgvs_p_input = inputs
    # Extract fields from row.
    transcript_id = row.get("transcript_id", "").strip()
    hgvs_c = row.get("hgvs_c", "").strip()
    hgvs_p = row.get("hgvs_p", "").strip()
    clinical_significance = row.get("clinical_significance", "").strip().lower()
    review_status = row.get("review_status", "").strip().lower()
    submitter = row.get("submitter", "").strip()
    
    db_codon = extract_codon(hgvs_p)
    
    match_type = None # for tracking: exact, codon_match, skip

    # Case 1: Exact protein match.
    if transcript_id_input == transcript_id and hgvs_p_input == hgvs_p:
        match_type = "exact protein"
        if hgvs_c_input == hgvs_c:
            # PP5.
            if "pathogenic" in clinical_significance:
                if review_status == "criteria provided, single submitter":
                    if submitter in reputable_submitters:
                        criteria.append("PP5")
                    else:
                        criteria.append("Not a reputable source")
                else:
                    criteria.append("Cannot apply PP5: review_status not single submitter")
            # BP6.
            elif "benign" in clinical_significance:
                if review_status == "criteria provided, single submitter":
                    if submitter in reputable_submitters:
                        criteria.append("BP6")
                    else:
                        criteria.append("Not a reputable source")
                else:
                    criteria.append("Cannot apply BP6: review_status not single submitter")
            elif "uncertain significance" in clinical_significance:
                criteria.append("VUS")
            elif "conflicting classifications" in clinical_significance:
                criteria.append("Conflicting classifications")
            else:
                criteria.append("Unclassified variant")
        else:
            # PS1.
            if "pathogenic" in clinical_significance:
                criteria.append("PS1")
            elif "uncertain significance" in clinical_significance:
                criteria.append("VUS")
            elif "conflicting classifications" in clinical_significance:
                criteria.append("Conflicting classifications")
            else:
                criteria.append("Likely benign")

    # Case 2: Same codon, different amino acid (PM5).
    elif transcript_id_input == transcript_id:
        match_type = "codon_match"
        # PM5.
        if input_codon == db_codon and hgvs_p_input != hgvs_p:
            if "pathogenic" in clinical_significance:
                criteria.append("PM5")

    return match_type, criteria
    
# --- Main Evaluation ---
def evaluate_all_variants():
    inputs = collect_input()
    transcript_id_input, hgvs_c_input, hgvs_p_input = inputs
    input_codon = extract_codon(hgvs_p_input)
    
    found_any = False
    all_matches = []
    
    if not os.path.exists(db_filename):
        print(f"Error: Database file '{db_filename}' not found.")
        return
        
    with open(db_filename, "r", encoding = "utf-8") as infile:
        reader = csv.DictReader(infile, delimiter = "\t")
        for row in reader:
            match_type, criteria = evaluate_variant(row, inputs, input_codon)
            
            if match_type:
                found_any = True
                match_info = {
                    "transcript_id": row.get("transcript_id", ""),
                    "hgvs_c": row.get("hgvs_c", ""),
                    "hgvs_p": row.get("hgvs_p", ""),
                    "submitter": row.get("submitter", ""),
                    "clinical_significance": row.get("clinical_significance", "")
                    "review_status": row.get("review_status", ""),
                    "evidence": ", ".join(criteria),
                    "match_type": match_type,
                }
                all_matches.append(match_info)
            else:
                print(f"Skipping row with different transcript or non-matching codon.")
        
        if found_any:
            print("\n--- Variant Evaluation Results ---\n")
            for match in all_matches:
                print(f"Match Type: {match['match_type'].upper()}")
                print(f"Transcript: {match['transcript_id']}")
                print(f"Variant: {match['hgvs_c']} / {match['hgvs_p']}")
                print(f"Submitter: {match['submitter']}")
                print(f"Significance: {match['clinical_significance']}")
                print(f"Review Status: {match['review_status']}")
                print(f"ACMG Evidence: {match['evidence']}")
                
            # Optionally export.
            export_results(all_matches)
            
        else:
            print("No matching variant found in the database.")
            
# --- Export Results ---
def export_results(all_matches):
    # Export to TSV.
    with open(output_filename, "w", encoding = "utf-8", newline = "") as outfile:
        fieldnames = ["transcript_id", "hgvs_c", "hgvs_p", "submitter", "clinical_significance:", "review_status", "evidence", "match_type"]
        writer = csv.DictWriter(outfile, fieldnames = fieldnames, delimimiter = "\t")
        writer.writeheader()
        for row in all_matches:
            writer.write(row)
            
    print(f"\nResults saved to TSV: {output_filename}")
    
    # Export to JSON.
    json_filename = output_filename.replace(".tsv", ".json")
    with open(json_filename, "w", encoding = "utf-8") as jsonfile:
        json.dump(all_matches, jsonfile, indent = 4, ensure_ascii = False)
        
    print(f"Results also saved to JSON: {json_filename}")

found = False
criteria = []

# Open and read variant database.
with open(db, "r", encoding = "utf-8") as infile:
    reader = csv.DictReader(infile, delimiter = "\t")

    for row in reader:
        review_status = row.get("review_status", "").strip().lower()
        clinical_significance = row.get("clinical_significance","").strip().lower()
        hgvs_c = row.get("hgvs_c", "").strip()
        hgvs_p = row.get("hgvs_p", "").strip()
        transcript_id = row.get("transcript_id", "").strip()
        submitter = row.get("submitter", "").strip()

        # Extract codon from this row's protein variant.
        db_codon_match = re.search(r"[A-Za-z]{3}(\d+)[A-Za-z]{3}", hgvs_p)
        db_codon = db_codon_match.group(1) if db_codon_match else None

        

# Final output.
if found:
    print(f"\nVariant: {hgvs_c_input} / {hgvs_p_input}\n")
    print(f"Transcript: {transcript_id_input}\n")
    print("ACMG Evidence: ", ",".join(criteria))
else:
    print("No matching variant found in the database.")