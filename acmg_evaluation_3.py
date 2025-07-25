"""Enhanced Evalutaion of Variants with ACMG criteria (Refractored)"""
import csv
import re
import os
import json

# --- Configuration ---
gene = "gene_symbol"
db_filename = "gene_variants.tsv"
output_filename = f"{gene}_evaluation.tsv"

# Set of reputable submitters.
REPUTABLE_SUBMITTERS = {
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
def evaluate_variant(row, inputs, input_codon):
    transcript_id_input, hgvs_c_input, hgvs_p_input = inputs

    # Extract fields from row.
    transcript_id = row.get("transcript_id", "").strip()
    hgvs_c = row.get("hgvs_c", "").strip()
    hgvs_p = row.get("hgvs_p", "").strip()
    significance = row.get("clinical_significance", "").strip().lower()
    review_status = row.get("review_status", "").strip().lower()
    submitter = row.get("submitter", "").strip()
    db_codon = extract_codon(hgvs_p)

    exact_protein = hgvs_p_input == hgvs_p
    same_transcript = transcript_id_input == transcript_id
    same_codon = input_codon == db_codon and hgvs_p_input != hgvs_p
    is_pathogenic = "pathogenic" in significance
    is_benign = "benign" in significance
    is_vus = "uncertain significance" in significance
    is_conflict = "conflicting classifications" in significance
    is_reputable = submitter in REPUTABLE_SUBMITTERS
    is_single = review_status == "criteria provided, single submitter"

    criteria = []
    match_type = None # for tracking: exact, codon_match, skip

    if same_transcript and exact_protein:
        match_type = "exact protein"
        if is_pathogenic and is_single and is_reputable:
            criteria.append("PP5")
        elif is_benign and is_single and is_reputable:
            criteria.append("BP6")
        elif is_pathogenic and hgvs_c_input != hgvs_c:
            criteria.append("PS1")
        elif is_vus:
            criteria.append("VUS")
        elif is_conflict:
            criteria.append("Conflicting classifications")
        else:
            criteria.append("Unclassified variant")
    elif same_codon and is_pathogenic:
        match_type = "codon match"
        criteria.append("PM5")

return match_type, criteria, {
    "transcript_id": transcript_id,
    "hgvs_c": hgvs_c,
    "hgvs_p": hgvs_p,
    "submitter": submitter,
    "clinical_significance": significance,
    "review_status": review_status
}

# --- Export Results ---
def export_results(all_matches):
    # Export to TSV.
    with open(output_filename, "w", encoding = "utf-8", newline = "") as outfile:
        fieldnames = ["transcript_id", "hgvs_c", "hgvs_p", "submitter", "significance", "review_status", "match_type", "criteria"]
        writer = csv.DictWriter(outfile, fieldnames = fieldnames, delimiter = "\t")
        writer.writeheader()
        for row in all_matches:
            writer.writerow(row)
        print(f"\nResults saned to TSV: {output_filename}")

    # Export to JSON.
    json_filename = output_filename.replace(".tsv", ".json")
    with open(json_filename, "w", encoding = "utf-8") as outfile:
        json.dump(all_matches, outfile, indent = 4, ensure_ascii = False)
        print(f"\nResults also saved to JSON: {json_filename}")

# --- Main Evaluation ---
def main():
    if not os.path.exists(db_filename):
        print(f"Error: Database file '{db_filename}' not found.")
        return

    inputs = collect_input()
    input_codon = extract_codon(inputs[2])

    all_matches = []

    with open(db_filename, "r", encoding = "utf-8") as infile:
        reader = csv.DictReader(infile, delimiter = "\t")
        for row in reader:
            match_type, criteria, info = evaluate_variant(row, inputs, input_codon)
            if match_type:
                info.update({
                    "match_type": match_type,
                    "criteria": ",".join(criteria)
                })
                all_matches.append(info)

    if all_matches:
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

if __name__ == "__main__":
    main()
