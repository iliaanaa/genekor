"""Enhanced Evalutaion of Variants with ACMG criteria (Refractored and Optimized)"""
import csv
import re
import os
import json

# Configuration.
GENE = "gene_symbol"
DB_FILE= "gene_variants.tsv"
OUTPUT_FILE = f"{gene}_evaluation.tsv"

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
CRITERIA = {
    "PS1": "PS1", # Pathogenic, same amino acid, different nucleotide.
    "PM5": "PM5", # Pathogenic, same codon, different AA change.
    "PP5": "PP5", # Pathogenic, single reputable submitter.
    "BP6": "BP6", # Benign, single reputable submitter.
    "VUS": "Uncertain significance", # Used when significance is uncertain.
    "CONFLICT": "Conflicting classifications", # Indicates internal conflict.
    "UNCLASSIFIED": "Unclassified variant" # Catch-all for unknown significance.
}

# Input Collection.
def collect_input():
    return tuple(input(f"{label}: ").strip() for label in ["Transcript_ID", "HGVS_c", "HGVS_p"])

# Codon Extraction.
def extract_codon(hgvs_p_str):
    match = re.search(r"[A-Za-z]{3}(\d+)[A-Za-z]{3}", hgvs_p_str or "")
    return match.group(1) if match else None

# Criteria Evaluation.
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
    codon_match = input_codon == db_codon and hgvs_p_input != hgvs_p
    
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
            criteria.append(CRITERIA["PP5"])
        elif is_benign and is_single and is_reputable:
            criteria.append(CRITERIA["BP6"])
        elif is_pathogenic and hgvs_c_input != hgvs_c:
            criteria.append(CRITERIA["PS1"])
        elif is_vus:
            criteria.append(CRITERIA["VUS"])
        elif is_conflict:
            criteria.append(CRITERIA["CONFLICT"])
        else:
            criteria.append(CRITERIA["UNCLASSIFIED"])
    elif codon_match and is_pathogenic:
        match_type = "codon match"
        criteria.append(CRITERIA["PM5"])

    result = {
        "transcript_id": transcript_id,
        "hgvs_c": hgvs_c,
        "hgvs_p": hgvs_p,
        "submitter": submitter,
        "clinical_significance": significance,
        "review_status": review_status
    }

    return match_type, criteria, result

# Export Results.
def export_results(all_matches):
    fieldnames = [
        "transcript_id", "hgvs_c", "hgvs_p", "submitter",
        "clinical_significance", "review_status", "conflicting_interpetations", "match_type", "criteria"
    ]
    
    # Export to TSV.
    with open(output_filename, "w", encoding = "utf-8", newline = "") as outfile:
        writer = csv.DictWriter(outfile, fieldnames = fieldnames, delimiter = "\t")
        writer.writeheader()
        writer.writerow(matches)
    print(f"\nResults saved to TSV: {OUTPUT_FILE}")

    # Export to JSON.
    json_path = OUTPUT_FILE.replace(".tsv", ".json")
    with open(json_path, "w", encoding = "utf-8") as f_json:
        json.dump(all_matches, outfile, indent = 4, ensure_ascii = False)
    print(f"\nResults also saved to JSON: {json_filename}")

# Main.
def main():
    if not os.path.exists(DB_FILE):
        print(f"Error: Database file '{DB_FILE}' not found.")
        return

    inputs = collect_input()
    input_codon = extract_codon(inputs[2])
    all_matches = []

    with open(DB_FILE, "r", encoding = "utf-8") as infile:
        reader = csv.DictReader(infile, delimiter = "\t")
        for row in reader:
            match_type, criteria, info = evaluate_variant(row, inputs, input_codon)
            if match_type:
                result.update({
                    "match_type": match_type,
                    "criteria": ",".join(criteria)
                })
                all_matches.append(result)

    if all_matches:
        print("\nVariant Evaluation Results:\n")
        for match in all_matches:
            print(f"Match Type: {entry['match_type']}")
            print(f"Transcript: {entry['transcript_id']}")
            print(f"Variant: {entry['hgvs_c']} / {match['hgvs_p']}")
            print(f"Submitter: {entry['submitter']}")
            print(f"Significance: {entry['clinical_significance']}")
            print(f"Review Status: {entry['review_status']}")
            print(f"ACMG Evidence: {entry['evidence']}")

        # Optionally export.
        export_results(all_matches)
    else:
        print("No matching variants found.")

if __name__ == "__main__":
    main()
