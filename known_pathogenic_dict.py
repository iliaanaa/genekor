"""Create a dictionary incuding the clinical significance of all the variants"""
# Filter known pathogenic variants to create a dictionary.
# This is tailored based on the variables in download_and_parse_data_new.py.

import csv
import re 

# --- Known Pathogenic Dictionary ---
# Create a dictionary of known pathogenic variants (from BRCA1).
# This will include HGVS_c, HGVS_p, clinical significance and codon(amino acid position - if available).
gene = "brca1"
INPUT_FILE = "{gene}_variants.tsv"
known_pathogenic = {}

with open(INPUT_FILE, "r", encoding = "utf-8") as infile:
	reader = csv.DictReader(infile, delimiter = "\t")
	for row in reader:
		clinical_significance = row.get("clinical_significance", "")
		hgvs_c = row.get("hgvs_c", "").strip()
		hgvs_p = row.get("hgvs_p", "").strip()
	
		if "Pathogenic" in clinical_significance:
			key = hgvs_p if hgvs_p else hgvs_c if hgvs_c else None
        		codon = None
        
        		# Extract codon number if HGVS_p is available and matches typical format.
       			codon_match = re.search(r"[A-Za-z]{3}(\d+)[A-Za-z]{3}", hgvs_p)
        		if codon_match:
            			codon = codon_match.group(1)
            
			if key:
            			if key in known_pathogenic:
                			print(f"Warning: Duplicate key '{key}' encountered.")
				known_pathogenic[key] = {
                			"Clinical Significance": clinical_significance,
                			"HGVS_c": hgvs_c,
                			"HGVS_p": hgvs_p or "N/A"
				}
