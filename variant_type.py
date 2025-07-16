# Annotates the variants based on the consequences of their mutations.
# Adds the column "Variant_Type" to the final database.
import csv
import re

# Input_file is the file that is used for the database. Alternatively it can be used variant_summary.txt as well. The first option is better because it contains hgvs_p column.
# The output_file can have the same name as the input_file.
input_file = "brca1_variants.tsv" 
output_file = "brca1_variants_full.tsv"

with open(input_file, "r", encoding = "utf-8") as infile, \
    open(output_file, "w", newline = '', encoding = "utf-8") as outfile:
        
        reader = csv.DictReader(infile, delimiter = "\t")
        fieldnames = reader.fieldnames + ["Variant_Type"] # Add new column.
        writer = csv.DictWriter(outfile, fieldnames = fieldnames, delimiter = "\t")
        writer.writeheader()
        
        # Classify based on HGVS_p (protein-level changes).
        for row in reader:
            # Selects the columns hgvs_c and hgvs_p from brca1_variant.tsv.
            hgvs_p = row.get("hgvs_p", "").strip()
            hgvs_c = row.get("hgvs_c", "").strip()
            
            # Default value.
            variant_type = "unknown"
            if hgvs_p:
                if "fs" in hgvs_p:
                    variant_type = "frameshift"
                elif "*" in hgvs_p:
                    variant_type = "nonsense"
                elif "del" in hgvs_p:
                    variant_type = "deletion"
                elif "dup" in hgvs_p:
                    variant_type = "duplication"
                elif "ins" in hgvs_p:
                    variant_type = "insertion"
                elif hgvs_p == "p.=":
                    variant_type = "synonymous"
                elif re.match(r"p\.[A-Z][a-z]{2]\d+[A-Z][a-z]{2}", hgvs_p):
                    variant_type = "missense"
                else:
                    variant_type = "protein_other"
                    
            # If no HGVS_p, use HGVS_c for splice site or UTR.
            else:
                if re.search(r"\+\d+|\-\d+", hgvs_c):
                    if re.search(r"\+1\+2|\-1\-2", hgvs_c):
                        variant_type = "splice_site_essential"
                    else:
                        variant_type = "splice_region"
                elif hgvs_c.startswith("c.-"):
                    variant_type = "5'UTR"
                elif "*" in hgvs_c:
                    variant_type = "3'UTR"
                else:
                    variant_type = "non_coding"
                
            # Add new column.
            row["Variant_Type"] = variant_type
            writer.writerow(row)