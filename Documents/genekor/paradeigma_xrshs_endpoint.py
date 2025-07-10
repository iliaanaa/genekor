
#Αίτημα 1 (DNA Change):
curl -X POST "http://localhost:8000/get_variant_info/" \
-H "Content-Type: application/json" \
-d '{"gene": "BRCA1", "variant": "c.123G>T"}'

#Απάντηση:
{
  "gene": "BRCA1",
  "variant": "c.123G>T",
  "clinical_significance": "Pathogenic",
  "review_status": "Expert panel",
  "acmg_criteria": [
    {
      "criterion": "PS1",
      "explanation": "Same amino acid change as a known pathogenic variant (different nucleotide change)."
    },
    {
      "criterion": "PP5",
      "explanation": "Pathogenic classification from reputable source without published evidence."
    }
  ],
  "conflicting_interpretations": null
}



#Αίτημα 2 (Protein Change):
curl -X POST "http://localhost:8000/get_variant_info/" \
-H "Content-Type: application/json" \
-d '{"gene": "BRCA1", "variant": "p.Val12Cys"}'

#Απάντηση:
{
  "gene": "BRCA1",
  "variant": "p.Val12Cys",
  "clinical_significance": "Likely pathogenic",
  "review_status": "Conflicting interpretations",
  "acmg_criteria": [
    {
      "criterion": "PM5",
      "explanation": "Different amino acid change at the same position as a known pathogenic variant."
    }
  ],
  "conflicting_interpretations": [
    {
      "submitter": "LabA",
      "classification": "Pathogenic"
    },
    {
      "submitter": "LabB",
      "classification": "Benign"
    }
  ]
}