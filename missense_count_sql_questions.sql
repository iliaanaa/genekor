--Πόσες Missense παραλλαγές έχει στην clinvar για το BRCA1
SELECT COUNT(*) AS missense_count
FROM gene_variants
WHERE gene_symbol = 'TP53'
    AND molecular_consequence = 'missense';

--Πόσες Missense παραλλαγές έχει στην clinvar για το BRCA1 που είναι παθογόνες
SELECT COUNT(*) AS pathogenic_missense_count
FROM gene_variants
WHERE gene_symbol = 'TP53'
    AND molecular_consequence = 'missense'
    AND clinical_significance = 'Pathogenic';

--Πόσες Missense παραλλαγές έχει στην clinvar για το BRCA1 που είναι παθογόνες και βρίσκονται μεταξύ των αμινοξέων 1756-1855
SELECT COUNT(*) AS pathogenic_missense_count_in_range
FROM gene_variants
WHERE gene_symbol  = 'TP53'
AND molecular_consequence = 'missense'
AND clinical_significance = 'Pathogenic'
AND protein_pos BETWEEN 1756 AND 1855;

--Πόσες Missense παραλλαγές έχει στην clinvar για το BRCA1 στην θέση 34 της πρωτεΐνης
SELECT COUNT(*) AS missense_count_at_position
FROM gene_variants
WHERE gene_symbol = 'TP53'
    AND molecular_consequence = 'missense'
    AND protein_pos = 34;