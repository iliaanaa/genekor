# Λεξικά γνωστών pathogenic variants και θέσεων
known_pathogenic = {
    'p.Arg504Gly': 'c.1510A>T',  # PS1: Ίδιο protein change, διαφορετικό DNA
    'p.Trp41*': 'c.123G>A'       # PM5: Θέση 41 είναι hotspot
}
pathogenic_positions = [504, 41]  # PM5: Γνωστές pathogenic θέσεις
trusted_submitters = ['ENIGMA', 'ClinVar']  # PP5/BP6: Αξιόπιστες πηγές

# Συνάρτηση ελέγχου ACMG κριτηρίων
def apply_acmg_criteria(row):
    criteria = []
    
    # PS1
    if row['ProteinChange'] in known_pathogenic and row['DNAChange'] != known_pathogenic[row['ProteinChange']]:
        criteria.append('PS1')
    
    # PM5
    protein_pos = int(''.join(filter(str.isdigit, row['ProteinChange']))) if pd.notna(row['ProteinChange']) else None
    if protein_pos in pathogenic_positions and row['ProteinChange'] not in known_pathogenic:
        criteria.append('PM5')
    
    # PP5/BP6
    if row['Submitter'] in trusted_submitters:
        if row['ClinicalSignificance'] == 'Pathogenic':
            criteria.append('PP5')
        elif row['ClinicalSignificance'] == 'Benign':
            criteria.append('BP6')
    
    return criteria

# Εφαρμογή στη στήλη 'acmg_criteria'
df_brca['acmg_criteria'] = df_brca.apply(apply_acmg_criteria, axis=1)