"""
Molecular consequence of variants.
DNA and protein included
Reads HGVS_c and HGVS_p and categorizes each variant by molecular consequence.
"""

import re

# Molecular consequence classifiers.
def consequence(hgvs_p):
    if not hgvs_p or hgvs_p.startswith('p.?'):
        return 'unknown'
    
    # Nonsense variant (premature stop codon, e.g., p.Arg127*)
    elif re.fullmatch(r'p\.[A-Z][a-z]{2}\d+(Ter|X|\*)', hgvs_p):
        return 'nonsense'

    # Stop lost (original stop codon changed to something else, e.g., p.Ter127Gln or p.X127Q)
    elif re.search(r'p\.(Ter|X|\*)\d+[A-Z][a-z]{2}', hgvs_p):
        return 'stop_lost'

    # Missense (e.g., p.Arg127Gly)
    elif re.fullmatch(r'p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}', hgvs_p):
        return 'missense_variant'
    
    # Synonymous (no amino acid change, e.g. p.Arg127=)
    elif '=' in hgvs_p:
        return 'synonymous_variant'

    # Frameshift (e.g., p.Arg127fs or p.Arg127AlafsTer5)
    elif re.search(r'fs', hgvs_p):
        return 'frameshift_variant'

    # Inframe deletion (no frameshift, with "del", e.g., p.Arg127del, p.Arg127_Lys130del)
    elif 'del' in hgvs_p and not re.search(r'fs', hgvs_p):
        return 'inframe_deletion'
    
    # Inframe insertion (no frameshift, with "ins", e.g., p.Arg127insGly, p.Arg127_Gly128insSer)
    elif 'ins' in hgvs_p and not re.search(r'fs', hgvs_p):
        return 'inframe_insertion'

    # Fallback
    return 'other'

def consequence_dna(hgvs_c):
    if not hgvs_c or hgvs_c.startswith('c.?'):
        return 'unknown'
    
    if re.search(r'[-+]\d+', hgvs_c):
        return 'intronic_variant'
    elif re.search(r'c\.\*\d+', hgvs_c):
        return '3_prime_UTR_variant'
    elif re.search(r'c\.-\d+', hgvs_c):
        return '5_prime_UTR_variant'
    elif re.search(r'=\s*$', hgvs_c):
        return 'synonymous_variant'
    elif 'del' in hgvs_c and 'ins' in hgvs_c:
        return 'indel'
    elif 'del' in hgvs_c:
        return 'deletion'
    elif 'ins' in hgvs_c:
        return 'insertion'
    elif 'dup' in hgvs_c:
        return 'duplication'
    elif '>' in hgvs_c:
        return 'substitution'
    
    return 'other'

def variant_consequence(hgvs_c, hgvs_p):
    return {
        'protein_consequence': consequence(hgvs_p),
        'dna_consequence': consequence_dna(hgvs_c)
    }

# Variant name parser.
def categorize_variant_name(name: str) -> dict:
    """
    Parses variant 'Name' field to extract HGVS_c and HGVS_p,
    then classifies moleucular consequence using variant_consequence().
    """

    result = {
        'molecular_consequence': None,
        'DNA_variant': None,
        'Protein_variant': None,
        'Other_variant': None
    }

    if not name or not isinstance(name, str):
        result['molecular_consequence'] = 'Other'
        result['Other_variant'] = name
        return result

    # Extract HGVS notation.
    dna_match = re.search(r'(c.[^()\s]+)', name)
    protein_match = re.search(r'(p\.[^()\s]+)', name)

    hgvs_c = dna_match.group(1) if dna_match else None
    hgvs_p = dna_match.group(1) if protein_match else None

    # Get detailed consequences.
    consequences = variant_consequence(hgvs_c, hgvs_p)

    # Fill result dictionary.
    if hgvs_p:
        result['DNA_variant'] = hgvs_c
        result['Protein_variant'] = hgvs_p
        result['molecular_consequence'] = consequences['protein_consequence']
    elif hgvs_c:
        result['DNA_variant'] = hgvs_c
        result['molecular_consequence'] = consequences['dna_consequence']
    else:
        result['molecular_consequence'] = 'Other'
        result['Other_variant'] = name

    return result
