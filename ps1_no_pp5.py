# PS1 — ίδιο p.HGVS, διαφορετικό c.HGVS από input
if gene and hgvs_p and input_c_hgvs:
    group = support_tables['same_p_groups'].get(f"{gene}:{hgvs_p}", [])
    for v in group:
        # Πρέπει να είναι pathogenic και διαφορετικό c.HGVS
        if (
            v.get('clinicalsignificance', '').lower() in ['pathogenic', 'likely pathogenic'] and
            v.get('hgvs_c') != input_c_hgvs
        ):
            criteria.append("PS1")
            break

# PS1 — ίδιο p.HGVS, διαφορετικό c.HGVS, variant είναι missense και (likely) pathogenic
if gene and hgvs_p and input_c_hgvs:
    group = support_tables['same_p_groups'].get(f"{gene}:{hgvs_p}", [])
    for v in group:
        clin_sign = v.get('clinicalsignificance', '').lower()
        v_c_hgvs = v.get('hgvs_c')
        consequence = v.get('consequence', '').lower()
        
        if (
            clin_sign in ['pathogenic', 'likely pathogenic'] and
            v_c_hgvs != input_c_hgvs and
            'missense' in consequence
        ):
            criteria.append("PS1")
            break

