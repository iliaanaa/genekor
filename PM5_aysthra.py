#PM5 aysthra
if gene and prot_pos and input_c_hgvs and input_p_hgvs:
    group = support_tables['same_pos_groups'].get(f"{gene}:{str(prot_pos)}", [])

    for v in group:
        clin_sign = v.get('clinicalsignificance', '').lower()
        v_c_hgvs = v.get('hgvs_c')
        v_hgvs_p = v.get('hgvs_p')
        consequence = v.get('consequence', '').lower()

        if(
            clin_sign in ['pathogenic', 'likely pathogenic'] and
            'missense' in consequence and
            v_hgvs_p != input_p_hgvs and
            v_c_hgvs != input_c_hgvs    
        ):
            criteria.append("PM5")
            break