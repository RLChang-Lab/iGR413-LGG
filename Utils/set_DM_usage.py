import cobra

def set_dm(model, dm):
    """
    Sets DM series formulations to cobra model in cobrapy. 

    Args:
        model (Cobra model): base CB model to constrain
        dm (str): media series formulation to set. Options: '52', '24', '16', 13' 
    Returns: 
        model (Cobra model): media-constrained CB model 
        media (dict): dictionary of exchanges set. Keys= cobra reaction object for each reaction in DM. Values: LB corresponding to concentration

    Notes:
        - All formulations are set in the same way: 
            1. Media conditions are defined based on formulation selected in arg DM 
            2. All model exchange reactions are closed to LB=0 
            3. Iterative and controlled opening of media-specific exchanges based on formulation selected in arg DM 
        -Returns model, media as tuple 
        -Echos DM formulation selected when a proper argument is passed. Also provides ValueError if an invalid DM str is passed. 
        - import: from set_DMs import set_dm
        - function namespace is based on iGR413 LGG model 
        - Concentrations of DM components are relative uptake rates, glucose = -10 
        - Cobrapy Version: 0.29.1, Python Version: 3.9.12
    """
    valid_dms = ['52', '24', '16', '13']
    if dm not in valid_dms:
        raise ValueError(f"Invalid DM '{dm}' passed. Must be one of {valid_dms}.")
    
    if dm == '52': 
        #Step 1: Media conditions are defined based on formulation selected in arg DM 
        media={model.exchanges.EX_4abut_e: -0.238239650186896,
model.exchanges.EX_adn_e: -0.909023633807176,
model.exchanges.EX_ala_L_e: -3.43941818181818,
model.exchanges.EX_arg_L_e: -0.282058240267196,
model.exchanges.EX_asn_L_e: -0.0929733300305507,
model.exchanges.EX_asp_L_e: -0.55373266853357,
model.exchanges.EX_btn_e: -0.502788930606048,
model.exchanges.EX_C00072_e: -0.697458344517168,
model.exchanges.EX_cbl1_e: -0.00906293953948838,
model.exchanges.EX_cit_e: -2.17226868802977,
model.exchanges.EX_cobalt2_e: -0.710036784025223,
model.exchanges.EX_csn_e: -0.44225513460437,
model.exchanges.EX_cu2_e: -0.000983950365558824,
model.exchanges.EX_cys_L_e: -1.16893013200794,
model.exchanges.EX_cytd_e: -0.101009685701545,
model.exchanges.EX_fe2_e: -0.0323445102064021,
model.exchanges.EX_fol_e: -0.000556576183218685,
model.exchanges.EX_glc_D_e: -10,
model.exchanges.EX_glu_L_e: -0.500929913558201,
model.exchanges.EX_gly_e: -0.837744,
model.exchanges.EX_gua_e: -0.487671661363185,
model.exchanges.EX_his_L_e: -0.65449154972054,
model.exchanges.EX_ile_L_e: -0.468006545454545,
model.exchanges.EX_inost_e: -0.00681833698681884,
model.exchanges.EX_k_e: -19.7463439075564,
model.exchanges.EX_leu_L_e: -0.468006545454545,
model.exchanges.EX_lys_L_e: -0.881965090909091,
model.exchanges.EX_met_L_e: -0.164648969420768,
model.exchanges.EX_mg2_e: -5.10261883074804,
model.exchanges.EX_mn2_e: -0.0221638272953635,
model.exchanges.EX_mops_e: -9.82706405984726,
model.exchanges.EX_na1_e: -0.0840769087175658,
model.exchanges.EX_nac_e: -0.0402313481163887,
model.exchanges.EX_NH4_e: -7.08581509092581,
model.exchanges.EX_phe_L_e: -0.371702836363636,
model.exchanges.EX_pnto_R_e: -0.00515545143585351,
model.exchanges.EX_pro_L_e: -0.426767062628509,
model.exchanges.EX_pydx_e: -0.0120652552437249,
model.exchanges.EX_ribflv_e: -0.0130548864158699,
model.exchanges.EX_ser_L_e: -1.81153897697883,
model.exchanges.EX_thm_e: -0.185161838463014,
model.exchanges.EX_thr_L_e: -0.515417920291105,
model.exchanges.EX_thym_e: -0.0584424852762019,
model.exchanges.EX_trp_L_e: -0.336817309090909,
model.exchanges.EX_tyr_L_e: -0.271098358335077,
model.exchanges.EX_ura_e: -0.00438361724922243,
model.exchanges.EX_val_L_e: -0.209711414000006,
model.exchanges.EX_xan_e: -0.0193811894502184,
model.exchanges.EX_zn2_e: -0.00854216715134657,
model.exchanges.EX_pi_e: -19.7463439075564,
model.exchanges.EX_h_e: -0.0120652552437249,
model.exchanges.EX_hco3_e:-1}
        #Step 2: All model exchange reactions are closed to LB=0 
        for q in model.reactions: 
            if "EX_" in q.id: 
                q.lower_bound=0
        #Step 3: Iterative and controlled opening of media-specific exchanges
        for name, concentration in media.items(): 
            name.lower_bound=concentration
        print("dm52 set")
        return model,media
    
    elif dm == '24': 
        #Step 1: Media conditions are defined based on formulation selected in arg DM 
        media={
model.exchanges.EX_arg_L_e: -0.282058240267196,
model.exchanges.EX_asn_L_e: -0.0929733300305507,
model.exchanges.EX_asp_L_e: -0.55373266853357,
model.exchanges.EX_C00072_e: -0.697458344517168,
model.exchanges.EX_cit_e: -2.17226868802977,
model.exchanges.EX_cobalt2_e: -0.710036784025223,
model.exchanges.EX_cys_L_e: -1.16893013200794,
model.exchanges.EX_cytd_e: -0.101009685701545,
model.exchanges.EX_glc_D_e: -10,
model.exchanges.EX_glu_L_e: -0.500929913558201,
model.exchanges.EX_his_L_e: -0.65449154972054,
model.exchanges.EX_ile_L_e: -0.468006545454545,
model.exchanges.EX_k_e: -19.7463439075564,
model.exchanges.EX_met_L_e: -0.164648969420768,
model.exchanges.EX_mg2_e: -5.10261883074804,
model.exchanges.EX_mops_e: -9.82706405984726,
model.exchanges.EX_NH4_e: -7.08581509092581,
model.exchanges.EX_pro_L_e: -0.426767062628509,
model.exchanges.EX_ser_L_e: -1.81153897697883,
model.exchanges.EX_ura_e: -0.00438361724922243,
model.exchanges.EX_val_L_e: -0.209711414000006,
model.exchanges.EX_xan_e: -0.0193811894502184,
model.exchanges.EX_pi_e: -19.7463439075564,
model.exchanges.EX_h_e: -0.0120652552437249,
model.exchanges.EX_hco3_e:-1}
        #Step 2: All model exchange reactions are closed to LB=0 
        for q in model.reactions: 
            if "EX_" in q.id: 
                q.lower_bound=0
        #Step 3: Iterative and controlled opening of media-specific exchanges
        for name, concentration in media.items(): 
            name.lower_bound=concentration
        print("dm24 set")
        return model,media
    
    elif dm == '16':
        #Step 1: Media conditions are defined based on formulation selected in arg DM 
        media={
model.exchanges.EX_arg_L_e: -0.282058240267196,
model.exchanges.EX_asn_L_e: -0.0929733300305507,
model.exchanges.EX_asp_L_e: -0.55373266853357,
model.exchanges.EX_cobalt2_e: -0.710036784025223,
model.exchanges.EX_cys_L_e: -1.16893013200794,
model.exchanges.EX_cytd_e: -0.101009685701545,
model.exchanges.EX_glc_D_e: -10,
model.exchanges.EX_glu_L_e: -0.500929913558201,
model.exchanges.EX_his_L_e: -0.65449154972054,
model.exchanges.EX_ile_L_e: -0.468006545454545,
model.exchanges.EX_k_e: -19.7463439075564,
model.exchanges.EX_mg2_e: -5.10261883074804,
model.exchanges.EX_mops_e: -9.82706405984726,
model.exchanges.EX_NH4_e: -7.08581509092581,
model.exchanges.EX_pro_L_e: -0.426767062628509,
model.exchanges.EX_val_L_e: -0.209711414000006,
model.exchanges.EX_pi_e: -19.7463439075564,
model.exchanges.EX_h_e: -0.0120652552437249,
model.exchanges.EX_hco3_e:-1}
        #Step 2: All model exchange reactions are closed to LB=0 
        for q in model.reactions: 
            if "EX_" in q.id: 
                q.lower_bound=0
        #Step 3: Iterative and controlled opening of media-specific exchanges
        for name, concentration in media.items(): 
            name.lower_bound=concentration
        print("dm16 set")
        return model,media

    elif dm == '13':
        #Step 1: Media conditions are defined based on formulation selected in arg DM 
        media={
model.exchanges.EX_arg_L_e: -0.282058240267196,
model.exchanges.EX_asn_L_e: -0.0929733300305507,
model.exchanges.EX_asp_L_e: -0.55373266853357,
model.exchanges.EX_cobalt2_e: -0.710036784025223,
model.exchanges.EX_cys_L_e: -1.16893013200794,
model.exchanges.EX_glc_D_e: -10,
model.exchanges.EX_glu_L_e: -0.500929913558201,
model.exchanges.EX_his_L_e: -0.65449154972054,
model.exchanges.EX_ile_L_e: -0.468006545454545,
model.exchanges.EX_k_e: -19.7463439075564,
model.exchanges.EX_mg2_e: -5.10261883074804,
model.exchanges.EX_pro_L_e: -0.426767062628509,
model.exchanges.EX_val_L_e: -0.209711414000006,
model.exchanges.EX_pi_e: -19.7463439075564,
model.exchanges.EX_h_e: -0.0120652552437249,
model.exchanges.EX_hco3_e:-1}
        #Step 2: All model exchange reactions are closed to LB=0 
        for q in model.reactions: 
            if "EX_" in q.id: 
                q.lower_bound=0
        #Step 3: Iterative and controlled opening of media-specific exchanges
        for name, concentration in media.items(): 
            name.lower_bound=concentration
        print("dm13 set")
        return model,media
