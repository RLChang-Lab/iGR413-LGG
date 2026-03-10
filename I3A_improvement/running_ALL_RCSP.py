#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 10:24:40 2025

@author: grichmond
"""

import cobra
import os 
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from run_fba import fba
ATPmodel = cobra.io.read_sbml_model('iGR_ATPmain.xml') #iGR413 with ATPM as objective function
from set_dms import set_dm
I3Amodel = cobra.io.read_sbml_model('iGR413.xml')

def UB_locked_sim(model, bm_val, objective_rxn, biomass_rxn, lac_check=0): 
    results={}
    model.reactions.get_by_id(biomass_rxn).bounds= bm_val, bm_val
    x=flux_variability_analysis(model)
    reduced_model=model.copy()
    to_remove=[]
    for rxn, results in x.iterrows():
        if results[0]<=1e-6 and results[1]<=1e-6:  
            to_remove.append(rxn)
    if lac_check == 1:
        if "LDH_L" not in to_remove: 
            print("lactate had non zero FVA flux")
            model.reactions.LDH_L.bounds=0,0
        else:
            print("lactate removed in FVA step") 
    to_remove = [reduced_model.reactions.get_by_id(rxn_id) for rxn_id in to_remove]
    reduced_model.remove_reactions(to_remove, remove_orphans=True) 
    y=fba(reduced_model, objective=objective_rxn)
    return y.shadow_prices, y.reduced_costs
def UB_LB_locked_sim(model, bm_val, objective_rxn, biomass_rxn, lac_check=0):
    results={}
    model.reactions.get_by_id(biomass_rxn).bounds= 0, bm_val
    x=flux_variability_analysis(model)
    reduced_model=model.copy()
    to_remove=[]
    for rxn, results in x.iterrows():
        if results[0]<=1e-6 and results[1]<=1e-6:  
            to_remove.append(rxn)
    if lac_check == 1:
        if "LDH_L" not in to_remove: 
            print("lactate had non zero FVA flux")
            model.reactions.LDH_L.bounds=0,0
        else:
            print("lactate removed in FVA step") 
    to_remove = [reduced_model.reactions.get_by_id(rxn_id) for rxn_id in to_remove]
    reduced_model.remove_reactions(to_remove, remove_orphans=True) 
    y=fba(reduced_model, objective=objective_rxn)
    return y.shadow_prices, y.reduced_costs
def BM_open_sim(model, objective_rxn, lac_check=0):
    results={}
    x=flux_variability_analysis(model)
    reduced_model=model.copy()
    to_remove=[]
    for rxn, results in x.iterrows():
        if results[0]<=1e-6 and results[1]<=1e-6:  
            to_remove.append(rxn)
    if lac_check == 1:
        if "LDH_L" not in to_remove: 
            print("lactate had non zero FVA flux")
            model.reactions.LDH_L.bounds=0,0
        else:
            print("lactate removed in FVA step") 
    to_remove = [reduced_model.reactions.get_by_id(rxn_id) for rxn_id in to_remove]
    reduced_model.remove_reactions(to_remove, remove_orphans=True) 
    y=fba(reduced_model, objective=objective_rxn)
    return y.shadow_prices, y.reduced_costs
def phase_minMax_sim(model, min_BM_val, max_BM_val, biomass_rxn, objective_rxn, lac_check=0):
    results={}
    model.reactions.get_by_id(biomass_rxn).bounds= min_BM_val, max_BM_val
    x=flux_variability_analysis(model)
    reduced_model=model.copy()
    to_remove=[]
    for rxn, results in x.iterrows():
        if results[0]<=1e-6 and results[1]<=1e-6:  
            to_remove.append(rxn)
    if lac_check == 1:
        if "LDH_L" not in to_remove: 
            print("lactate had non zero FVA flux")
            model.reactions.LDH_L.bounds=0,0
        else:
            print("lactate removed in FVA step") 
    to_remove = [reduced_model.reactions.get_by_id(rxn_id) for rxn_id in to_remove]
    reduced_model.remove_reactions(to_remove, remove_orphans=True) 
    y=fba(reduced_model, objective=objective_rxn)
    return y.shadow_prices, y.reduced_costs

dm24_i3a_bmVals=[0.15]
dm24_i3a_bm_Range={0.00:0.24}
dm52_i3a_bmVals=[0.143, 0.381]
dm52_i3a_bm_Range={0.00:0.237, 0.332:0.437}

dm24_lac_bmVals=[0.14, 0.30]
dm24_lac_bm_Range={0.0696:0.209, 0.278:0.313}
dm52_lac_bmVals=[0.333]
dm52_lac_bm_Range={0.237:0.427}

I3A_tmodelc24, media24 = set_dm(I3Amodel, '24')
I3A_tmodelc52, media52 = set_dm(I3Amodel, '52')
lac_tmodelc24, media24 = set_dm(ATPmodel, '24')
lac_tmodelc52, media52 = set_dm(ATPmodel, '52')

i3a_24_open_SP, i3a_24_open_RC=BM_open_sim(I3A_tmodelc24, 'IAA_I3A', lac_check=0)
i3a_52_open_SP, i3a_52_open_RC=BM_open_sim(I3A_tmodelc52, 'IAA_I3A', lac_check=0)
lac_24_open_SP, lac_24_open_RC=BM_open_sim(lac_tmodelc24, 'ATPM', lac_check=1)
lac_52_open_SP, lac_52_open_RC=BM_open_sim(lac_tmodelc52, 'ATPM', lac_check=1)

UB_I3A_24_results={}
UBLB_I3A_24_results={}
for i in dm24_i3a_bmVals: 
    i3a_24_UB_SP, i3a_24_UB_RC=UB_locked_sim(I3A_tmodelc24, i, 'IAA_I3A', 'curated_biomass', lac_check=0)
    i3a_24_UBLB_SP, i3a_24_UBLB_RC=UB_LB_locked_sim(I3A_tmodelc24, i, 'IAA_I3A', 'curated_biomass', lac_check=0)
    UB_I3A_24_results[i]=i3a_24_UB_SP, i3a_24_UB_RC
    UBLB_I3A_24_results[i]=i3a_24_UBLB_SP, i3a_24_UBLB_RC
UB_lac_24_results={}
UBLB_lac_24_results={}
for i in dm24_lac_bmVals: 
    lac_24_UB_SP, lac_24_UB_RC=UB_locked_sim(lac_tmodelc24, i, 'ATPM', 'BM_noATP', lac_check=1)
    lac_24_UBLB_SP, lac_24_UBLB_RC=UB_LB_locked_sim(lac_tmodelc24, i, 'ATPM', 'BM_noATP', lac_check=1)
    UB_lac_24_results[i]=lac_24_UB_SP, lac_24_UB_RC
    UBLB_lac_24_results[i]=lac_24_UBLB_SP, lac_24_UBLB_RC

UB_I3A_52_results={}
UBLB_I3A_52_resuls={}
for i in dm52_i3a_bmVals: 
    i3a_52_UB_SP, i3a_52_UB_RC=UB_locked_sim(I3A_tmodelc52, i, 'IAA_I3A', 'curated_biomass', lac_check=0)
    i3a_52_UBLB_SP, i3a_52_UBLB_RC=UB_LB_locked_sim(I3A_tmodelc52, i, 'IAA_I3A', 'curated_biomass', lac_check=0)
    UB_I3A_52_results[i]=i3a_52_UB_SP, i3a_52_UB_RC
    UBLB_I3A_52_resuls[i]=i3a_52_UBLB_SP, i3a_52_UBLB_RC
UB_lac_52_results={}
UBLB_lac_52_results={}
for i in dm52_lac_bmVals: 
    lac_52_UB_SP, lac_52_UB_RC=UB_locked_sim(lac_tmodelc52, i, 'ATPM', 'BM_noATP', lac_check=1)
    lac_52_UBLB_SP, lac_52_UBLB_RC=UB_LB_locked_sim(lac_tmodelc52, i, 'ATPM', 'BM_noATP', lac_check=1)
    UB_lac_24_results[i]=lac_52_UB_SP, lac_52_UB_RC
    UBLB_lac_24_results[i]=lac_52_UBLB_SP, lac_52_UBLB_RC

phase_24_i3a_results={}
for bm_min, bm_max in dm24_i3a_bm_Range.items():
    phase_I3A_24_SP, phase_I3A_24_RC=phase_minMax_sim(I3A_tmodelc24, bm_min, bm_max, 'curated_biomass', 'IAA_I3A', lac_check=0)
    phase_24_i3a_results[str(str(bm_min) + "-" + str(bm_max))]= phase_I3A_24_SP, phase_I3A_24_RC

phase_52_i3a_results={}
for bm_min, bm_max in dm52_i3a_bm_Range.items():
    phase_I3A_52_SP, phase_I3A_52_RC=phase_minMax_sim(I3A_tmodelc52, bm_min, bm_max, 'curated_biomass', 'IAA_I3A', lac_check=0)
    phase_52_i3a_results[str(str(bm_min) + "-" + str(bm_max))]= phase_I3A_52_SP, phase_I3A_52_RC

phase_24_lac_results={}
for bm_min, bm_max in dm24_lac_bm_Range.items():
    phase_lac_24_SP, phase_lac_24_RC=phase_minMax_sim(lac_tmodelc24, bm_min, bm_max, 'BM_noATP', 'ATPM', lac_check=1)
    phase_24_lac_results[str(str(bm_min) + "-" + str(bm_max))]= phase_I3A_24_SP, phase_I3A_24_RC

phase_52_lac_results={}
for bm_min, bm_max in dm52_lac_bm_Range.items():
    phase_lac_52_SP, phase_lac_52_RC=phase_minMax_sim(lac_tmodelc52, bm_min, bm_max, 'BM_noATP', 'ATPM', lac_check=1)
    phase_52_lac_results[str(str(bm_min) + "-" + str(bm_max))]= phase_lac_52_SP, phase_lac_52_RC


def save_to_excel(filename, data_dict, kind="pair"):
    """
    Save shadow prices (SP) and reduced costs (RC) from result dicts to one Excel workbook.
    
    Parameters
    ----------
    filename : str
        Output Excel file name (no extension needed).
    data_dict : dict
        Dictionary where key -> BM_val or range, and value -> tuple of (SP_df, RC_df)
    kind : str
        'open', 'pair', or 'phase'
        'open' means flat outputs with 24_SP etc.
        'pair' means each BM_val or range gets two sheets (SP and RC)
    """
    with pd.ExcelWriter(f"{filename}.xlsx", engine="openpyxl") as writer:
        if kind == "open":
            # expects keys like "24_SP", "24_RC"
            for sheet_name, df in data_dict.items():
                df.to_excel(writer, sheet_name=str(sheet_name)[:31])
        else:
            # expects dict: {BM_value: (SP_df, RC_df)}
            for key, (SP, RC) in data_dict.items():
                sheet_base = str(key).replace(":", "-").replace(".", "_")
                SP.to_excel(writer, sheet_name=f"{sheet_base}_SP"[:31])
                RC.to_excel(writer, sheet_name=f"{sheet_base}_RC"[:31])
    print(f"✅ Saved: {filename}.xlsx")

# ---- create open sets ----
i3a_open_sets = {"24_SP": i3a_24_open_SP,"24_RC": i3a_24_open_RC,"52_SP": i3a_52_open_SP,"52_RC": i3a_52_open_RC}
lac_open_sets = {"24_SP": lac_24_open_SP,"24_RC": lac_24_open_RC,"52_SP": lac_52_open_SP,"52_RC": lac_52_open_RC,}
os.chdir("/Users/grichmond/Desktop/RC_SPouts")
save_to_excel("I3A_open_results", i3a_open_sets, kind="open")
save_to_excel("LAC_open_results", lac_open_sets, kind="open")
save_to_excel("I3A_UB_results", UB_I3A_24_results, kind="pair")
save_to_excel("I3A_UB_52_results", UB_I3A_52_results, kind="pair")
save_to_excel("LAC_UB_results", UB_lac_24_results, kind="pair")
save_to_excel("LAC_UB_52_results", UB_lac_52_results, kind="pair")
save_to_excel("I3A_UBLB_24_results", UBLB_I3A_24_results, kind="pair")
save_to_excel("I3A_UBLB_52_results", UBLB_I3A_52_resuls, kind="pair")
save_to_excel("LAC_UBLB_24_results", UBLB_lac_24_results, kind="pair")
save_to_excel("LAC_UBLB_52_results", UBLB_lac_52_results, kind="pair")
save_to_excel("I3A_phase_24_results", phase_24_i3a_results, kind="pair")
save_to_excel("I3A_phase_52_results", phase_52_i3a_results, kind="pair")
save_to_excel("LAC_phase_24_results", phase_24_lac_results, kind="pair")
save_to_excel("LAC_phase_52_results", phase_52_lac_results, kind="pair")
