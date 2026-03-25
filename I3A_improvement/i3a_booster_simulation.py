#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 14:17:51 2026

@author: grichmond
"""


import cobra
import numpy as np
import matplotlib.pyplot as plt
from cobra.flux_analysis import flux_variability_analysis
import os
from set_DM_usage import set_dm 
from run_fba_usage import fba


def pe_Data(model, target, objective='curated_biomass', num_points=10):
    # Validate objective
    if objective not in [rxn.id for rxn in model.reactions]:
        raise KeyError(f"Objective reaction '{objective}' not in model.")

    # Validate target
    if target not in [rxn.id for rxn in model.reactions]:
        raise KeyError(f"Target reaction '{target}' not in model.")

    # 1. Baseline FBA to get max biomass
    base_fba = fba(model, objective=objective)

    if base_fba.objective_value is None:
        raise ValueError(f"FBA returned None objective value for {objective}")

    max_bm = float(base_fba.objective_value)

    # 2. Evaluate production envelope
    bm_values = np.linspace(0, max_bm, num_points)
    results = {}

    for b in bm_values:
        mcopy = model.copy()

        # Constrain biomass to fixed point
        mcopy.reactions.get_by_id(objective).bounds = (b, b)

        try:
            sol = flux_variability_analysis(mcopy, processes=1)
        except Exception:
            # Infeasible or solver failure → assign zeros
            results[b] = {
                'min_target': 0,
                'max_target': 0
            }
            continue

        # Extract min and max for target
        try:
            target_min = sol.loc[target, 'minimum']
            target_max = sol.loc[target, 'maximum']
        except KeyError:
            target_min, target_max = 0, 0

        results[b] = {
            'min_target': target_min,
            'max_target': target_max
        }

    return results

def calculate_yields(PE_data): 
    result = {}
    
    for biomass, secretion in PE_data.items():
        max_flux = secretion["max_target"]
        
        if biomass == 0:
            point_yield = None  # or float('inf') if you prefer
        else:
            point_yield = max_flux / biomass
        
        result[biomass] = {
            "max_target": max_flux,
            "yield": point_yield
        }
    
    return result


model=cobra.io.read_sbml_model('iGR413.xml')
model_52, media_52= set_dm(model,"52")
unsup_52=pe_Data(model_52, "EX_i3a_e")
unsup_52_Yield=calculate_yields(unsup_52)

model_I= model_52.copy()
model_I.add_boundary(model_52.metabolites.indole_c, type="sink")
model_I.reactions.get_by_id('SK_indole_c').lower_bound= -0.125 #insert here, gotta check
I_data=pe_Data(model_I, "EX_i3a_e")
I_yields=calculate_yields(I_data)

model_Q= model_52.copy()
model_Q.add_boundary(model_52.metabolites.quln_c, type="sink")
model_Q.reactions.get_by_id('SK_quln_c').lower_bound= -0.025 #insert here, gotta check
Q_data=pe_Data(model_Q, "EX_i3a_e")
Q_yields=calculate_yields(Q_data)

model_R25= model_52.copy()
model_R25.add_boundary(model_52.metabolites.rib_D_c, type="sink")
model_R25.reactions.get_by_id('SK_rib_D_c').lower_bound= -2.5 #insert here, gotta check
R25_data=pe_Data(model_R25, "EX_i3a_e")
R25_yields=calculate_yields(R25_data)


model_R100= model_52.copy()
model_R100.add_boundary(model_52.metabolites.rib_D_c, type="sink")
model_R100.reactions.get_by_id('SK_rib_D_c').lower_bound= -10 #insert here, gotta check
model_R100.reactions.get_by_id("EX_glc_D_e").lower_bound=0
R100_data=pe_Data(model_R100, "EX_i3a_e")
R100_yields=calculate_yields(R100_data)

#######################
model_24, media_25= set_dm(model,"24")
unsup_24=pe_Data(model_24, "EX_i3a_e")
unsup_24_Yield=calculate_yields(unsup_24)

model_I_24= model_24.copy()
model_I_24.add_boundary(model_24.metabolites.indole_c, type="sink")
model_I_24.reactions.get_by_id('SK_indole_c').lower_bound= -0.125 #insert here, gotta check
I_data_24=pe_Data(model_I_24, "EX_i3a_e")
I_24_Yield=calculate_yields(I_data_24)


model_Q_24= model_24.copy()
model_Q_24.add_boundary(model_24.metabolites.quln_c, type="sink")
model_Q_24.reactions.get_by_id('SK_quln_c').lower_bound= -0.025 #insert here, gotta check
Q_data_24=pe_Data(model_Q_24, "EX_i3a_e")
Q_24_Yield=calculate_yields(Q_data_24)

model_R25_24= model_24.copy()
model_R25_24.add_boundary(model_24.metabolites.rib_D_c, type="sink")
model_R25_24.reactions.get_by_id('SK_rib_D_c').lower_bound= -2.5 #insert here, gotta check
R25_data_24=pe_Data(model_R25_24, "EX_i3a_e")
R25_24_Yield=calculate_yields(R25_data_24)



model_R100_24= model_24.copy()
model_R100_24.add_boundary(model_24.metabolites.rib_D_c, type="sink")
model_R100_24.reactions.get_by_id('SK_rib_D_c').lower_bound= -10 #insert here, gotta check
model_R100_24.reactions.get_by_id("EX_glc_D_e").lower_bound=0
R100_data_24=pe_Data(model_R100_24, "EX_i3a_e")
R100_24_Yield=calculate_yields(R100_data_24)

yield_dict={
    "52_unsup": unsup_52_Yield,
    "52_indole": I_yields,
    "52_quln": Q_yields,
    "52_Ribose_25": R25_yields,
    "52_Ribose_100": R100_yields,
    "24_unsup": unsup_24_Yield, 
    "24_indole": I_24_Yield, 
    "24_quln": Q_24_Yield, 
    "24_Ribose_25": R25_24_Yield,
    "24_Ribose_100": R100_24_Yield
}

import pandas as pd
def export_yields_to_excel(results_dict, filename="yields.xlsx"):
    with pd.ExcelWriter(filename, engine="openpyxl") as writer:
        
        for sheet_name, result in results_dict.items():
            
            # Convert nested dict → DataFrame
            df = pd.DataFrame.from_dict(result, orient="index")
            
            # Make biomass a column instead of index
            df.index.name = "biomass_flux"
            df.reset_index(inplace=True)
            
            # Write to sheet (Excel sheet names max length = 31)
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    print(f"Saved to {filename}")

export_yields_to_excel(yield_dict, "figureC_yields.xlsx")

result_dict={
    "52_unsup": [unsup_52, "#010101"], 
    "52_indole": [I_data, '#7d277d'],
    "52_quln": [Q_data, '#3a53a4'],
    "52_Ribose_25": [R25_data,'#818181'],
    "52_Ribose_100": [R100_data,'#faa51a'],
    "24_unsup": [unsup_24, "#a0a0a0"], 
    "24_indole": [I_data_24, '#c541c9'],
    "24_quln": [Q_data_24, '#4d75d6'],
    "24_Ribose_25": [R25_data_24,'#dddddd'],
    "24_Ribose_100": [R100_data_24,'#f9e65d']}

def big_PE(PE_result_dict, save_dir, title, xlab, ylab, filename):

    os.makedirs(save_dir, exist_ok=True)

    plt.rcParams.update({
        "font.family": "Helvetica",
        "axes.linewidth": 1.2,
        "axes.edgecolor": "gray",
        "grid.color": "lightgray",
        "grid.linestyle": "--",
        "grid.linewidth": 0.8,
        "legend.frameon": False,
    })

    plt.figure(figsize=(8, 6))

    for label, (results, color) in PE_result_dict.items():
        bm_vals = np.array(sorted(results.keys()))
        maxs = np.array([results[b]['max_target'] for b in bm_vals])

        # optional: append final drop to zero
        bm_vals_plot = np.append(bm_vals, bm_vals[-1])
        maxs_plot = np.append(maxs, 0)

        plt.plot(
            bm_vals_plot,
            maxs_plot,
            marker='o',
            linestyle='-',
            color=color,
            label=label,
            alpha=0.9
        )

    plt.title(title)
    plt.xlabel(xlab, fontsize=12)
    plt.ylabel(ylab, fontsize=12)
    plt.grid(alpha=0.3)
    plt.legend(fontsize=10)
    plt.tight_layout()

    outpath = os.path.join(save_dir, filename)
    if not outpath.endswith(".svg"):
        outpath += ".svg"

    plt.savefig(outpath, format='svg')
    plt.show()
    plt.close()
    
big_PE(
    result_dict,
    "directory",
    "Checking Growth Coupling",
    "Biomass Flux",
    "I3A Secretion Flux",
    "I3A supplementation"
)