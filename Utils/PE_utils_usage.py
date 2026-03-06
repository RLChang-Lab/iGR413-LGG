from cobra
from run_fba import fba
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt

def pe_Data(model, objective='curated_biomass', target, num_points=10):
    """
    Generates the data needed to plot a production envlope between the objective function and a target function across the full range of objective flux values

    Args:
        model (Cobra model): CB model to run simulation on 
        objective (str): desired objective function for simulation. Default is iGR413 LGG model biomass function 'curated_biomass', but could be any model rxn. X-axis on the production envelope 
        target (str): desired target reaction to interrogate. Y-axis on the production envelope
        num_points (int): number of objective value points to be assessed in the FVA range 
    Returns: 
        results (dict): Dictionary of objective value, target value min and max value pairs for plotting
    Notes:
        - Model is constrained previously (example, set_dm must be run first to simulate DM condition)
        - KeyError if objective reaction input does not exist in model 
        - Use  plot_production_envelope to visualize 
        - Cobrapy Version: 0.29.1, Python Version: 3.9.12
    """
    base_fba = fba(model, objective=objective) #Identify max objective flux
    max_bm = float(base_fba.objective_value) 
    bm_values = np.linspace(0, max_bm, num_points) #determine set of objective values to be assessed between 0 and max flux
    results={}
    for b in bm_values:
        mcopy = model.copy()
        mcopy.reactions.get_by_id(objective).bounds = b,b   #constrain objective to specific value
        sol= flux_variability_analysis(mcopy, processes=1)
        target_min=sol.loc[target, 'minimum'] #extract minimum and maximum flux through target at that objective value
        target_max=sol.loc[target, 'maximum']
        results[b]= target_min, target_max 

def plot_production_envelope_dual(results1, results2, label1, label2, save_dir,title, xlab,ylab,filename):
    """
    Plot a single production envelope comparing two treatments.
    Args:
        results1/2 (dict): results dicts from pe_Data function for the two conditions 
        label1/2 (str): label of each condition being compared. 
        save_dir (str): directory to save figure file to 
        title, xlab, ylab (str): naming features on the plot
        filename (str): base name of the desired output file of figure. No extension

    Returns: 
        Figure as svg named as filename at save_dir
    Notes:
        - Model is constrained previously (example, set_dm must be run first to simulate DM condition)
        - Must run pe_Data prior to running this
        - This PE script was used to generate iNF561's lactate production envelop 
        - Only used on comparing two conditions 
        - Cobrapy Version: 0.29.1, Python Version: 3.9.12
    """
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

    # Treatment 1
    bm_vals = np.array(sorted(results1.keys()))
    maxs = np.array([results1[b][1] for b in bm_vals])
    bm_vals = np.append(bm_vals, bm_vals[-1])
    maxs = np.append(maxs, 0)
    plt.plot(bm_vals, maxs, 'o-', color="#1f77b4", label=label1, alpha=0.9)

    # Treatment 2
    bm_vals_c = np.array(sorted(results2.keys()))
    maxs_c = np.array([results2[b][1] for b in bm_vals_c])
    bm_vals_c = np.append(bm_vals_c, bm_vals_c[-1])
    maxs_c = np.append(maxs_c, 0)
    plt.plot(bm_vals_c, maxs_c, 's--', color="#ff7f0e", label=label2, alpha=0.9)

    plt.title(title) 
    plt.xlabel(xlab, fontsize=12)
    plt.ylabel(ylab, fontsize=12)
    plt.grid(alpha=0.3)
    plt.legend(fontsize=10)
    plt.tight_layout()
    filename= str(filename + ".svg")
    plt.savefig(filename, format='svg') 
    plt.show()
    plt.close()

def plot_production_envelope_single(results1, label1, save_dir,title, xlab,ylab,filename, color_pick):
    """
    """
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

    # Treatment 1
    bm_vals = np.array(sorted(results1.keys()))
    maxs = np.array([results1[b][1] for b in bm_vals])
    bm_vals = np.append(bm_vals, bm_vals[-1])
    maxs = np.append(maxs, 0)
    plt.plot(bm_vals, maxs, 'o-', color=color_pick, label=label1, alpha=0.9)
    plt.title(title) 
    plt.xlabel(xlab, fontsize=12)
    plt.ylabel(ylab, fontsize=12)
    plt.grid(alpha=0.3)
    plt.legend(fontsize=10)
    plt.tight_layout()
    filename= str(filename + ".svg")
    plt.savefig(filename, format='svg') 
    plt.show()
    plt.close()
