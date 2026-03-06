from set_dms import set_dm 
from pe_utils import pe_Data, plot_production_envelope_single
from cobra.io import read_sbml_model 
import numpy as np
import matplotlib.pyplot as plt

############################ Generating single PEs for different metabolites in iGR413 ############################
#   Uses two functions from pe_utils to generate seperate production envelopes for lactate and i3a against biomass
#   A dictionary of media conditions and corresponding colors is looped through to:
        # Set proper media conditions using set_DM
        # Generate lactate and I3A production envelope-ready data dicts using pe_Data
        # Single curve production envelopes for each metabolite and media condition are generated using plot_production_envelope_single 
model= read_sbml_model('iGR413.xml')
media={'52': '#6e9869', '24': '#5c67a8'}
for media_name, color in media.items():
    model_run= model.copy()
    model_set, media_dict= set_dm(model_run, media_name) 
    results_lac=pe_Data(model_set, target='EX_lac_L_e', num_points=10)
    results_i3a=pe_Data(model_set, target='EX_I3A_e', num_points=10)
    title_lac= str( media_name + " Lactate Flux PE")
    filename_lac= str( media_name + "_lac_PE")
    title_i3a= str( media_name + " I3A Flux PE")
    filename_i3a= str( media_name + "_i3a_PE")
    plot_production_envelope_single(results_lac, '/Users/grichmond/Desktop/Lactate_PEs',title_lac, "biomass flux",'lactate secretion flux',filename_lac, color)
    plot_production_envelope_single(results_i3a, '/Users/grichmond/Desktop/I3A_PEs',title_i3a, "biomass flux",'I3A secretion flux',filename_i3a, color)

############################ Identifying Phases of Changing Slope ############################
#   Uses two functions from pe_utils to generate identify regions where the slope changes for lactate PEs
#   A dictionary of media conditions and corresponding colors is looped through to:
        # Set proper media conditions using set_DM
        # Generate lactate production envelope-ready data dicts using pe_Data
        # identify periods of slope changes using detect_slope_changes_pe
        # Using the results, shading was performed in Adobe Illustrator to visually identify detected phases 

def detect_slope_changes_pe(pe_dict, which='max', threshold=0.1, return_absolute_diff=False):
    """
    Detect slope changes in production envelope data from the result of .

    Parameters
    ----------
    pe_dict : dict
        Mapping objective_value -> [target_min, target_max].
    which : str
        'max' or 'min' — which curve to analyze.
    threshold : float
        Minimum absolute slope difference to call a change.
    return_absolute_diff : bool
        If True, 'slope_diff' is abs(slope_AB - slope_BC). Otherwise signed.

    Returns
    -------
    change_points : list of dict
        Each dict includes:
            'A','B','C' : objective flux values
            'flux_A','flux_B','flux_C' : corresponding target flux values
            'slope_AB','slope_BC','slope_diff','index_B'
    """
    if not isinstance(pe_dict, dict) or len(pe_dict) < 3:
        return []
    xs = np.array(sorted(pe_dict.keys()), dtype=float)
    idx = 1 if which == 'max' else 0
    ys = np.array([pe_dict[x][idx] for x in xs], dtype=float)
    dx = np.diff(xs)
    dy = np.diff(ys)
    slopes = np.divide(dy, dx, out=np.zeros_like(dy), where=dx != 0)
    change_points = []
    for i in range(1, len(slopes)):
        slope_AB = slopes[i - 1]
        slope_BC = slopes[i]
        if np.isnan(slope_AB) or np.isnan(slope_BC):
            continue
        diff = abs(slope_AB - slope_BC) if return_absolute_diff else slope_AB - slope_BC
        if abs(slope_AB - slope_BC) > threshold:
            record = {
                'A': xs[i - 1], 'B': xs[i], 'C': xs[i + 1],
                'flux_A': ys[i - 1], 'flux_B': ys[i], 'flux_C': ys[i + 1],
                'slope_AB': slope_AB, 'slope_BC': slope_BC,
                'slope_diff': diff,
                'index_B': i
            }
            change_points.append(record)
    return change_points

model= read_sbml_model('iGR413.xml')
media={'52': '#6e9869', '24': '#5c67a8'}
changed={}
for media_name, color in media.items():
    model_run= model.copy()
    model_set, media_dict= set_dm(model_run, media_name) 
    results_lac=pe_Data(model_set, target='EX_lac_L_e', num_points=10)
    x=detect_slope_changes_pe(results_lac)
    changed[media]= x
print(changed)

############################ Characterizing Changes in Rxns Across Slope Change Phases ############################
#   Using the slope phases identified, this step creates figures that allow for the visual assessment of the reactions changing 
#   Uses additional scripts to perform flux value extaction, calculate differences, filter by meaningful differences, and sort reactions into categories 

def get_flux_results(model, biomass_values, objective='curated_biomass'):
    results = {}
    for bm_val in biomass_values:
        mcopy = model.copy()
        mcopy.reactions.get_by_id(objective).bounds = (0, bm_val)
        sol = fba(mcopy)
        results[bm_val] = sol.fluxes
    return results

def compare_fluxes(results_dict, val1, val2, abs_tol=1e-6):
    """
    Compare flux distributions between two conditions, normalized by biomass.

    Parameters:
        results_dict (dict): Dictionary of {value: flux_series}
        val1 (float): First biomass value (also key in results_dict)
        val2 (float): Second biomass value (also key in results_dict)
        abs_tol (float): Tolerance for considering flux changes as nonzero

    Returns:
        DataFrame with normalized fluxes and their differences
    """
    # Normalize all fluxes by their biomass value
    flux1 = results_dict[val1] / val1
    flux2 = results_dict[val2] / val2

    # Calculate difference
    diff = flux2 - flux1
    df = pd.DataFrame({
        f'flux_{val1}_norm': flux1,
        f'flux_{val2}_norm': flux2,
        'flux_difference': diff
    })

    # Keep only reactions above tolerance
    changed_df = df[df['flux_difference'].abs() > abs_tol]
    changed_df = changed_df.sort_values(by='flux_difference', key=abs)

    return changed_df

def filter_flux_changes(df, threshold=1e-3):
    """
    Filter reactions with flux difference above a threshold.
    
    Parameters:
        df (DataFrame): Output from compare_fluxes()
        threshold (float): Minimum absolute difference to include
        
    Returns:
        Filtered DataFrame
    """
    return df[df['flux_difference'].abs() > threshold].copy()

def exclude_reactions_from_df(df, exclude_exact=None, exclude_prefix=None):
    """
    Exclude reactions from the DataFrame by exact name or by prefix.

    Parameters:
        df (DataFrame): Flux comparison DataFrame with reaction IDs as index.
        exclude_exact (list): List of exact reaction IDs to exclude.
        exclude_prefix (list): List of string prefixes to exclude.

    Returns:
        Filtered DataFrame
    """
    exclude_exact = exclude_exact or []
    exclude_prefix = exclude_prefix or []

    # Build exclusion mask
    mask = df.index.isin(exclude_exact)
    for prefix in exclude_prefix:
        mask |= df.index.str.startswith(prefix)

    return df[~mask].copy()

def plot_flux_changes_4_panels(df, title_prefix="Filtered Flux Changes (>0.05) Normalized to Biomass: ",
                                extracellular_containing=None, baseline_value=None):
    extracellular_containing = extracellular_containing or []

    # Classify reactions
    is_extracellular_related = df.index.isin(extracellular_containing)
    df_extracellular = df[is_extracellular_related]
    df_metabolic = df[~is_extracellular_related]

    # Split into up/down
    up_ex = df_extracellular[df_extracellular['flux_difference'] > 0].sort_values(by='flux_difference', ascending=False)
    down_ex = df_extracellular[df_extracellular['flux_difference'] < 0].sort_values(by='flux_difference')
    up_met = df_metabolic[df_metabolic['flux_difference'] > 0].sort_values(by='flux_difference', ascending=False)
    down_met = df_metabolic[df_metabolic['flux_difference'] < 0].sort_values(by='flux_difference')

    # Subplot layout
    fig, axs = plt.subplots(4, 1, figsize=(16, 14))
    plt.suptitle(title_prefix, fontsize=16)

    def plot_on_ax(ax, data, title, color):
        if data.empty:
            ax.set_visible(False)
            return

        bars = ax.bar(data.index, data['flux_difference'], color=color)
        ax.set_title(title, fontsize=12)
        label = f"Flux Change relative to {baseline_value}" if baseline_value is not None else "Flux Change"
        ax.set_ylabel(label)
        ax.axhline(0, color='gray', linewidth=0.8)
        ax.tick_params(axis='x', rotation=90, labelsize=8)

        # Annotate small bars with rotated vertical text
        for bar, val in zip(bars, data['flux_difference']):
            if abs(val) < 0.5:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    val + 0.01 if val > 0 else val - 0.01,
                    f"{val:.3f}",
                    ha='center',
                    va='bottom' if val > 0 else 'top',
                    fontsize=7,
                    rotation=90
                )

    plot_on_ax(axs[0], up_ex, "↑ Extracellular-Related Reactions", 'skyblue')
    plot_on_ax(axs[1], down_ex, "↓ Extracellular-Related Reactions", 'salmon')
    plot_on_ax(axs[2], up_met, "↑ Internal Metabolic Reactions", 'skyblue')
    plot_on_ax(axs[3], down_met, "↓ Internal Metabolic Reactions", 'salmon')

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for title
    plt.show()

extracellular_containing = []
for rxn in model.reactions:
    if any("_e" in met.id for met in rxn.metabolites):
        extracellular_containing.append(rxn.id)

media={'52': '#6e9869', '24': '#5c67a8'}
for media_name, color in media.items():
    model_set, media_type = set_dm(model.copy(), media_name)
    bm_values = sorted({cp['A'] for cp in changed[media_name]} | {cp['B'] for cp in changed[media_name]} | {cp['C'] for cp in changed[media_name]})
    bm_results = get_flux_results(model_set, bm_values)
    for i in range(len(bm_values) - 1):
    val1, val2 = bm_values[i], bm_values[i + 1]
    desc = f"{val1:.5f} vs {val2:.5f}"
    df = compare_fluxes(bm_results, val1, val2)
    df_filtered = filter_flux_changes(df, threshold=0.05)
    df_clean = exclude_reactions_from_df(df_filtered, exclude_exact=['curated_biomass'], exclude_prefix=['e_', 'DM_'])
    plot_flux_changes_4_panels(df_clean, title_prefix=f"Flux Changes {desc} (Medium {media_name}, Lactate)",extracellular_containing=extracellular_containing,baseline_value=val1)
