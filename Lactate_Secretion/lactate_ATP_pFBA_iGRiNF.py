
import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from set_dms import set_dm

########## Functions Required  ##########
#   phase_minMax_pFBA works by running FVA on a range-locked biomass objective and generating a reduced model where only reactions carrying flux and their metabolites are retained. 
        # Then, if lactate dehydrogenase reaction is still retained in the model, this reaction is set to 0 to prevent any lactate related fermentation 
        # fluxes for mets and rxns respectively are calculated from pFBA solution vector of the reduced model \
# make_cytoscape_edges converts the fluxes from phase_minMax_pFBA into a csv table readable by cytoscape, while excluding common secondary metabolites from the file. 
def phase_minMax_pFBA(model, objective_rxn):
    results={}
    x=flux_variability_analysis(model)
    reduced_model=model.copy()
    to_remove=[]
    for rxn, results in x.iterrows():
        if results[0]<=1e-6 and results[1]<=1e-6:  
            to_remove.append(rxn)
    if "LDH_L" not in to_remove: 
        print("lactate had non zero FVA flux")
        model.reactions.LDH_L.bounds=0,0
    else:
        print("lactate removed in FVA step") 
    to_remove = [reduced_model.reactions.get_by_id(rxn_id) for rxn_id in to_remove]
    reduced_model.remove_reactions(to_remove, remove_orphans=True) 
    y=cobra.flux_analysis.parsimonious.pfba(model, objective= model.reactions.get_by_id(objective_rxn))
    return y.fluxes
import cobra
def make_cytoscape_edges(model, flux_dict=None, exclude=None):
    """
    Create Cytoscape-compatible metabolite→metabolite edge list from a COBRA model,
    filtering out common cofactors like ATP, NADH, etc.
    
    Args:
        model (cobra.Model): COBRA model
        flux_dict (dict): Optional {rxn_id: flux_value}
        exclude (list): metabolite IDs or names to ignore (case-insensitive substrings)
    
    Returns:
        pd.DataFrame: edge list with source, target, rxn_id, flux, reaction_name
    """
    # Default cofactors to exclude
    if exclude is None:
        exclude = [
            'atp', 'adp', 'amp', 'nad', 'nadh', 'nadp', 'nadph',
            'h', 'h2o', 'pi', 'ppi', 'coa', 'co2', 'o2', 'nh4', 'trdox', 'trdrd'
        ]

    edges = []

    for rxn in model.reactions:
        if flux_dict and rxn.id not in flux_dict:
            continue

        flux_value = flux_dict.get(rxn.id, 0) if flux_dict else None

        reactants = [met for met, coeff in rxn.metabolites.items() if coeff < 0]
        products  = [met for met, coeff in rxn.metabolites.items() if coeff > 0]

        # Filter out cofactors by substring match (case-insensitive)
        def is_excluded(met):
            return any(ex in met.id.lower() or ex in met.name.lower() for ex in exclude)

        reactants = [m for m in reactants if not is_excluded(m)]
        products  = [m for m in products if not is_excluded(m)]

        # skip if reaction now has no valid connections
        if not reactants or not products:
            continue

        for r in reactants:
            for p in products:
                edges.append({
                    "source": r.id,
                    "target": p.id,
                    "reaction": rxn.id,
                    "flux": flux_value,
                    "reaction_name": rxn.name
                })

    edges_df = pd.DataFrame(edges)
    return edges_df

########## Asessing iGR632  ##########
#   This code does the following steps
        # Replace the biomass function with demands for each pseudometabolite to better assess each pathway, rather than how they might work together to support biomass 
        # Sets DM conditions and creates a seperate model foe each condition
        # runs the specific pFBA simulation to get the minimum path to ATP maintaince in each condition
        # Filters out all of the 0 flux reactions and generates a cytoscape input file for visualization 

ATPmodel = cobra.io.read_sbml_model('iGR_ATPmain.xml')
ATPmodel.add_boundary(ATPmodel.metabolites.e_Protein_c, type='demand')
ATPmodel.add_boundary(ATPmodel.metabolites.e_Cofactor_c, type='demand')
ATPmodel.add_boundary(ATPmodel.metabolites.e_RNA_c, type='demand')
ATPmodel.add_boundary(ATPmodel.metabolites.e_DNA_c, type='demand')
ATPmodel.add_boundary(ATPmodel.metabolites.M8801_c, type='demand')
ATPmodel.add_boundary(ATPmodel.metabolites.M8807_c, type='demand')
ATPmodel.add_boundary(ATPmodel.metabolites.M8811_c, type='demand')
ATPmodel.add_boundary(ATPmodel.metabolites.e_Lipid_c, type='demand')
ATPmodel.remove_reactions(['BM_noATP'])
lac_tmodelc24, media24 = set_dm(ATPmodel, '24')
lac_tmodelc52, media52 = set_dm(ATPmodel, '52')

phase_24_lac_results={}
fluxes_24 = phase_minMax_pFBA(lac_tmodelc24, 'ATPM')

phase_52_lac_results={}
fluxes_52=phase_minMax_pFBA(lac_tmodelc52, 'ATPM')
    
nonzero_rxns_52={}
for rxn, flux in fluxes_52.items():
    if flux != 0:
        nonzero_rxns_52[rxn]=flux
            
nonzero_rxns_24={}
for rxn, flux in fluxes_24.items():
    if flux != 0:
        nonzero_rxns_24[rxn]=flux

edges_52 = make_cytoscape_edges(lac_tmodelc52, nonzero_rxns_52)
edges_24 = make_cytoscape_edges(lac_tmodelc24, nonzero_rxns_24)

edges_52["condition"] = "52"
edges_24["condition"] = "24"

edges_52.to_csv("cleaned_edges_52.csv", index=False)
edges_24.to_csv("cleaned_edges_24.csv", index=False)

rxns_52 = set(nonzero_rxns_52.keys())
rxns_24 = set(nonzero_rxns_24.keys())

shared_rxns = rxns_52 & rxns_24
unique_52 = rxns_52 - rxns_24
unique_24 = rxns_24 - rxns_52

edges_combined = pd.concat([edges_52, edges_24], ignore_index=True)

#this also categorizes if a reaction from the pFBA solutions are unique to media type. In practice, all reactions where shared
def categorize_reaction(rxn_id):
    if rxn_id in shared_rxns:
        return "shared"
    elif rxn_id in unique_52:
        return "unique_52"
    elif rxn_id in unique_24:
        return "unique_24"
    else:
        return "unknown"

edges_combined["rxn_category"] = edges_combined["reaction"].apply(categorize_reaction)
edges_combined.to_csv("cleaned_cytoscape_merged.csv", index=False)



########## Asessing iNF517  ##########
#   This code does the following steps to the Lactococcus lactis model in BiGG 
        # Generates a version of the model where biomass precusors are decoupled from each other 
        # Sets DM conditions and creates a seperate model foe each condition
        # runs the specific pFBA simulation to get the minimum path to ATP maintaince in each condition
        # Filters out all of the 0 flux reactions and generates a cytoscape input file for visualization import cobra
import os 
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
model = cobra.io.read_sbml_model('iNF517.xml')
for i in model.reactions.BIOMASS_LLA_noATPnoH.metabolites:
    model.add_boundary(i, type='demand')
model.remove_reactions(['BIOMASS_LLA','LDH_L', 'L_LACt2r', 'MALLAC'])
for i in model.exchanges:
    i.lower_bound=-10

fluxes = phase_minMax_pFBA(model, 'ATPM') 
nonzero_rxns={}
for rxn, flux in fluxes.items():
    if flux != 0:
        nonzero_rxns[rxn]=flux          
edges = make_cytoscape_edges(model, nonzero_rxns)
edges.to_csv("iNF_edges.csv", index=False)
