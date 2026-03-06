import cobra
def fba(model, objective='curated_biomass'):
     """
    Runs FBA simulation on a cobra model object for maximization of specified objective function 

    Args:
        model (Cobra model): CB model to run simulation on 
        objective (str): desired objective function for simulation. Default is iGR413 LGG model biomass function 'curated_biomass', but could be any model rxn  
    Returns: 
        result (cobra.core.solution.Solution): Cobrapy solution object 
    Notes:
        - Model is constrained previously (example, set_dm must be run first to simulate DM condition)
        - KeyError if objective reaction input does not exist in model 
        - Solution object contains: 
            fluxes (pandas.core.series.Series) 
            get_primal_by_id function (method)
            objective value (float)
            shadow prices (pandas.core.series.Series)
            reduced costs (pandas.core.series.Series)
            to_frame function (method)
            status (str) 
        - Cobrapy Version: 0.29.1, Python Version: 3.9.12
    """
        #Set objective function 
        model.objective=model.reactions.get_by_id(objective)
        #designation optimization directionality. 
        model.objective_direction='max'
        #perform optimization 
        result=model.optimize()
        return result
