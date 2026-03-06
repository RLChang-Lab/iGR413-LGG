from cobra.io import read_sbml_model
from PE_utils import pe_Data, plot_production_envelope_dual

##################################### Generating a Two Condition PE for Lactococcus Lactis Lactate #####################################
#   Since the model is preconstrained per the citation, we compared both all exchanges open to -10 and the publication reported constraints on the same PE
#   The script functions by generating 2 pe_Data result dicts from PE_utils and then plotting them together using plot_production_envelope_dual
model = read_sbml_model('iNF517.xml')
paper_constrained=pe_Data(model, objective='BIOMASS_LLA', target='EX_lac__L_e', num_points=10)
LB_constrained=model.copy()
for i in LB_constrained.exchanges:
    i.lower_bound=-10
unconstrained= pe_Data(LB_constrained, objective='BIOMASS_LLA', target='EX_lac__L_e', num_points=10)
plot_production_envelope_dual(paper_constrained, LB_constrained, 'Paper Constrained', 'All LBs = -10', '/Users/grichmond/Desktop',"Lactococcus Lactis Lactate PE", "Biomass Flux","Lactate Secretion Flux","iNF_lactate_PE")
