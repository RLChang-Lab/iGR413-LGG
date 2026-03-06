################## Rxn Confidence Figure ##################
#   Generates a waffle chart of counts of rxns with condifence score
#   Data used to generate was identified from counts of each category within the model file 

import matplotlib.pyplot as plt
from pywaffle import Waffle

#Read in data 
data = {'Reaction w/ GPR' : 5,
'Reaction w/o GPR' : 5,
'Pathway w/ GPR' : 129,
'Pathway w/o GPR' : 99,
'GPR only' : 723,
'Model' : 216}

# Basic waffle
plt.figure(
    FigureClass=Waffle,
    rows=10,
    columns=30,
    values=data,
    figsize=(10, 5), 
    colors=["#000000", "#666766", "#999999", "#cccccb", "#eeeeee", "#ffffff"],
    icons='square',
    legend={'loc': 'upper left', 'bbox_to_anchor': (1.05, 1), 'fontsize': 8},
    labels=data,
    font_size=15)

#plt.show()
plt.savefig('rxn_confidence.svg', format='svg') 


