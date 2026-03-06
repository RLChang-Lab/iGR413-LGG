################################## Figure Generation of Memote Barchart (Fig 1C) ##################################
#   Using data from the memote HTML in the supplemental materials, a bar chart is made from MEMOTE Diff function 
#   Data from HTML is copied into this script via a dict called "data" 

import matplotlib.pyplot as plt
import numpy as np

models=["iGR413", "MERLIN", "AGORA"]

data={'SBO Annotation' : [33,0,0],
'Gene Annotation' : [33,40,0],
'Metabolite Annotation' : [54,30,24],
'Reaction Annotation' : [60,50,22],
'Consistency' : [82,86,43],
'Total Score' : [61,41,20]}

scores = np.array(list(data.values()))  # shape: (num_categories, 3)
categories = list(data.keys())
x = np.arange(len(categories))
width = 0.25

hatches = ['', '\\', 'x']

fig, ax = plt.subplots(figsize=(10, 6))

for i, (model, hatch) in enumerate(zip(models, hatches)):
    bars = ax.bar(x + i*width, scores[:, i], width, color='lightgray', edgecolor='black', hatch=hatch, label=model)
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + 1, f'{int(height)}', ha='center', va='bottom', fontsize=12)

ax.set_xticks(x + width) 
ax.set_xticklabels(categories, rotation=45, ha='right', fontsize=15)
ax.tick_params(axis='y', labelsize=12)
ax.set_ylabel('MEMOTE Score (%)', fontsize=18)
ax.set_ylim(0, 100)

# Legend inside plot, aligned to left
ax.legend(fontsize=20, loc='upper left', bbox_to_anchor=(0.02, 0.98))

plt.tight_layout()
plt.savefig('memote_scores.png', dpi=1200, bbox_inches='tight') 
plt.show()
