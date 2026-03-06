############## Summarizing model performance from Validation Script ################ 
# Hard codes the model performance metrics from the validation script into a table using a heatmap structure

import matplotlib.pyplot as plt
import numpy as np

metrics = ['Accuracy', 'Recall', 'Precision', 'Specificity', 'FPR']
media = ['DM52', 'DM24', 'DM16', 'DM13', 'Sun et al.']
data_dict = {
    'DM52': [75, 14, 100, 100, 0],
    'DM24': [67, 38, 60, 85, 15],
    'DM16': [53, 40, 33, 60, 40],
    'DM13': [50, 57, 57, 40, 60],
    'Sun et al.': [57, 15, 75, 96, 5]
}

data_pct = np.array([data_dict[m] for m in media]).T  # transpose so rows=metrics
media_colors = ['#6f9969', '#5c66a8', '#808fe1', '#454a74']

fig, ax = plt.subplots(figsize=(8, 5))

im = ax.imshow(data_pct, cmap='Greys', vmin=0, vmax=100)

ax.set_xticks(np.arange(len(media)))
ax.set_yticks(np.arange(len(metrics)))
ax.set_xticklabels(media)
ax.set_yticklabels(metrics)

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

for tick_label, color in zip(ax.get_xticklabels(), media_colors):
    tick_label.set_color(color)
    tick_label.set_rotation(45)
    tick_label.set_ha('left')

for i in range(len(metrics)):
    for j in range(len(media)):
        text_color = "white" if data_pct[i, j] > 50 else "black"
        ax.text(j, i, f"{data_pct[i, j]}%", ha="center", va="center", color=text_color, fontsize=12)

plt.tight_layout()
plt.savefig('model_metrics.svg', format='svg') 
plt.show()
