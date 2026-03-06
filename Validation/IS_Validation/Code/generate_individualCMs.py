##### Generating final confusion matrix figure from data in Validation script #####

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec

data = {
    'DM52' : [33, 0, 12, 2],
    'DM24' : [11, 2, 5, 3],
    'DM16' : [6, 4, 3, 2],
    'DM13' : [2, 3, 3, 4],
    'Sun et al.' : [21, 1, 17, 3]
}

media_colors = {
    'DM52': '#6f9969',
    'DM24': '#5c66a8',
    'DM16': '#808fe1',
    'DM13': '#454a74'
}

def plot_confusion_matrix(ax, vals, title, text_color='black'):
    for i in range(2):
        for j in range(2):
            rect = Rectangle((j-0.5, i-0.5), 1, 1, facecolor='white', edgecolor='black', linewidth=2)
            ax.add_patch(rect)
            ax.text(j, i, str(vals[i*2 + j]), ha='center', va='center', fontsize=30, color=text_color)

    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.5, 1.5)
    
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['0', '-'], fontsize=20, rotation=0, ha='center', va='center')
    ax.set_xlabel("IS", fontsize=20)
    ax.xaxis.set_label_coords(0.5, -0.2)
    ax.tick_params(axis='x', which='major', pad=15)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(['0', '-'], fontsize=20, rotation=0, ha='right', va='center')
    ax.set_ylabel("IV", fontsize=20, rotation=0)
    ax.yaxis.set_label_coords(-0.25, 0.45)
    ax.tick_params(axis='y', which='major', pad=10)

    ax.set_title(title, fontsize=30, color=text_color)
    ax.set_aspect('equal')
    ax.invert_yaxis()

# ---- Figure layout ----
fig = plt.figure(figsize=(20, 6))  # wider to fit 5 matrices
gs = GridSpec(1, 5, figure=fig, wspace=0.4)  # 1 row, 5 columns

# Plot all matrices in a single row
all_labels = ['DM52', 'DM24', 'DM16', 'DM13', 'Sun et al.']
for col, label in enumerate(all_labels):
    ax = fig.add_subplot(gs[0, col])
    color = media_colors.get(label, 'black')
    plot_confusion_matrix(ax, data[label], label, text_color=color)

plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
plt.savefig('individual_CMs.svg', format='svg')
plt.show()

