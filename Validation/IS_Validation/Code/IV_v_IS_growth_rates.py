################ Comparing in silico and in vitro LGG growth rates in various DMs ################
# This script hard codes the growth rates from both in vitro data (IV) and the IS data from the full media control in the the DO simulations

import matplotlib.pyplot as plt
import numpy as np

models=["iGR413", "AGORA", "MERLIN"]
IS_data={"DM52": [0.428, 0.0, 0.0],
      "DM24": [0.313, 0.0, 0.0],
      "DM16": [0.134, 0.0, 0.0],
      "DM13": [0.134, 0.0, 0.0]}

IV_data={"DM52": [0.224,0.214,0.226,0.219,0.227,0.22,0.23,0.234,0.251,0.233,0.231,0.227,0.215,0.21,0.221,0.218,0.214,0.178,0.236,0.246,0.248,0.236,0.228,0.256],
      "DM24": [0.276,0.276,0.296,0.278,0.273,0.268,0.275,0.287,0.409,0.289,0.294,0.283],
      "DM16": [0.147,0.152,0.146,0.137,0.15,0.171,0.21,0.22,0.209,0.25,0.27,0.24],
      "DM13": [0.239,0.243,0.235,0.248,0.236,0.211]}

colors = {"iGR413": "#ec2727", "AGORA": "#ee9c4b", "MERLIN": "#3c54a4"}
media = list(IS_data.keys())

# Plot setup
fig, ax = plt.subplots(figsize=(8,6))
x_positions = np.arange(len(media))  

rng = np.random.default_rng(42) 

for i, medium in enumerate(media):
    # Plot IV replicates with jitter
    y_vals = IV_data[medium]
    jitter = rng.uniform(-0.2, 0.2, size=len(y_vals))
    ax.scatter(x_positions[i] + jitter, y_vals,
               facecolors="none", edgecolors="black", alpha=0.7, label="IV replicates" if i==0 else "")
    
    # Plot IV mean line
    mean_val = np.mean(y_vals)
    ax.hlines(mean_val, x_positions[i]-0.3, x_positions[i]+0.3,
              colors="black", linestyles="dashed", label="IV mean" if i==0 else "")
    
    # Plot IS stars (1 per model, slightly spaced)
    for j, model in enumerate(models):
        ax.scatter(x_positions[i] + (j-1)*0.15, IS_data[medium][j],
                   marker="*", s=200, color=colors[model], edgecolor="black", linewidth=1.0,
                   label=model if i==0 else "")

# Formatting
ax.set_xticks(x_positions)
ax.set_xticklabels(media)
ax.set_ylabel("Growth Rate")
ax.legend(
    loc="upper left",
    bbox_to_anchor=(0.75, 0.95),
    frameon=True,
    facecolor="white",
    edgecolor="black"
)


plt.tight_layout()
plt.savefig('ISvIV_growthrates.svg', format='svg')
plt.show()
