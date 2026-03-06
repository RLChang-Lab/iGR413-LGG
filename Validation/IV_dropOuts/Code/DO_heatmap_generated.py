############## In Vitro Media Dropout Heatmap Figure ##############
#   Takes the average fold change for each component and each media formulation, sorts, and plots as a heatmap. 

import matplotlib.pyplot as plt
import numpy as np

components = ['Cys',
'Magnesium Sulfate',
'Iso',
'Potassium Phosphate Dibasic',
'Glu',
'Asp',
'Asn ',
'Arg',
'His',
'Val',
'Cobalt Sulfate',
'Pro',
'Cytidine',
'MOPS',
'Ammonium Chloride',
'Citrate',
'Uracil',
'Myo-inositol',
'Xanthine',
'Aminobutyric Acid',
'Methionine',
'Ser',
'Ascorbic Acid',
'Ala',
'Lys',
'Trp',
'Gly',
'Leu',
'Tyro',
'Thr',
'Phe',
'Adenine',
'Guanine',
'Cytosine',
'Thiamin',
'Nicotinic Acid',
'Pantothenate',
'Pyroxidal',
'Aminobenzoic acid',
'Folate',
'Zinc Sulfate',
'Copper Sulfate',
'Borate',
'Biotin',
'B12',
'Riboflavin',
'Mn SO4',
'Thymine',
'Sodium Chloride',
'Calcium Chloride',
'Iron Sulfate']
media = ['DM52', 'DM24', 'DM16', 'DM13']
data_dict = {
   'Cys' : [0.580689655172414,0.472705458908218,0.206433163596381,0.100160884929866],
'Magnesium Sulfate' : [0.43560411311054,1.14480130647795,0.137763012181617,0.103330100226318],
'Iso' : [0.332165605095541,0.548127853881278,0.242530375374834,0.186323381122545],
'Potassium Phosphate Dibasic' : [0.540231362467866,0.599346761023408,0.356677740863787,0.255221467830585],
'Glu' : [0.875,0.73972602739726,0.808434584582779,0.570394115785339],
'Asp' : [1.10344827586207,1.00547945205479,0.904074325888859,0.67739329406383],
'Asn ' : [0.626349614395887,0.857534246575343,1.15968547305668,0.78519697744216],
'Arg' : [0.684713375796178,0.557808219178082,0.705503921293725,0.843925849730354],
'His' : [0.855172413793103,1.08310502283105,1.1129377968129,0.881737589422753],
'Val' : [0.886942675159236,1.17936412717457,1.0874861572536,0.93566117038474],
'Cobalt Sulfate' : [0.590745501285347,0.486502699460108,0.880629008139686,0.990345777900921],
'Pro' : [0.630573248407643,0.711415525114155,0.887776970868092,1.07803683378329],
'Cytidine' : [0.746863722289053,0.874625074985003,0.8859357696567],
'MOPS' : [0.633161953727506,1.29994556341862,0.960132890365449],
'Ammonium Chloride' : [0.682390745501285,0.0786608600979859,1.35658914728682],
'Citrate' : [0.83624678663239,0.634186173108329],
'Uracil' : [0.757690324798075,0.923815236952609],
'Myo-inositol' : [0.617429305912596,0.929214157168566],
'Xanthine' : [0.838803918199003,0.93881223755249],
'Aminobutyric Acid' : [1.32706766917293,0.95380923815237],
'Methionine' : [0.863057324840764,0.959208158368326],
'Ser' : [0.787420382165605,0.980821917808219],
'Ascorbic Acid' : [0.66413881748072,1.01735159817352],
'Ala' : [0.957006369426751],
'Lys' : [0.969745222929936],
'Trp' : [1.2351724137931],
'Gly' : [1.23448275862069],
'Leu' : [1.07103448275862],
'Tyro' : [1.25793103448276],
'Thr' : [0.988275862068965],
'Phe' : [1.18137931034483],
'Adenine' : [1.49166523457639],
'Guanine' : [0.977315689981096],
'Cytosine' : [1.00446812167039],
'Thiamin' : [1.09957325746799],
'Nicotinic Acid' : [1.22901849217639],
'Pantothenate' : [0.995732574679943],
'Pyroxidal' : [1.06543385490754],
'Aminobenzoic acid' : [1.15149359886202],
'Folate' : [1.19985775248933],
'Zinc Sulfate' : [1.19132290184922],
'Copper Sulfate' : [1.0348506401138],
'Borate' : [1],
'Biotin' : [1.23157894736842],
'B12' : [1.23007518796992],
'Riboflavin' : [1.07293233082707],
'Mn SO4' : [1.10977443609023],
'Thymine' : [1.11278195488722],
'Sodium Chloride' : [1.09172932330827],
'Calcium Chloride' : [1.0812030075188],
'Iron Sulfate' : [1.19097744360902],}

media_colors = ['#6f9969', '#5c66a8', '#808fe1', '#454a74']

data_array = []
for comp in components:
    values = data_dict.get(comp, [])
    row = [values[i] if i < len(values) else np.nan for i in range(len(media))]
    data_array.append(row)
data_array = np.array(data_array)

fig, ax = plt.subplots(figsize=(8, len(components)*0.25))

cmap = plt.get_cmap('inferno')
im = ax.imshow(data_array, aspect='auto', cmap=cmap, interpolation='none')

for i in range(data_array.shape[0]):
    for j in range(data_array.shape[1]):
        val = data_array[i, j]
        if not np.isnan(val):
            norm_val = (val - np.nanmin(data_array)) / (np.nanmax(data_array) - np.nanmin(data_array))
            text_color = 'white' if norm_val < 0.5 else 'black'
            ax.text(j, i, f"{val:.2f}", ha='center', va='center', color=text_color, fontsize=12)


ax.set_xticks(np.arange(len(media)))
ax.set_xticklabels(media, fontsize=20)
ax.set_yticks(np.arange(len(components)))
ax.set_yticklabels(components, fontsize=15)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

for tick_label, color in zip(ax.get_xticklabels(), media_colors):
    tick_label.set_color(color)

plt.setp(ax.get_xticklabels(), ha='right')

for tick_label, color in zip(ax.get_xticklabels(), media_colors):
    tick_label.set_color(color)
    tick_label.set_horizontalalignment('center')

ax.set_xticks(np.arange(-0.5, len(media), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(components), 1), minor=True)
for x in range(len(media)+1):
    ax.axvline(x-0.5, color='gray', linewidth=4) 
for y in range(len(components)+1):
    ax.axhline(y-0.5, color='gray', linewidth=0)

ax.tick_params(which="minor", bottom=False, left=False)

# Colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Fold Change', fontsize=20)

plt.tight_layout()
plt.savefig('DO_heatmap.png', dpi=1200, bbox_inches='tight') 
plt.show()
