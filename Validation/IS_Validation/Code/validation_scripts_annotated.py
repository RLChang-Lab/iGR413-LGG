import cobra
import pandas as pd
import matplotlib.pyplot as plt
from run_fba import fba
from set_dms import set_dm
import numpy as np
########################## NOTES ON VALIDATION SCRIPT ##########################
#   This script compares the fold changes in silico to fold change in vitro for single componenent dropout experiments
#   Each DM is done individually in this script
#   Takes the in vitro data from an excel sheet called iv_fc_EX.xlsx or sun_paper_fc.xlsx
#   Steps of analysis:
        #   load in model and set DM conditions
        #   Simulate dropout of each single component's effect on biomass 
        #   Convert raw flux to fold change of full media 
        #   Classify the in silico fold change change as no effect or deleterious based on 0.8 threshold
        #   Generate confusion matrix values (TN,TP, FN, FP) based on comparing to the in vitro fold changes. Calculate metrics. 
        #   Generate graphic depiction of confiusion matrix. 
#   Each formulation has some components that are not perfectly analogous between the in silico and in vitro conditions. these are listed and excluded from the confusion matrix. 

############# DM52 Validation #############
tmodelc = cobra.io.read_sbml_model('iGR413.xml')
tmodelc, media = set_dm(tmodelc, '52')
x = fba(tmodelc)  
print(f"Full growth rate: {x.objective_value}")
result = {}
is_fc = {}
for name, concentration in media.items():
    name.lower_bound = 0
    y = fba(tmodelc)
    result[name.id] = y.objective_value
    if x.objective_value < 0.43:
        is_fc[name.id] = y.objective_value / 0.43
    else:
        is_fc[name.id] = y.objective_value / x.objective_value
    name.lower_bound = concentration  
is_fc = pd.DataFrame.from_dict(is_fc, orient='index', columns=['IS'])
is_fc.index.name = 'exchange'
iv_data = pd.read_excel('iv_fc_EX.xlsx', na_values=["#N/A"])
iv_data = iv_data[['exchange', 'DM52']].rename(columns={'DM52': 'IV'})
iv_data.set_index('exchange', inplace=True)
merged_fcs = iv_data.join(is_fc, how='inner').reset_index()
excluded=['EX_inost_e','EX_ala_D_e','EX_lys_L_e','EX_gua_e','EX_thr_L_e','EX_pnto_R_e','EX_csn_e','EX_cu2_e','EX_pydx_e','EX_leu_L_e','EX_ribflv_e','EX_na1_e','EX_thm_e','EX_mn2_e','EX_thym_e','EX_phe_L_e','EX_fe2_e','EX_zn2_e','EX_fol_e','EX_nac_e','EX_cbl1_e','EX_btn_e','EX_gly_e','EX_trp_L_e','EX_tyr_L_e','EX_4abut_e','EX_adn_e']
TP = 0
FN = 0
TN = 0
FP = 0
matrix = {}
tp_list = []
fn_list = []
tn_list = []
fp_list = []

for _, results in merged_fcs.iterrows():
    exchange = results['exchange']
    
    IV_label = "no_effect" if results['IV'] > 0.8 else "deleterious"
    IS_label = "no_effect" if results['IS'] > 0.8 else "deleterious"
    
    matrix[exchange] = (IV_label, IS_label)
    
    if IV_label == "deleterious" and IS_label == "deleterious":
        TP += 1
        tp_list.append(exchange)
        
    elif IV_label == "deleterious" and IS_label == "no_effect":
        if exchange in excluded:
            TN += 1
            tn_list.append(exchange)
        else:
            FN += 1
            fn_list.append(exchange)
            
    elif IV_label == "no_effect" and IS_label == "no_effect":
        TN += 1
        tn_list.append(exchange)
        
    elif IV_label == "no_effect" and IS_label == "deleterious":
        FP += 1
        fp_list.append(exchange)


total = TP + TN + FP + FN
accuracy = (TP + TN) / total if total > 0 else 0
recall = TP / (TP + FN) if (TP + FN) > 0 else 0  # same as TPR
precision = TP / (TP + FP) if (TP + FP) > 0 else 0
specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
fpr = FP / (FP + TN) if (FP + TN) > 0 else 0

cm = np.array([[TN, FP],
               [FN, TP]])

metrics_text = (f"Accuracy={accuracy:.3f} | Recall={recall:.3f} | Precision={precision:.3f} | "
                f"Specificity={specificity:.3f} | FPR={fpr:.3f} ")

fig, ax = plt.subplots(figsize=(6,6))
im = ax.imshow(cm, cmap='Blues')
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(['Predicted No Effect', 'Predicted Deleterious'])
ax.set_yticklabels(['Actual No Effect', 'Actual Deleterious'])
ax.set_title("DM52 0.8 Threshold\n" + metrics_text, fontsize=14, fontweight='bold', pad=20)
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        color = "white" if cm[i, j] > cm.max() / 2 else "black"
        ax.text(j, i, cm[i, j], ha="center", va="center", color=color, fontsize=20)

plt.tight_layout()
plt.show()

############# DM24 Validation #############
model = cobra.io.read_sbml_model('iGR413.xml')
tmodelc, media = set_dm(tmodelc, '24')
x = fba(tmodelc)  # simulate full biomass 
result = {}
is_fc = {}
for name, concentration in media.items():
    name.lower_bound = 0
    y = fba(tmodelc)
    result[name.id] = y.objective_value
    if x.objective_value < 0.31:
        is_fc[name.id] = y.objective_value / 0.31
    else:
        is_fc[name.id] = y.objective_value / x.objective_value
    name.lower_bound = concentration  # restore original value
is_fc = pd.DataFrame.from_dict(is_fc, orient='index', columns=['IS'])
is_fc.index.name = 'exchange'
iv_data = pd.read_excel('iv_fc_EX.xlsx', na_values=["#N/A"])
iv_data = iv_data[['exchange', 'DM24']].rename(columns={'DM24': 'IV'})
iv_data.set_index('exchange', inplace=True)
merged_fcs = iv_data.join(is_fc, how='inner').reset_index()
DM24_excluded_exchanges=['EX_cit_e','EX_ura_e','EX_xan_e','EX_met_L_e']
TP = 0
FN = 0
TN = 0
FP = 0
matrix = {}
tp_list = []
fn_list = []
tn_list = []
fp_list = []

for _, results in merged_fcs.iterrows():
    exchange = results['exchange']
    IV_label = "no_effect" if results['IV'] > 0.8 else "deleterious"
    IS_label = "no_effect" if results['IS'] > 0.8 else "deleterious"
    matrix[exchange] = (IV_label, IS_label)
    if IV_label == "deleterious" and IS_label == "deleterious":
        TP += 1
        tp_list.append(exchange)
    elif IV_label == "deleterious" and IS_label == "no_effect":
        if exchange in DM24_excluded_exchanges:
            TN += 1
            tn_list.append(exchange)
        else:
            FN += 1
            fn_list.append(exchange)
    elif IV_label == "no_effect" and IS_label == "no_effect":
        TN += 1
        tn_list.append(exchange)
    elif IV_label == "no_effect" and IS_label == "deleterious":
        FP += 1
        fp_list.append(exchange)
total = TP + TN + FP + FN
accuracy = (TP + TN) / total if total > 0 else 0
recall = TP / (TP + FN) if (TP + FN) > 0 else 0  # same as TPR
precision = TP / (TP + FP) if (TP + FP) > 0 else 0
specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
fpr = FP / (FP + TN) if (FP + TN) > 0 else 0
cm = np.array([[TN, FP],[FN, TP]])
metrics_text = (f"Accuracy={accuracy:.3f} | Recall={recall:.3f} | Precision={precision:.3f} | "
                f"Specificity={specificity:.3f} | FPR={fpr:.3f} ")

fig, ax = plt.subplots(figsize=(6,6))
im = ax.imshow(cm, cmap='Blues')
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(['Predicted No Effect', 'Predicted Deleterious'])
ax.set_yticklabels(['Actual No Effect', 'Actual Deleterious'])
ax.set_title("DM24 0.8 Threshold\n" + metrics_text, fontsize=14, fontweight='bold', pad=20)
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        color = "white" if cm[i, j] > cm.max() / 2 else "black"
        ax.text(j, i, cm[i, j], ha="center", va="center", color=color, fontsize=20)

plt.tight_layout()
plt.show()

############# DM16 Validation #############
tmodelc = cobra.io.read_sbml_model('iGR413.xml')
tmodelc, media = set_dm(tmodelc, '16')
x = fba(tmodelc)  

result = {}
is_fc = {}
for name, concentration in media.items():
    name.lower_bound = 0
    y = fba(tmodelc)
    result[name.id] = y.objective_value
    if x.objective_value < 0.07:
        is_fc[name.id] = y.objective_value / 0.07
    else:
        is_fc[name.id] = y.objective_value / x.objective_value
    name.lower_bound = concentration 

is_fc = pd.DataFrame.from_dict(is_fc, orient='index', columns=['IS'])
is_fc.index.name = 'exchange'

iv_data = pd.read_excel('iv_fc_EX.xlsx', na_values=["#N/A"])
iv_data = iv_data[['exchange', 'dm16']].rename(columns={'dm16': 'IV'})
iv_data.set_index('exchange', inplace=True)
merged_fcs = iv_data.join(is_fc, how='inner').reset_index()

dm16_excluded_exchanges=['EX_cytd_e','EX_mops_e','EX_NH4_e']

TP = 0
FN = 0
TN = 0
FP = 0
matrix = {}
tp_list = []
fn_list = []
tn_list = []
fp_list = []

for _, results in merged_fcs.iterrows():
    exchange = results['exchange']
    IV_label = "no_effect" if results['IV'] > 0.8 else "deleterious"
    IS_label = "no_effect" if results['IS'] > 0.8 else "deleterious"
    matrix[exchange] = (IV_label, IS_label)
    if IV_label == "deleterious" and IS_label == "deleterious":
        TP += 1
        tp_list.append(exchange)
    elif IV_label == "deleterious" and IS_label == "no_effect":
        if exchange in dm16_excluded_exchanges:
            TN += 1
            tn_list.append(exchange)
        else:
            FN += 1
            fn_list.append(exchange) 
    elif IV_label == "no_effect" and IS_label == "no_effect":
        TN += 1
        tn_list.append(exchange)
    elif IV_label == "no_effect" and IS_label == "deleterious":
        FP += 1
        fp_list.append(exchange)

total = TP + TN + FP + FN
accuracy = (TP + TN) / total if total > 0 else 0
recall = TP / (TP + FN) if (TP + FN) > 0 else 0  # same as TPR
precision = TP / (TP + FP) if (TP + FP) > 0 else 0
specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
fpr = FP / (FP + TN) if (FP + TN) > 0 else 0

cm = np.array([[TN, FP],
               [FN, TP]])

metrics_text = (f"Accuracy={accuracy:.3f} | Recall={recall:.3f} | Precision={precision:.3f} | "
                f"Specificity={specificity:.3f} | FPR={fpr:.3f} ")

fig, ax = plt.subplots(figsize=(6,6))
im = ax.imshow(cm, cmap='Blues')

ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(['Predicted No Effect', 'Predicted Deleterious'])
ax.set_yticklabels(['Actual No Effect', 'Actual Deleterious'])
ax.set_title("DM16 0.8 Threshold\n" + metrics_text, fontsize=14, fontweight='bold', pad=20)

for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        color = "white" if cm[i, j] > cm.max() / 2 else "black"
        ax.text(j, i, cm[i, j], ha="center", va="center", color=color, fontsize=20)

plt.tight_layout()
plt.show()

############# DM13 Validation #############
tmodelc = cobra.io.read_sbml_model('iGR413.xml')
tmodelc, media = set_dm(tmodelc, '13')
x = fba(tmodelc)  

result = {}
is_fc = {}
for name, concentration in media.items():
    name.lower_bound = 0
    y = fba(tmodelc)
    result[name.id] = y.objective_value
    if x.objective_value < 0.11:
        is_fc[name.id] = y.objective_value / 0.11
    else:
        is_fc[name.id] = y.objective_value / x.objective_value
    name.lower_bound = concentration  

is_fc = pd.DataFrame.from_dict(is_fc, orient='index', columns=['IS'])
is_fc.index.name = 'exchange'
iv_data = pd.read_excel('iv_fc_EX.xlsx', na_values=["#N/A"])
iv_data = iv_data[['exchange', 'dm13']].rename(columns={'dm13': 'IV'})
iv_data.set_index('exchange', inplace=True)
merged_fcs = iv_data.join(is_fc, how='inner').reset_index()
TP=0
FN=0 
TN=0
FP=0
matrix={}
tp_list=[]
fn_list=[]
tn_list=[]
fp_list=[]
for _, results in merged_fcs.iterrows():
    exchange=results['exchange']
    if results['IV'] > 0.8: 
        IV_label="no_effect" 
    else: 
        IV_label="deleterious"
    if results['IS'] > 0.8: 
        IS_label="no_effect" 
    else: 
        IS_label="deleterious"
    matrix[results['exchange']]= IV_label, IS_label 
    if "deleterious" in IV_label and "deleterious" in IS_label: 
        TP=TP + 1
        tp_list.append(exchange)
    if "deleterious" in IV_label and "no_effect" in IS_label: 
        FN=FN + 1 
        fn_list.append(exchange)
    if "no_effect" in IV_label and "no_effect" in IS_label: 
        TN= TN + 1
        tn_list.append(exchange)
    if "no_effect" in IV_label and "deleterious" in IS_label: 
        FP= FP + 1
        fp_list.append(exchange)

total = TP + TN + FP + FN
accuracy = (TP + TN) / total if total > 0 else 0
recall = TP / (TP + FN) if (TP + FN) > 0 else 0  # same as TPR
precision = TP / (TP + FP) if (TP + FP) > 0 else 0
specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
fpr = FP / (FP + TN) if (FP + TN) > 0 else 0

cm = np.array([[TN, FP],
               [FN, TP]])
metrics_text = (f"Accuracy={accuracy:.3f} | Recall={recall:.3f} | Precision={precision:.3f} | "
                f"Specificity={specificity:.3f} | FPR={fpr:.3f} ")
fig, ax = plt.subplots(figsize=(6,6))
im = ax.imshow(cm, cmap='Blues')
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(['Predicted No Effect', 'Predicted Deleterious'])
ax.set_yticklabels(['Actual No Effect', 'Actual Deleterious'])
ax.set_title("DM13 0.8 Threshold\n" + metrics_text, fontsize=14, fontweight='bold', pad=20)

for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        color = "white" if cm[i, j] > cm.max() / 2 else "black"
        ax.text(j, i, cm[i, j], ha="center", va="center", color=color, fontsize=20)

plt.tight_layout()
plt.show()

############# Sun et al Validation #############
tmodelc = cobra.io.read_sbml_model('iGR413.xml')

media={tmodelc.reactions.EX_adn_e: -0.004,
tmodelc.reactions.EX_ala_D_e: -0.2,
tmodelc.reactions.EX_NH4_e: -1.2,
tmodelc.reactions.EX_arg_L_e: -0.4,
tmodelc.reactions.EX_asn_L_e: -0.44,
tmodelc.reactions.EX_asp_L_e: -0.32,
tmodelc.reactions.EX_btn_e: -0.002,
tmodelc.reactions.EX_cbl1_e: -0.002,
tmodelc.reactions.EX_cys_L_e: -0.6,
tmodelc.reactions.EX_k_e: -0.24,
tmodelc.reactions.EX_fol_e: -0.002,
tmodelc.reactions.EX_glc_D_e: -10,
tmodelc.reactions.EX_glu_L_e: -0.48,
tmodelc.reactions.EX_gln_L_e: -0.4,
tmodelc.reactions.EX_gly_e: -0.16,
tmodelc.reactions.EX_gua_e: -0.004,
tmodelc.reactions.EX_his_L_e: -0.176,
tmodelc.reactions.EX_fe2_e: -0.004,
tmodelc.reactions.EX_ile_L_e: -0.2,
tmodelc.reactions.EX_leu_L_e: -0.2,
tmodelc.reactions.EX_lys_L_e: -0.42,
tmodelc.reactions.EX_mg2_e: -0.08,
tmodelc.reactions.EX_mn2_e: -0.004,
tmodelc.reactions.EX_met_L_e: -0.08,
tmodelc.reactions.EX_k_e: -0.48,
tmodelc.reactions.EX_inost_e: -0.002,
tmodelc.reactions.EX_nac_e: -0.0016,
tmodelc.reactions.EX_pnto_R_e: -0.002,
tmodelc.reactions.EX_phe_L_e: -0.2,
tmodelc.reactions.EX_pro_L_e: -0.16,
tmodelc.reactions.EX_pydx_e: -0.0016,
tmodelc.reactions.EX_ribflv_e: -0.002,
tmodelc.reactions.EX_ser_L_e: -0.62,
tmodelc.reactions.EX_na1_e: -8,
tmodelc.reactions.EX_thm_e: -0.002,
tmodelc.reactions.EX_thr_L_e: -0.2,
tmodelc.reactions.EX_trp_L_e: -0.224,
tmodelc.reactions.EX_tyr_L_e: -0.16,
tmodelc.reactions.EX_ura_e: -0.004,
tmodelc.reactions.EX_val_L_e: -0.4,
tmodelc.reactions.EX_xan_e: -0.004,
tmodelc.reactions.EX_ac_e: -8,
tmodelc.reactions.EX_hco3_e: -1,
tmodelc.exchanges.EX_pi_e: -0.48, 
tmodelc.exchanges.EX_h_e: -.48} 
for q in tmodelc.reactions: 
    if "EX_" in q.id: 
        q.lower_bound=0
for name, concentration in media.items(): 
    name.lower_bound=concentration
print("sun_et_al_DM set")

x = fba(tmodelc)
result = {}
is_fc = {}
for name, concentration in media.items():
    name.lower_bound = 0
    y = fba(tmodelc)
    result[name.id] = y.objective_value
    if x.objective_value < 0.11:
        is_fc[name.id] = y.objective_value / 0.11
    else:
        is_fc[name.id] = y.objective_value / x.objective_value
    name.lower_bound = concentration  # restore original value

is_fc = pd.DataFrame.from_dict(is_fc, orient='index', columns=['IS'])
is_fc.index.name = 'exchange'
iv_data = pd.read_excel('sun_paper_fc.xlsx')
iv_data = iv_data[['exchange', 'FC-prelim']].rename(columns={'FC-prelim': 'IV'})
iv_data.set_index('exchange', inplace=True)
merged_fcs = iv_data.join(is_fc, how='inner').reset_index()
TP=0
FN=0 
TN=0
FP=0
matrix={}
tp_list=[]
fn_list=[]
tn_list=[]
fp_list=[]
for _, results in merged_fcs.iterrows():
    exchange=results['exchange']
    if results['IV'] > 0.8: 
        IV_label="no_effect" 
    else: 
        IV_label="deleterious"
    if results['IS'] > 0.8: 
        IS_label="no_effect" 
    else: 
        IS_label="deleterious"
    matrix[results['exchange']]= IV_label, IS_label 
    if "deleterious" in IV_label and "deleterious" in IS_label: 
        TP=TP + 1
        tp_list.append(exchange)
    if "deleterious" in IV_label and "no_effect" in IS_label: 
        FN=FN + 1 
        fn_list.append(exchange)
    if "no_effect" in IV_label and "no_effect" in IS_label: 
        TN= TN + 1
        tn_list.append(exchange)
    if "no_effect" in IV_label and "deleterious" in IS_label: 
        FP= FP + 1
        fp_list.append(exchange)

total = TP + TN + FP + FN
accuracy = (TP + TN) / total if total > 0 else 0
recall = TP / (TP + FN) if (TP + FN) > 0 else 0  
precision = TP / (TP + FP) if (TP + FP) > 0 else 0
specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
fpr = FP / (FP + TN) if (FP + TN) > 0 else 0
cm = np.array([[TN, FP],
               [FN, TP]])
metrics_text = (f"Accuracy={accuracy:.3f} | Recall={recall:.3f} | Precision={precision:.3f} | "
                f"Specificity={specificity:.3f} | FPR={fpr:.3f} ")
fig, ax = plt.subplots(figsize=(6,6))
im = ax.imshow(cm, cmap='Blues')
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(['Predicted No Effect', 'Predicted Deleterious'])
ax.set_yticklabels(['Actual No Effect', 'Actual Deleterious'])
ax.set_title("Sun 2019 Dropout 0.8 Threshold\n" + metrics_text, fontsize=14, fontweight='bold', pad=20)
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        color = "white" if cm[i, j] > cm.max() / 2 else "black"
        ax.text(j, i, cm[i, j], ha="center", va="center", color=color, fontsize=20)

plt.tight_layout()
plt.show()
