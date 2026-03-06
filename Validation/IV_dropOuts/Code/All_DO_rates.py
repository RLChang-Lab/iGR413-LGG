############### All component Dropout Rates for supplemental figure ############ 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === CONFIG ===
xlsx_path = "All_rates.xlsx"  # path to your file
output_svg = "media_component_rates.svg"
palette = ["#7FC97F", "#BEAED4", "#FDC086", "#FFFF99"] 
formula_order = ["DM13", "DM16", "DM24", "DM52"]  

# === READ DATA ===
df = pd.read_excel(xlsx_path, sheet_name='Full_Data_Fig', index_col=0)
df = df.reset_index().rename(columns={"index": "Component"}) 
df.columns = df.columns.str.strip() 

# === DETERMINE RATE AND STAT COLUMNS ===
formula_cols = {}
for formula in formula_order:
    rate_cols = [c for c in df.columns if c.startswith(formula + "_") and "STAT" not in c]
    stat_cols = [c for c in df.columns if c.startswith(formula + "_STAT")]
    formula_cols[formula] = {"rates": rate_cols, "stat": stat_cols[0] if stat_cols else None}

# Optional: check that columns are detected correctly
for formula in formula_order:
    print(f"{formula} rate_cols:", formula_cols[formula]["rates"])
    print(f"{formula} stat_col:", formula_cols[formula]["stat"])

# === PREPARE LONG DATA ===
data_long = []
for formula, cols in formula_cols.items():
    for _, row in df.iterrows():
        comp = row["Component"]
        stat_p = row[cols["stat"]] if cols["stat"] else np.nan
        rates = [row[c] for c in cols["rates"] if not pd.isna(row[c])]
        for r in rates:
            data_long.append({"Component": comp, "Formula": formula, "Rate": r, "p_adj": stat_p})
data_long = pd.DataFrame(data_long)

# === PARSE p-values ===
def parse_pvalue(p):
    if pd.isna(p):
        return np.nan
    if isinstance(p, str):
        p = p.strip()
        if p.startswith("<") or p.startswith(">"):
            try:
                return float(p[1:])
            except ValueError:
                return np.nan
        else:
            try:
                return float(p)
            except ValueError:
                return np.nan
    return float(p)

data_long["p_adj"] = data_long["p_adj"].apply(parse_pvalue)

summary = (
    data_long.groupby(["Component", "Formula"])
    .agg(mean=("Rate", "mean"), std=("Rate", "std"), p_adj=("p_adj", "first"))
    .reset_index()
)

def p_to_stars(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

summary["stars"] = summary["p_adj"].apply(p_to_stars)

# === PLOT ===
plt.figure(figsize=(80, 20))  # very wide and tall
sns.set(style="whitegrid")

x = np.arange(len(df["Component"]))
bar_width = 0.2  # wider bars
offsets = np.linspace(-0.3, 0.3, len(formula_order))  # spread bars per group


for i, formula in enumerate(formula_order):
    sub = summary[summary["Formula"] == formula].set_index("Component")
    components_in_sub = sub.index.tolist()
    x_sub = [xi for xi, comp in enumerate(df["Component"]) if comp in components_in_sub]

    means = sub.loc[components_in_sub, "mean"]
    stds = sub.loc[components_in_sub, "std"]
    stars = sub.loc[components_in_sub, "stars"]


    for xi, comp in zip(x_sub, components_in_sub):
        plt.bar(
            xi + offsets[i],
            means.loc[comp],
            yerr=stds.loc[comp],
            width=bar_width,
            color=palette[i],
            capsize=3,
            edgecolor="black",
        )


        vals = data_long.query("Component == @comp and Formula == @formula")["Rate"]
        plt.scatter(
            np.repeat(xi + offsets[i], len(vals)),
            vals,
            color="black",
            s=10,
            alpha=0.7,
            zorder=3,
        )

        # Add stars above bar
        if stars.loc[comp] != "":
            plt.text(
                xi + offsets[i],
                means.loc[comp] + stds.loc[comp] + 0.05,
                stars.loc[comp],
                ha="center",
                va="bottom",
                fontsize=20,
                color="black",
            )

plt.xticks(x, df["Component"], rotation=45, ha="right")
plt.ylabel("Rate")
plt.title("Component rates by media formulation")
plt.xlim(x[0] - 0.5, x[-1] + 0.5)  # shrink left/right padding
plt.tight_layout()
plt.savefig(output_svg, format="svg", dpi=300)
plt.show()

print(f"Saved plot to {output_svg}")
