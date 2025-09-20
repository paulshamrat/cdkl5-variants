# A) folding per method ddg with picked max
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- User parameters: adjust the Excel path if needed ---
folding_path = "/project/ealexov/compbio/shamrat/250519_energy/02_folding/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx"

# Load and filter to positions 1–302
df_fold = pd.read_excel(folding_path)
df_fold = df_fold[df_fold['position'].between(1, 302)]

# Identify per-method columns (exclude FoldX)
str_cols = [c for c in df_fold.columns if c.endswith('_str') and 'foldx' not in c.lower()]

# Compute the max absolute ΔΔG across methods
df_fold['ddG_Fmax'] = df_fold[str_cols].abs().max(axis=1)

# Split into benign and pathogenic
benign_df = df_fold[df_fold['Germline classification'] == 'Benign']
patho_df = df_fold[df_fold['Germline classification'] == 'Pathogenic']

# Build data dicts keyed by mutation, including position and each method's |ddG|
def build_data_dict(subdf):
    data = {}
    for _, row in subdf.iterrows():
        mut = row['mutation']
        entry = {'pos': row['position']}
        for m in str_cols:
            entry[m.replace('_str','')] = abs(row[m])
        # include picked max
        entry['max'] = row['ddG_Fmax']
        data[mut] = entry
    return data

benign_data = build_data_dict(benign_df)
patho_data  = build_data_dict(patho_df)

# Methods list (without '_str' suffix)
methods = [m.replace('_str','') for m in str_cols]
colors = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#B07AA1']
width = 0.15

# Plotting function matching the user style
def plot_side(ax, data, title):
    items = sorted(data.items(), key=lambda kv: kv[1]['pos'])
    variants = [k for k, _ in items]
    x = np.arange(len(variants))
    vals = np.array([[d[m] for m in methods] for _, d in items])
    max_vals = np.array([d['max'] for _, d in items])
    # grey overlay
    ax.bar(x, max_vals, width * len(methods), color='grey', alpha=0.3, zorder=0, label='Picked Max')
    # per-method bars
    for i, m in enumerate(methods):
        ax.bar(x + (i-(len(methods)-1)/2)*width, vals[:, i], width, color=colors[i], zorder=1,
               label=m if ax is plt.gca() and i==0 else "")
    # annotate max
    for xi, y in zip(x, max_vals):
        ax.text(xi, y + 0.05, f'{y:.2f}', ha='center', va='bottom')
    ax.set_xticks(x)
    ax.set_xticklabels(variants, rotation=45, ha='right')
    ax.set_title(title)
    ax.set_ylabel('|ddG|')

# Create figure
plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=False)

plot_side(axes[0], benign_data, 'Folding: Benign Variants')
plot_side(axes[1], patho_data,  'Folding: Pathogenic Variants')

# Legend centered above
handles = [plt.Rectangle((0,0),1,1,color='grey', alpha=0.3)] + [plt.Rectangle((0,0),1,1,color=c) for c in colors]
labels = ['Picked Max'] + methods
fig.suptitle('Folding: Per-Method |ddG| with Picked Max Overlay', fontsize=16, y=0.95)
fig.legend(handles, labels, loc='upper center', ncol=len(labels), frameon=False, fontsize=12)

plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.show()




## B) Folding ddg_Fmax for variants with threshold
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# === 1. Load & filter your folding dataset ===
folding_path = "/project/ealexov/compbio/shamrat/250519_energy/02_folding/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx"
df_fold = pd.read_excel(folding_path)
df_fold = df_fold[df_fold['position'].between(1, 302)]

# === 2. Identify per-method columns and compute overall max ΔΔG (ddG_Fmax) ===
str_cols = [c for c in df_fold.columns if c.endswith('_str') and 'foldx' not in c.lower()]
df_fold['ddG_Fmax'] = df_fold[str_cols].abs().max(axis=1)

# === 3. Subset to Benign vs Pathogenic variants ===
sub = df_fold[df_fold['Germline classification'].isin(['Benign', 'Pathogenic'])]
mutations = sub['mutation'].tolist()
ddG_vals   = sub['ddG_Fmax'].tolist()
groups     = sub['Germline classification'].tolist()

# === 4. Sort by ddG_Fmax ascending ===
sorted_idx    = np.argsort(ddG_vals)
mut_sorted    = [mutations[i] for i in sorted_idx]
vals_sorted   = [ddG_vals[i] for i in sorted_idx]
groups_sorted = [groups[i]   for i in sorted_idx]

# === 5. Compute threshold = (max Benign + min Pathogenic) / 2 ===
ben_max = sub[sub['Germline classification']=='Benign']['ddG_Fmax'].max()
pat_min = sub[sub['Germline classification']=='Pathogenic']['ddG_Fmax'].min()
threshold = (ben_max + pat_min) / 2

# === 6. Plotting ===
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(12, 6))

# Bar colors by group
color_map = {'Benign': '#0072B2', 'Pathogenic': '#D55E00'}
colors = [color_map[g] for g in groups_sorted]

# Draw bars
bars = ax.bar(mut_sorted, vals_sorted, color=colors)

# Draw threshold line
ax.axhline(threshold, color='gray', linestyle='--')

# Annotate each bar with its ddG_Fmax value
for bar, val in zip(bars, vals_sorted):
    ax.text(
        bar.get_x() + bar.get_width()/2,
        val + 0.05,
        f"{val:.3f}",
        ha='center',
        va='bottom'
    )

# Y‐axis ticks include the threshold
yticks = sorted(set(list(ax.get_yticks()) + [threshold]))
ax.set_yticks(yticks)
ax.set_yticklabels([
    f"{y:.1f}" if y != threshold else f"{threshold:.3f} (threshold)"
    for y in yticks
])

# Legend
legend_handles = [
    Patch(color=color_map['Benign'],    label='Benign'),
    Patch(color=color_map['Pathogenic'],label='Pathogenic'),
    Patch(color='gray',                 label='Threshold')
]
ax.legend(handles=legend_handles, loc='upper left')

# Labels & styling
ax.set_title('Folding ddG_Fmax for Variants with Threshold')
ax.set_xlabel('Mutation')
ax.set_ylabel('ddG_Fmax')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()








## C) Binding methods & Partners Averages (Benign vs Pathogenic)
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# — 0. Uniform font size —
plt.rcParams.update({'font.size': 14})

# — 1. Load & filter your binding data —
binding_path = "/project/ealexov/compbio/shamrat/250519_energy/04_binding/clinvar_1kgp_hector_gaf_final_binding_znf219_111_only.xlsx"
df_bind = pd.read_excel(binding_path)
df_bind = df_bind[df_bind['position'].between(1,302)].copy()

# — 2. Identify raw ddG columns and partners present —
ddg_cols = [c for c in df_bind.columns 
            if c.startswith('ddg_') and '_str_' in c and 'foldx' not in c.lower()]
all_partners = ['SOX9','AMPH1','GATAD2A','ZNF219']
partners = [p for p in all_partners if any(f"_{p}_" in c for c in ddg_cols)]

# — 3. Compute partner averages and overall max —
for p in partners:
    cols = [c for c in ddg_cols if f"_{p}_" in c]
    df_bind[f"avg_{p}"] = df_bind[cols].abs().mean(axis=1)
df_bind['ddG_Bmax'] = df_bind[[f"avg_{p}" for p in partners]].max(axis=1)

# — 4. Build sorted summary tables & variant lists —
benign_bind = (df_bind[df_bind['Germline classification']=='Benign']
               [['position','mutation','ddG_Bmax']]
               .sort_values('ddG_Bmax', ascending=False))
patho_bind = (df_bind[df_bind['Germline classification']=='Pathogenic']
               [['position','mutation','ddG_Bmax']]
               .sort_values('ddG_Bmax', ascending=True))

benign_idx = benign_bind.index
patho_idx  = patho_bind.index
benign_vars = benign_bind['mutation'].tolist()
patho_vars  = patho_bind['mutation'].tolist()

# — 5. Extract method suffixes (the part after "_str_") —
methods = sorted({ re.search(r'_str_(.+)$', c).group(1) 
                   for c in ddg_cols })

# — 6. Plotting colors & bar width —
colors    = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2']
avg_color = '#CCCCCC'
width     = 0.15

# — 7. Gather per‐method |ddG| values for benign & pathogenic —
benign_methods = {p: {} for p in partners}
patho_methods  = {p: {} for p in partners}
for p in partners:
    for m in methods:
        # find the one raw column matching this partner+method
        col = [c for c in ddg_cols if f"_{p}_" in c and c.endswith(f"_str_{m}")][0]
        benign_methods[p][m] = df_bind.loc[benign_idx, col].abs().values
        patho_methods[p][m]  = df_bind.loc[patho_idx,  col].abs().values

# — 8. Create a 2×N grid of plots (N = number of partners) —
fig, axes = plt.subplots(2, len(partners), figsize=(20, 10), sharey='row')
for row, (vars_list, methods_dict, grp_label) in enumerate([
    (benign_vars, benign_methods, 'Benign'),
    (patho_vars,  patho_methods,  'Pathogenic')
]):
    for col, p in enumerate(partners):
        ax = axes[row, col]
        x  = np.arange(len(vars_list))
        vals = np.vstack([methods_dict[p][m] for m in methods])
        avg  = vals.mean(axis=0)

        # faded average bar
        ax.bar(x, avg, width*len(methods), color=avg_color, alpha=0.3, zorder=0)
        # per-method bars
        for i, m in enumerate(methods):
            ax.bar(x + (i-(len(methods)-1)/2)*width,
                   vals[i], width, color=colors[i], zorder=1)
        # annotate the average
        for xi, yi in zip(x, avg):
            ax.text(xi, yi + np.nanmax(avg)*0.01,
                    f'{yi:.2f}', ha='center', va='bottom')

        ax.set_xticks(x)
        ax.set_xticklabels(vars_list, rotation=45, ha='right')
        ax.set_ylabel('|ΔΔG|')
        ax.set_title(f'{p} {grp_label}', fontsize=16, pad=10)
        ax.grid(axis='y', linestyle='--', alpha=0.5)

# — 9. Add a supertitle and shared legend —
fig.suptitle('Binding: Methods & Partner Averages — Benign vs Pathogenic', 
             fontsize=18, y=0.95)
handles = ([plt.Rectangle((0,0),1,1,color=avg_color, alpha=0.3)] +
           [plt.Rectangle((0,0),1,1,color=c) for c in colors])
labels  = ['Average'] + methods
fig.legend(handles, labels, loc='upper center', 
           bbox_to_anchor=(0.5,0.90), ncol=len(labels), frameon=False, fontsize=14)

plt.tight_layout(rect=[0,0,1,0.88])
plt.show()








## D) Binding ddG_Bmax for variants with threshold
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# 1. Load & filter the binding dataset
binding_path = "/project/ealexov/compbio/shamrat/250519_energy/04_binding/clinvar_1kgp_hector_gaf_final_binding_znf219_111_only.xlsx"
df = pd.read_excel(binding_path)
df = df[df['position'].between(1, 302)].copy()

# 2. Compute per-partner averages & overall ddG_Bmax if not already present
ddg_cols = [c for c in df.columns if c.startswith('ddg_') and '_str_' in c and 'foldx' not in c.lower()]
partners = ['SOX9','AMPH1','GATAD2A','ZNF219']
for p in partners:
    cols = [c for c in ddg_cols if f"_{p}_" in c]
    if cols:
        df[f"avg_{p}"] = df[cols].abs().mean(axis=1)
df['ddG_Bmax'] = df[[f"avg_{p}" for p in partners if f"avg_{p}" in df]].max(axis=1)

# 3. Build summary DataFrame and sort ascending by ddG_Bmax
ben = (df[df['Germline classification']=='Benign']
       [['mutation','ddG_Bmax']].assign(group='Benign'))
pat = (df[df['Germline classification']=='Pathogenic']
       [['mutation','ddG_Bmax']].assign(group='Pathogenic'))
summary = pd.concat([ben, pat], ignore_index=True)
summary = summary.sort_values('ddG_Bmax', ascending=True).reset_index(drop=True)

# 4. Extract aligned lists for plotting
mutations = summary['mutation'].tolist()
ddG_vals  = summary['ddG_Bmax'].tolist()
groups    = summary['group'].tolist()

# 5. Colors & threshold
color_map = {'Benign':'#0072B2', 'Pathogenic':'#D55E00'}
colors    = [color_map[g] for g in groups]
threshold = (ben['ddG_Bmax'].max() + pat['ddG_Bmax'].min()) / 2

# 6. Create the bar plot
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(12, 6))

bars = ax.bar(mutations, ddG_vals, color=colors)
ax.axhline(threshold, color='gray', linestyle='--')

# Annotate each bar with its value
for bar, val in zip(bars, ddG_vals):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.03,
            f'{val:.3f}', ha='center', va='bottom')

# Y-ticks including the threshold
yticks = sorted(set(list(ax.get_yticks()) + [threshold]))
ax.set_yticks(yticks)
ax.set_yticklabels([
    f'{y:.1f}' if y != threshold else f'{threshold:.3f} (threshold)'
    for y in yticks
])

# Legend
legend_handles = [
    Patch(color=color_map['Benign'],    label='Benign'),
    Patch(color=color_map['Pathogenic'],label='Pathogenic'),
    Patch(color='gray',                 label='Threshold')
]
ax.legend(handles=legend_handles, loc='upper left')

# Titles and labels
ax.set_title('Binding ddG_Bmax for Variants (Ascending) with Threshold')
ax.set_xlabel('Mutation')
ax.set_ylabel('ddG_Bmax')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()