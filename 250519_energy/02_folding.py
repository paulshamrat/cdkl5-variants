import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
from matplotlib.patches import Patch

# ——— Configuration ———
file_path           = "00_data/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx"
output_dir          = "02_folding/10_output/250522_seq_str"
fig_width           = 12
fig_height          = 10
font_size           = 12
panel_letter_size   = 16
title_fontsize      = 16
legend_fontsize     = 16   # legend text size
hspace              = 0.6   # vertical space between rows
wspace              = 0.5   # horizontal space between columns
pad_inches          = 0.1
rect_params         = [0, 0.03, 1, 0.94]  # leave more room at top for legend

os.makedirs(output_dir, exist_ok=True)

# ——— Helper function ———
def add_germline_legend(fig, palette, fontsize, y=1.02):
    header = Patch(facecolor='none', edgecolor='none',
                   label='Germline classification:')
    ben_handle  = Patch(facecolor=palette['Benign'], edgecolor='black',
                        label='Benign')
    path_handle = Patch(facecolor=palette['Pathogenic'], edgecolor='black',
                        label='Pathogenic')

    fig.legend(
        handles=[header, ben_handle, path_handle],
        labels=[h.get_label() for h in (header, ben_handle, path_handle)],
        loc='upper center',
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, y),
        fontsize=fontsize,
        handlelength=1
    )

# Load data
df = pd.read_excel(file_path, sheet_name=0)

# Identify sequence- and structure-based ΔΔG columns
ddg_seq_cols = [c for c in df.columns if re.match(r'(?i)^ddg.*seq$', c)]
ddg_str_cols = [c for c in df.columns if re.match(r'(?i)^ddg.*str$', c)]

# Compute global y-limits
all_vals = df[ddg_seq_cols + ddg_str_cols].to_numpy().flatten()
all_vals = all_vals[~pd.isna(all_vals)]
y_min, y_max = all_vals.min(), all_vals.max()

# Classes + palette
classes = ['Benign', 'Pathogenic']
palette = {'Benign': 'tab:blue', 'Pathogenic': 'tab:orange'}

# Apply base font size
plt.rcParams.update({'font.size': font_size})

# Create 2×2 grid
fig, axes = plt.subplots(
    nrows=2, ncols=2,
    figsize=(fig_width, fig_height),
    sharey=True
)
fig.subplots_adjust(hspace=hspace, wspace=wspace)

# Panel definitions
panels = [
    (axes[0,0], ddg_seq_cols, "Folding (Sequence): Positions 1–960", 1, 895, 'A'),
    (axes[0,1], ddg_seq_cols, "Folding (Sequence): Positions 1–302", 1, 302, 'B'),
    (axes[1,0], ddg_str_cols, "Folding (Structure): Positions 1–960", 1, 895, 'C'),
    (axes[1,1], ddg_str_cols, "Folding (Structure): Positions 1–302", 1, 302, 'D'),
]

for ax, cols, title, lo, hi, letter in panels:
    sub = df[df.position.between(lo, hi) &
             df['Germline classification'].isin(classes)]
    counts    = sub['Germline classification'].value_counts()
    count_str = " | ".join(f"{c}: {counts.get(c,0)}"
                           for c in classes)
    long = sub.melt(
        id_vars=['Germline classification'],
        value_vars=cols,
        var_name='Method',
        value_name='ΔΔG'
    )
    sns.violinplot(
        data=long, x='Method', y='ΔΔG',
        hue='Germline classification',
        palette=palette, inner='box',
        cut=0, width=0.4, dodge=True,
        ax=ax
    )
    # Subplot letter
    ax.text(-0.1, 1.05, letter,
            transform=ax.transAxes,
            fontsize=panel_letter_size,
            fontweight='bold')
    # Title
    ax.set_title(f"{title}\nCounts: {count_str}",
                 pad=16, fontsize=title_fontsize)
    # X-label on every panel now
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.tick_params(axis='x', rotation=45)
    ax.set_ylim(y_min, y_max)
    if ax.get_legend():
        ax.get_legend().remove()

# Add the Germline classification legend
add_germline_legend(fig, palette, legend_fontsize, y=1.02)

# Final layout and save
plt.tight_layout(rect=rect_params)
print(f"Figure size (inches): {fig.get_size_inches()}")
print(f"Base font size: {plt.rcParams['font.size']}")

out_path = os.path.join(output_dir,
                        "folding_sequence_structure_composite.png")
fig.savefig(out_path,
            dpi=300,
            bbox_inches='tight',
            pad_inches=pad_inches)
plt.show()
plt.close(fig)
