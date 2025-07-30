#!/usr/bin/env python3
import os
import pandas as pd

# Base working directory
BASE_DIR = "/project/ealexov/compbio/shamrat/250519_energy/03_haddock"

# Complexes of interest
complexes = [
    "P49418_AMPH1_290-294",
    "Q86YP4_GATAD2A_97-101",
    "P48436_SOX9_197-202",
    "Q9P2Y4_ZNF219_111-115"
]

# Columns to extract from capri_ss.tsv
energy_cols = [
    'model', 'score', 'irmsd', 'fnat', 'lrmsd', 'dockq',
    'air', 'angles', 'bonds', 'bsa', 'cdihcoup',
    'dani', 'desolv', 'dihe', 'elec', 'improper',
    'rdcsrg', 'total', 'vdw', 'vean', 'xpcs'
]

# Terms to summarize
summary_terms = ['vdw', 'elec', 'desolv', 'air', 'total', 'score']

def extract_components():
    all_dfs = []
    for comp in complexes:
        ss_path = os.path.join(BASE_DIR, comp, "10_caprieval", "capri_ss.tsv")
        if not os.path.isfile(ss_path):
            print(f"⚠️ Missing {ss_path}, skipping {comp}")
            continue
        df = pd.read_csv(ss_path, sep=r"\s+", engine='python')
        df['complex'] = comp
        present = [c for c in energy_cols if c in df.columns]
        all_dfs.append(df[present + ['complex']])
    if not all_dfs:
        raise RuntimeError("No energy data found for any complex.")
    combined = pd.concat(all_dfs, ignore_index=True)
    cols = ['complex'] + [c for c in energy_cols if c in combined.columns]
    return combined[cols]

def save_components(combined):
    out_dir = os.path.join(BASE_DIR, "07_energies")
    os.makedirs(out_dir, exist_ok=True)
    csv_path   = os.path.join(out_dir, "haddock_energy_components.csv")
    excel_path = os.path.join(out_dir, "haddock_energy_components.xlsx")
    combined.to_csv(csv_path, index=False)
    combined.to_excel(excel_path, index=False)
    print(f"✅ Saved full components table to:\n  {csv_path}\n  {excel_path}")
    return excel_path

def summarize_components(components_file):
    df = pd.read_excel(components_file)
    summary = df.groupby('complex')[summary_terms].agg(['mean', 'std'])
    # flatten columns
    summary.columns = [f"{term}_{stat}" for term, stat in summary.columns]
    # Print summary to console
    print("\n===== Summary Statistics (Mean ± SD) =====")
    print(summary.round(2).to_string())
    # Save to Excel
    out_dir = os.path.join(BASE_DIR, "07_energies")
    summary_path = os.path.join(out_dir, "haddock_energy_summary.xlsx")
    summary.to_excel(summary_path)
    print(f"\n✅ Saved summary table to: {summary_path}")

def main():
    combined = extract_components()
    comp_file = save_components(combined)
    summarize_components(comp_file)

if __name__ == "__main__":
    main()

    
    

    
# plotting    
import pandas as pd
import matplotlib.pyplot as plt
import os

# User-configurable font size
FONT_SIZE = 16

# Load the summary table
base_dir = "/project/ealexov/compbio/shamrat/250519_energy/03_haddock/07_energies"
summary_path = os.path.join(base_dir, "haddock_energy_summary.xlsx")
df_summary = pd.read_excel(summary_path, index_col=0)

# Prepare data for plotting
terms = ['vdw', 'elec', 'desolv', 'air', 'total', 'score']
means = df_summary[[f"{t}_mean" for t in terms]]
stds = df_summary[[f"{t}_std" for t in terms]]

# Update global font size
plt.rcParams.update({'font.size': FONT_SIZE})

# Plot grouped bar chart
x = range(len(df_summary))
width = 0.15

plt.figure(figsize=(12, 7))
for i, term in enumerate(terms):
    plt.bar(
        [xi + i*width for xi in x],
        means[f"{term}_mean"],
        width=width,
        yerr=stds[f"{term}_std"],
        label=term,
        capsize=5
    )

# Tweak tick label size
plt.xticks([xi + width*2.5 for xi in x], df_summary.index, rotation=45, fontsize=FONT_SIZE)
plt.yticks(fontsize=FONT_SIZE)

# Axis labels and title with adjusted sizes
plt.ylabel("Energy (kcal/mol)", fontsize=FONT_SIZE)
plt.title("Mean ± SD of HADDOCK Energy Components by Complex", fontsize=FONT_SIZE + 2)

# Legend with adjusted text size
plt.legend(title="Component", fontsize=FONT_SIZE - 2, title_fontsize=FONT_SIZE)

plt.tight_layout()

# Save plot
plot_path = os.path.join(base_dir, "haddock_energy_components_plot.png")
plt.savefig(plot_path, dpi=300)
plt.show()
