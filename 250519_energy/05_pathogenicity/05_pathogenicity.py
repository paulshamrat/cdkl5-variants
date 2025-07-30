# 1) Polyphen2

#!/usr/bin/env python3
import pandas as pd
import subprocess

# ─── 0. Paths ────────────────────────────────────────────────────────────────
cdkl5_variants    = "/project/ealexov/compbio/shamrat/250519_energy/00_data/" \
                    "01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af_noddg.xlsx"
batch_out         = "/project/ealexov/compbio/shamrat/250519_energy/" \
                    "05_pathogenicity/01_polyphen2/cdkl5_mutation_polyphen2.txt"
tsv_in            = "/project/ealexov/compbio/shamrat/250519_energy/" \
                    "05_pathogenicity/01_polyphen2/" \
                    "cdkl5_mutation_polyphen2_results.tsv"
polyphen2_results = "/project/ealexov/compbio/shamrat/250519_energy/" \
                    "05_pathogenicity/01_polyphen2/" \
                    "cdkl5_mutation_polyphen2_results.xlsx"
merged_out        = "/project/ealexov/compbio/shamrat/250519_energy/" \
                    "05_pathogenicity/01_polyphen2/" \
                    "01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af_noddg_polyphen2.xlsx"

uniprot_id = "O76039"


# ─── 1.1 Prepare batch submission for webserver ──────────────────────────────
df = pd.read_excel(cdkl5_variants)
for col in ("wild", "position", "mutant"):
    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in {cdkl5_variants}")

lines = []
for _, row in df.iterrows():
    wt, pos, mt = row["wild"], row["position"], row["mutant"]
    if pd.isna(wt) or pd.isna(pos) or pd.isna(mt):
        continue
    lines.append(f"{uniprot_id} {int(pos)} {wt} {mt}")

with open(batch_out, "w") as fo:
    fo.write("\n".join(lines))
print(f"Wrote {len(lines)} variants to {batch_out}")


# ─── 1.2 Preview batch file head ────────────────────────────────────────────
subprocess.run(["head", batch_out])


# ─── 1.3 Convert returned TSV → Excel ────────────────────────────────────────
df_pp2 = pd.read_csv(tsv_in, sep="\t", engine="python")
df_pp2.columns = [c.strip().lstrip("#") for c in df_pp2.columns]
df_pp2.to_excel(polyphen2_results, index=False)
print(f"Wrote {len(df_pp2)} rows to {polyphen2_results}")


# ─── 1.4 Preview the new Excel & original columns ───────────────────────────
df_check = pd.read_excel(polyphen2_results)
print("PolyPhen-2 results head:\n", df_check.head().to_string(), "\n")
df_orig = pd.read_excel(cdkl5_variants)
print("cdkl5_variants columns:", df_orig.columns.tolist())


# ─── 1.5 Merge PolyPhen-2 into your variants (1:1 mapping) ───────────────────
variants_df = pd.read_excel(cdkl5_variants)
pp2_df       = pd.read_excel(polyphen2_results)

# Normalize key columns
variants_df['position'] = variants_df['position'].astype(int)
pp2_df['pos']          = pp2_df['pos'].astype(int)

# Select + rename
cols_to_pull = ['pos', 'prediction', 'pph2_prob', 'pph2_FPR', 'pph2_TPR']
pp2_small    = pp2_df[cols_to_pull].rename(columns={'pos': 'position'})

# Drop any duplicate positions (keep the first)
pp2_unique = pp2_small.drop_duplicates(subset='position', keep='first')

# Map back onto the original DataFrame
lookup = pp2_unique.set_index('position')
for col in ['prediction', 'pph2_prob', 'pph2_FPR', 'pph2_TPR']:
    variants_df[col] = variants_df['position'].map(lookup[col])

# Optional diagnostics
missing = variants_df['prediction'].isna().sum()
print(f"{missing} out of {len(variants_df)} variants had NO PolyPhen-2 result.")
extra_pp2 = set(pp2_unique['position']) - set(variants_df['position'])
print(f"{len(extra_pp2)} PolyPhen-2 positions didn’t match any original variants.")

# Write merged file
variants_df.to_excel(merged_out, index=False)
print(f"Wrote {len(variants_df)} rows (1:1 mapping) to {merged_out}")






#======================================================================================

# 02 MutPred2
## 2.1 Prepare fasta file with mutation for Mutpred2
#!/usr/bin/env python3
import os
import requests
import pandas as pd
from textwrap import wrap

# ─────── CONFIG ───────
UNIPROT_ID = "O76039"
OUT_DIR    = "/project/ealexov/compbio/shamrat/250519_energy/05_pathogenicity/02_mutpred2"
BASENAME   = "cdkl5_mutation_mutpred2"
XLSX_PATH  = "/project/ealexov/compbio/shamrat/250519_energy/00_data/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx"
WILD_COL   = "wild"
POS_COL    = "position"
MUT_COL    = "mutant"
SPECIES_TAG= "CDKL5_HUMAN"
MAX_PER    = 100
# ──────────────────────

# 1) fetch WT sequence from UniProt
url   = f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ID}.fasta"
resp  = requests.get(url)
resp.raise_for_status()
lines = resp.text.strip().splitlines()
seq   = "".join(lines[1:])

# 2) read your spreadsheet & build mutation tags
df   = pd.read_excel(XLSX_PATH, usecols=[WILD_COL,POS_COL,MUT_COL])
tags = [f"{wt}{int(pos)}{mt}" for wt,pos,mt in df.itertuples(index=False)]

# 3) chunk into ≤MAX_PER and write one FASTA per chunk
os.makedirs(OUT_DIR, exist_ok=True)
for i in range(0, len(tags), MAX_PER):
    chunk = tags[i : i+MAX_PER]
    idx   = i//MAX_PER + 1
    out_f = os.path.join(OUT_DIR, f"{BASENAME}_part{idx}.fasta")
    with open(out_f, "w") as fh:
        header = f">{SPECIES_TAG}_part{idx} " + " ".join(chunk)
        fh.write(header + "\n")
        for line in wrap(seq, 60):
            fh.write(line + "\n")
    print(f"→ wrote {out_f}")



## 2.2 Submit the fasta files and retrived the result

## 2.3 Read Mutpred2 part1 and part2 csv













# === 2.3a View the first few lines of your MutPred2 result CSVs ===S
print('=== 2.3a View the first few lines of your MutPred2 result CSVs ===')

import pandas as pd

# Updated paths with your new filenames
part1_csv = "/project/ealexov/compbio/shamrat/250519_energy/05_pathogenicity/02_mutpred2/" \
            "cdkl5_mutation_mutpred2_part1_result.csv"
part2_csv = "/project/ealexov/compbio/shamrat/250519_energy/05_pathogenicity/02_mutpred2/" \
            "cdkl5_mutation_mutpred2_part2_result.csv"

# Load and display the first 5 rows of part 1
df1 = pd.read_csv(part1_csv)
print(f"Head of {part1_csv}:")
print(df1.head().to_string(), "\n")

# Load and display the first 5 rows of part 2
df2 = pd.read_csv(part2_csv)
print(f"Head of {part2_csv}:")
print(df2.head().to_string())










# === 2.3b Combine your two MutPred2 CSVs into one table with proper headers ===

import pandas as pd
import os

# paths to your part-1 and part-2 result CSVs
dir02 = "/project/ealexov/compbio/shamrat/250519_energy/05_pathogenicity/02_mutpred2"
part1_csv = os.path.join(dir02, "cdkl5_mutation_mutpred2_part1_result.csv")
part2_csv = os.path.join(dir02, "cdkl5_mutation_mutpred2_part2_result.csv")

# 1) Define the column names in the correct order
column_names = [
    "ID",
    "Substitution",
    "MutPred2_score",
    "Molecular_mechanisms",
    "Affected_PROSITE_and_ELM_Motifs",
    "Remarks"
]

# 2) Read each CSV without a header row, assigning your names
df1 = pd.read_csv(part1_csv, header=None, names=column_names)
df2 = pd.read_csv(part2_csv, header=None, names=column_names)

# 3) Concatenate into a single DataFrame
combined = pd.concat([df1, df2], ignore_index=True)

# 4) (Optional) inspect the first few rows
print(combined.head().to_string())

# 5) Write out to a new Excel file, matching your filename style
out_xlsx = os.path.join(
    dir02,
    "cdkl5_mutation_mutpred2_results_combined.xlsx"
)
combined.to_excel(out_xlsx, index=False)
print(f"Wrote combined MutPred2 results to {out_xlsx}")












# === 2.3c Merge MutPred2 into your CDKL5 variants and save ===

import pandas as pd
import os

# paths
variants_xlsx = "/project/ealexov/compbio/shamrat/250519_energy/00_data/" \
                "01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af_noddg.xlsx"
mutpred2_dir  = "/project/ealexov/compbio/shamrat/250519_energy/05_pathogenicity/02_mutpred2"

# use the exact combined filename you generated
combined_xlsx = os.path.join(
    mutpred2_dir,
    "cdkl5_mutation_mutpred2_results_combined.xlsx"
)

# 1) Read your original variants and the combined MutPred2 results
variants_df = pd.read_excel(variants_xlsx, engine="openpyxl")
mutpred_df  = pd.read_excel(combined_xlsx, engine="openpyxl")

# 2) Merge on the variant string: 'mutation' in variants matches 'Substitution' in MutPred2
merged = variants_df.merge(
    mutpred_df,
    left_on="mutation",
    right_on="Substitution",
    how="left"
)

# 3) Write out to Excel, appending '_mutpred2' to the original base name
base = os.path.splitext(os.path.basename(variants_xlsx))[0]
out_path = os.path.join(
    mutpred2_dir,
    f"{base}_mutpred2.xlsx"
)
merged.to_excel(out_path, index=False)
print(f"Wrote merged MutPred2 annotations to {out_path}")









# 3 use esm code for esm prediction this need py11 version and 'haddock3_env' has python11 version
# the code is loaded in '03_esm/score_cdkl5_variants.py'




# === 3.2 Merge ESM-1v Δscore into your CDKL5 variants and save ===

import pandas as pd
import os

# 1) Paths to your original variants and the ESM results
variants_xlsx = "/project/ealexov/compbio/shamrat/250519_energy/00_data/" \
                "01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af_noddg.xlsx"
esm_folder    = "/project/ealexov/compbio/shamrat/250519_energy/05_pathogenicity/03_esm"
esm_results   = os.path.join(esm_folder, "cdkl5_esm1v_scores.xlsx")

# 2) Read in both tables
variants_df = pd.read_excel(variants_xlsx, engine="openpyxl")
esm_df      = pd.read_excel(esm_results,    engine="openpyxl")

# 3) Merge on the mutation string ('mutation' in variants ↔ 'mutation' in esm_df)
merged = variants_df.merge(
    esm_df,
    on="mutation",
    how="left"
)

# 4) Write out a new Excel with “_esm1v” appended
base = os.path.splitext(os.path.basename(variants_xlsx))[0]
out_path = os.path.join(
    esm_folder,
    f"{base}_esm1v.xlsx"
)
merged.to_excel(out_path, index=False)
print(f"Wrote merged ESM-1v annotations to {out_path}")












# === 08 Merge AlphaMissense results into your CDKL5 variants and save ===
## 08 alphamissense
#!/usr/bin/env python3
import pandas as pd
import os

# ─── Paths ───────────────────────────────────────────────────────────────────
variants_xlsx     = "/project/ealexov/compbio/shamrat/250519_energy/00_data/" \
                    "01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af_noddg.xlsx"
alphamissense_csv = "/project/ealexov/compbio/shamrat/250519_energy/" \
                    "05_pathogenicity/08_alphamissense/AF-O76039-F1-hg38.csv"
out_dir           = "/project/ealexov/compbio/shamrat/250519_energy/" \
                    "05_pathogenicity/08_alphamissense"

# ─── 1. Load data ─────────────────────────────────────────────────────────────
variants_df = pd.read_excel(variants_xlsx, engine="openpyxl")
am_df       = pd.read_csv(alphamissense_csv)

# ─── 2. Normalize merge keys ─────────────────────────────────────────────────
variants_df['mutation']        = variants_df['mutation'].astype(str)
am_df['protein_variant']       = am_df['protein_variant'].astype(str)

# ─── 3. Prepare AlphaMissense lookup ─────────────────────────────────────────
am_small = (
    am_df[['protein_variant', 'am_pathogenicity', 'am_class']]
    .rename(columns={'protein_variant': 'mutation'})
).drop_duplicates('mutation', keep='first')

# ─── 4. Map back onto variants_df ────────────────────────────────────────────
lookup = am_small.set_index('mutation')
variants_df['am_pathogenicity'] = variants_df['mutation'].map(lookup['am_pathogenicity'])
variants_df['am_class']         = variants_df['mutation'].map(lookup['am_class'])

# ─── 5. Quick confirmation ────────────────────────────────────────────────────
total   = len(variants_df)
matched = variants_df['am_pathogenicity'].notna().sum()
print(f"Pulled AlphaMissense results for {matched}/{total} variants.")

# Show a few example rows to confirm correct columns
print(variants_df[['mutation','am_pathogenicity','am_class']].head().to_string(index=False))

# ─── 6. Write out merged Excel ────────────────────────────────────────────────
base     = os.path.splitext(os.path.basename(variants_xlsx))[0]
out_path = os.path.join(out_dir, f"{base}_alphamissense.xlsx")
variants_df.to_excel(out_path, index=False)
print(f"Wrote {total} rows to {out_path}")
