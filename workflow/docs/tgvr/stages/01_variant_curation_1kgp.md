# 01 Variant Curation: 1KGP

This page documents the `1kgp` part of TGVR `01_variant_curation`.

TGVR now uses one public entrypoint for all `1kgp` workflows:

- [run_variant_curation.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_variant_curation.py)

There are three supported `1kgp` modes under that one script:

- `cached`
  Reuse the existing local workbook under `workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/cache/`.
- `vcf`
  Build `1kgp` from the public `phase3` chrX crossmap VCF without Palmetto.
- `palmetto`
  Use the legacy Palmetto-backed raw workflow through the same main script.

## 1KGP Modes

The main runner supports three explicit `1kgp` commands:

Recommended public path:

```bash
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage 1kgp --1kgp-mode vcf
```

Trusted cached workbook path:

```bash
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage 1kgp --1kgp-mode cached
```

Legacy raw backend:

```bash
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action setup
```

## Recommended Local Path

For normal local TGVR use, prefer the public VCF-backed path:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode vcf
```

This is the main public-friendly route because it:

- uses a public source VCF
- avoids Palmetto bridge/setup steps
- still keeps the full legacy final variant core set for `CDKL5`

Use `cached` when you already have a trusted local workbook and want the fastest rerun. Use `palmetto` only when you explicitly need the legacy raw backend for comparison or regeneration.

## Public VCF-Backed 1KGP

The recommended public backend is:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode vcf
```

This mode:

- downloads or reuses the public `phase3` chrX crossmap VCF
- stores that source under `workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/source/`
- subsets the `CDKL5` region locally with `bcftools`
- annotates the region through Ensembl
- writes the standard TGVR `1kgp` stage outputs into `workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/`

Current verified local VCF-backed counts for `CDKL5` are:

- `4480` raw rows
- `22` missense rows
- `22` gene-only missense rows
- `21` unique missense rows
- `21` protein-change rows

Current comparison table:

| Route | Raw rows | Missense rows | Gene-only missense | Unique missense | Protein-change rows |
|---|---:|---:|---:|---:|---:|
| Legacy workbook | 4480 | 19 | 13 | 12 | 12 |
| Cached / fetched Palmetto workbook | 4480 | 19 | 13 | 12 | 12 |
| Public VCF-backed rebuild | 4480 | 22 | 22 | 21 | 21 |

Important comparison against the legacy workbook:

- all `12` legacy protein-change variants are present in the VCF-backed output
- the current VCF-backed run also includes `9` additional protein changes from the public source + current annotation path

## Raw Palmetto-Backed 1KGP

The raw `1kgp` path is different.

This path is for the heavier `VCF -> VEP -> parsed table` style workflow and is still available as a Palmetto-backed mode for `CDKL5`, but it is now considered an optional legacy/advanced path rather than the main documented route.

Before using the Palmetto-backed raw `1kgp` flow, establish the Palmetto bridge first.

## Palmetto Bridge

Run this exact block first:

```bash
rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock $USER@slogin.palmetto.clemson.edu "hostname && whoami"
```

If this bridge is not active, the Palmetto-backed `1kgp` commands should not be expected to work.

## Full Palmetto Sequence

After the bridge is active, use this sequence through the same main script:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action setup
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action submit
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action status
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action log
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action fetch
```

What these steps do:

- `setup`
  prepares TGVR-owned runtime assets under the Palmetto home workspace
- `submit`
  launches the raw `1kgp` job
- `status`
  checks the submitted job state
- `log`
  inspects the run log
- `fetch`
  pulls the completed result back into local TGVR and rebuilds the cached workbook

## After Fetch

Once the Palmetto result has been fetched, the normal local cached route should be used again:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached
```

This is the bridge between the heavy remote run and the normal local TGVR workflow.

## Current Verified Palmetto Rerun

The latest verified successful raw `1kgp` rerun used:

- remote asset root:
  `/home/YOUR_USERNAME/tgvr/resources/1kgp`
- remote script:
  `/home/YOUR_USERNAME/tgvr/jobs/1kgp_cdkl5_grch38_allvar_noid_<timestamp>.sh`
- remote run dir:
  `/home/YOUR_USERNAME/tgvr/outputs/cdkl5/01_variant_curation/1kgp/grch38_allvar_noid_<timestamp>`
- batch job id:
  `12138904`

The completed remote run produced:

- `cdkl5.GRCh38.all_variants_noid.tsv`
- `cdkl5.GRCh38.all_variants_noid.xlsx`
- `cdkl5.GRCh38.raw.af.csq.tsv`
- `cdkl5.GRCh38.vep.vcf.gz`
- `cdkl5.region.GRCh38.vcf.gz`
- `logs/tgvr_1kgp_cdkl5_grch38_allvar_noid_12138904.log`

After `fetch`, the local TGVR paths were:

- local fetched run dir:
  `workflow/tgvr/outputs/cdkl5/01_variant_curation/1kgp/palmetto_runs/grch38_allvar_noid_20260329_234509`
- rebuilt cached workbook:
  `workflow/tgvr/outputs/cdkl5/01_variant_curation/1kgp/cache/1kgp_cdkl5_grch38.xlsx`

## Current Verified CDKL5 Counts

The current verified legacy-style cached counts for the rebuilt workbook are:

- `4480` raw rows
- `19` missense rows
- `13` gene-only rows
- `12` unique rows
- `12` protein-change rows

## Practical Recommendation

Use this decision rule:

- if you want the main public TGVR workflow, use `--1kgp-mode vcf`
- if you already trust the local workbook and want the fastest rerun, use `--1kgp-mode cached`
- if you need to regenerate the heavy raw `1kgp` source path for legacy comparison, bridge Palmetto first and use `--1kgp-mode palmetto`

Current storage rule:

- keep stable workflow inputs in `workflow/tgvr/00_data/`
- keep generated `1kgp` caches, fetched Palmetto runs, and source VCF downloads in `workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/`

## Canonical End-To-End Command Block

To avoid ambiguity later, keep the following exact block as the canonical begin-to-finish Palmetto-backed `CDKL5` `1kgp` regeneration sequence:

```bash
cd ~/cdkl5-variants
conda activate tgvr

rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock $USER@slogin.palmetto.clemson.edu "hostname && whoami"

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action setup
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action submit
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action status
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action log
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action fetch

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage clinvar
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage gnomad
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --use-manual
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-full
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-range --range 1-302 --label "Kinase domain"
```
