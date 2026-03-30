# 01 Variant Curation

`01_variant_curation` is the most mature TGVR stage.

Its job is to build a master variant dataset for a target gene by combining:

- ClinVar
- 1000 Genomes Phase 3
- gnomAD allele frequency
- optional manual curated variants

It then performs a reference-sequence check and writes a final filtered master table.

## Main Runner

- `workflow/tgvr/scripts/run_variant_curation.py`

## Core Workflow

The master dataset is built in this order:

1. build the ClinVar core set
2. build the 1KGP missense set
3. keep only 1KGP variants that are unique versus ClinVar
4. load manual curated variants if available
5. keep only manual curated variants unique versus ClinVar plus 1KGP
6. use gnomAD as an AF lookup source
7. merge everything into one master table
8. run a sequence check against the configured reference FASTA
9. write mismatches separately and keep the final filtered table

## Per-Gene Layout

```text
workflow/tgvr/00_data/<gene>/01_variant_curation/
workflow/tgvr/outputs/<gene>/01_variant_curation/
```

On the first run for a new gene, TGVR also creates:

```text
workflow/tgvr/config/genes/<GENE>.yaml
workflow/tgvr/00_data/<gene>/reference/canonical_uniprot_<UNIPROT>.fasta
```

## Common Commands

Optional prerequisite when manual curated variants should be included:

```bash
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage master --manual-file /path/to/curated_variants.csv
```

Fastest run:

```bash
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage all
```

Recommended reproducible pattern when manual curated variants should be included:

```bash
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage all
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage master --manual-file /path/to/curated_variants.csv
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage master --use-manual
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage summary-full
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage summary-range --range START-END --label "Feature label"
```

Important:

- `--stage all` currently means the non-manual path
- if you want manual curated variants included, rerun `master` with `--use-manual`
- the manual curated file is a master-stage input, not a `1kgp`-stage input
- by default TGVR looks for it at `workflow/tgvr/00_data/<gene>/01_variant_curation/manual/curated_variants.csv`
- you can also install a manual file with `--manual-file <path>` before running the master stage
- the required minimal columns are `Mutation`, `Consequence`, and `Source`

## 1KGP Modes

TGVR currently supports:

- `cached`
- `vcf`
- `palmetto`

Mode summary:

- `cached`: reuse the local cached 1KGP workbook already stored under `workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/cache/`
- `vcf`: rebuild 1KGP from the public `phase3` chrX crossmap VCF without Palmetto; this is the main recommended public route
- `palmetto`: run the heavier legacy-style raw 1KGP reproduction path through Palmetto, then fetch and cache the rebuilt workbook locally

Current `CDKL5` comparison:

- `cached` and `palmetto` reproduce the legacy 1KGP benchmark exactly at the 1KGP stage: `4480 / 19 / 13 / 12 / 12`
- `vcf` is Palmetto-free and keeps the full legacy final variant core set, but currently expands the 1KGP-derived portion of the final table
- final full-length totals are:
  - legacy: `156`
  - cached/palmetto-backed current workflow: `162`
  - public VCF-backed current workflow: `163`

Comparison table:

| Route | 1KGP raw | 1KGP missense | 1KGP gene-only | 1KGP unique | 1KGP prot | Final total | Shared vs legacy final | Extra vs legacy final | Missing vs legacy final | Kinase total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Legacy | 4480 | 19 | 13 | 12 | 12 | 156 | 156 | 0 | 0 | 112 |
| Cached/Palmetto | 4480 | 19 | 13 | 12 | 12 | 162 | 156 | 6 | 0 | 114 |
| Public VCF | 4480 | 22 | 22 | 21 | 21 | 163 | 156 | 7 | 0 | 114 |

Full-length germline classification comparison:

| Germline classification | Legacy | Cached/Palmetto | Public VCF |
|---|---:|---:|---:|
| Benign | 20 | 22 | 23 |
| Benign/Likely benign | 15 | 15 | 15 |
| Conflicting classifications of pathogenicity | 13 | 7 | 7 |
| Likely benign | 10 | 10 | 10 |
| Likely pathogenic | 22 | 26 | 26 |
| Pathogenic | 9 | 10 | 10 |
| Pathogenic/Likely pathogenic | 24 | 24 | 24 |
| Uncertain significance | 43 | 48 | 48 |
| Total final variants | 156 | 162 | 163 |

Important:

- prefer `vcf` for the main public workflow and for most local reruns
- the `palmetto` route requires the Palmetto bridge first
- all three routes are accessed through `workflow/tgvr/scripts/run_variant_curation.py`
- see the dedicated `1KGP` page for the exact bridge block and the full `setup -> submit -> status -> log -> fetch` flow

Detailed guide:

- [01 Variant Curation: 1KGP](01_variant_curation_1kgp.md)

## CDKL5 Notes

Current documented TGVR choices for `CDKL5`:

- canonical UniProt accession: `O76039`
- displayed isoform: `O76039-2`
- sequence length: `960`
- working kinase-domain default: `1-302`

The older `13-297` range is historical context, but current TGVR documentation should use `1-302`.

## Stable User-Facing Outputs

The main stable outputs are:

- `10_master_dataset_final.csv`
- `11_full_length_classification_summary.csv`
- `range_summary_<label>_<start>_<end>.csv`
