# 01 Variant Curation

`01_variant_curation` is the most mature TGVR stage.

Its job is to build a master variant dataset for a target gene by combining:

- ClinVar
- 1000 Genomes Phase 3
- gnomAD allele frequency
- optional manual curated variants

It then performs a reference-sequence check and writes a final filtered master table.

## Main Runner

- [run_variant_curation.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_variant_curation.py)

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

- `auto`
- `cached`
- `live`

`auto` is the default and prefers the cached workbook backend when available.

For `CDKL5`, this cached path is the preferred stable local route.

The heavier raw `1kgp` reproduction path is currently Palmetto-backed and exposed through:

- [manage_1kgp_palmetto.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/manage_1kgp_palmetto.py)

Important:

- the Palmetto-backed raw `1kgp` path requires the Palmetto bridge first
- see the dedicated `1KGP` page for the exact bridge block and the full `setup -> submit -> status -> log -> fetch` flow

Detailed guide:

- [01 Variant Curation: 1KGP](/home/paul/cdkl5-variants/workflow/docs/tgvr/stages/01_variant_curation_1kgp.md)

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
