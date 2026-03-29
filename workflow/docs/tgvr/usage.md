# TGVR Usage

TGVR is currently driven through explicit stage scripts under [workflow/tgvr/scripts](/home/paul/cdkl5-variants/workflow/tgvr/scripts).

## Current Public Entry Points

- [run_variant_curation.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_variant_curation.py)
- [run_folding.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_folding.py)
- [manage_1kgp_palmetto.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/manage_1kgp_palmetto.py)

At the moment, there is no public `run_binding.py` script in the tree, so the binding stage should be treated as planned/documented work rather than a runnable public CLI.

## General Pattern

Most TGVR commands take:

```text
GENE UNIPROT_ID
```

Example:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage all
```

## Typical Variant-Curation Flow

The canonical begin-to-finish `CDKL5` variant-curation workflow is the following sequence.

This is the preferred command block to keep in mind when you want the full reproducible path from Palmetto-backed raw `1kgp` regeneration through the final TGVR summaries:

```bash
cd ~/cdkl5-variants
conda activate tgvr

rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN YOUR_USERNAME@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check YOUR_USERNAME@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock YOUR_USERNAME@slogin.palmetto.clemson.edu "hostname && whoami"

python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 setup
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 submit
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 status
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 log
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 fetch

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage clinvar
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage gnomad
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --use-manual
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-full
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-range --range 1-302 --label "Kinase domain"
```

If the cached `1kgp` workbook is already trusted and does not need regeneration, you can skip the Palmetto block and start from the local TGVR stage commands.
The `all` shortcut is still available, but the explicit staged form above is preferred because the master depends on a correct `1kgp` stage and the raw Palmetto-backed regeneration path has its own required setup.

Important note:

- `--stage all` currently means the non-manual path
- if you want manual curated variants included, rerun the `master` step explicitly with `--use-manual`
- for the heavier raw Palmetto-backed `1kgp` path, do not start with `manage_1kgp_palmetto.py` until the Palmetto bridge is active

## Typical Folding Prep Flow

Show stage status:

```bash
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage status
```

Prepare DDGun input:

```bash
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage prepare-ddgun --method ddgun_seq --range 1-302 --label "Kinase domain"
```

## Palmetto 1KGP Flow

For the heavier raw `1kgp` reproduction path, TGVR currently provides one public manager.

Use this when you need to regenerate the cached `1kgp` workbook from the raw Palmetto-backed path before the local TGVR master build.

First establish the Palmetto bridge:

```bash
rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN YOUR_USERNAME@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check YOUR_USERNAME@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock YOUR_USERNAME@slogin.palmetto.clemson.edu "hostname && whoami"
```

Then run:

```bash
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 setup
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 submit
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 status
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 log
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 fetch
```

This flow is currently specialized around the `CDKL5` raw `1kgp` path.

Detailed stage notes:

- [01 Variant Curation: 1KGP](/home/paul/cdkl5-variants/workflow/docs/tgvr/stages/01_variant_curation_1kgp.md)

## Current Verified CDKL5 Rerun

The current verified staged rerun for `CDKL5` uses exactly this sequence:

```bash
cd ~/cdkl5-variants
conda activate tgvr

rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN YOUR_USERNAME@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check YOUR_USERNAME@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock YOUR_USERNAME@slogin.palmetto.clemson.edu "hostname && whoami"

python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 setup
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 submit
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 status
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 log
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 fetch

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage clinvar
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage gnomad
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --use-manual
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-full
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-range --range 1-302 --label "Kinase domain"
```

Verified current outputs:

- full-length total: `162`
- kinase-domain total (`1-302`): `114`
- full-length counts:
  - `Benign` 22
  - `Benign/Likely benign` 15
  - `Conflicting classifications of pathogenicity` 7
  - `Likely benign` 10
  - `Likely pathogenic` 26
  - `Pathogenic` 10
  - `Pathogenic/Likely pathogenic` 24
  - `Uncertain significance` 48
