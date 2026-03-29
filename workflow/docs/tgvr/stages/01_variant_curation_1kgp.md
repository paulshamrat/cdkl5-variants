# 01 Variant Curation: 1KGP

This page documents the `1kgp` part of TGVR `01_variant_curation`.

There are two different ways TGVR can handle `1kgp`, and they should not be mixed up:

- `cached` or `auto` local use
  This is the normal stable local path.
- Colab-independent API-backed `1kgp`
  This is the safest path for users who do not have Palmetto and only need a Colab-friendly `1kgp` workbook build.
- raw Palmetto-backed `1kgp`
  This is the heavier legacy-style reproduction path and requires the Palmetto bridge first.

## 1KGP Modes

The main runner supports:

- `auto`
- `cached`
- `live`

`auto` is the default.
It prefers the cached workbook and only falls back to the live Ensembl-style path when needed.

Main command:

```bash
python workflow/tgvr/scripts/run_variant_curation.py GENE UNIPROT_ID --stage 1kgp --1kgp-mode auto
```

## Recommended Local Path

For normal local TGVR use, prefer the cached path:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached
```

For `CDKL5`, TGVR can seed the cached workbook from the legacy local export on first use.

This is the path to trust first when you want reproducible local work without rerunning the heavy raw workflow.

## Colab-Independent 1KGP

For users without Palmetto access, the recommended Colab path is to build the `1kgp` workbook directly from public Ensembl APIs.

This path is different from the raw Palmetto backend:

- it does not require the Palmetto bridge
- it does not require `bcftools`, `samtools`, or a local `VEP` cache
- it does not require downloading a full human offline `VEP` cache into Colab
- it produces the same practical TGVR-style `1kgp` workbook shape for `CDKL5`

What it does instead:

- looks up `CDKL5` from Ensembl
- queries Ensembl overlap for `1kg_3` variants in the gene region
- annotates missense `rs` ids through Ensembl `VEP`
- writes a workbook with the same useful sheets used by the cached TGVR path

Minimal Colab setup:

```python
import sys
import subprocess
from pathlib import Path

subprocess.run([sys.executable, "-m", "pip", "install", "-U", "pip"], check=True)
subprocess.run([sys.executable, "-m", "pip", "install", "pandas", "openpyxl", "requests"], check=True)

Path("/content/out").mkdir(parents=True, exist_ok=True)
print("Colab setup ready. Output dir: /content/out")
```

Practical recommendation:

- if you only need a usable `CDKL5` `1kgp` workbook in Colab, use the Colab-independent API-backed path
- if you need the heavier legacy-style `VCF -> VEP -> parsed table` reproduction, use Palmetto
- if you already trust the cached workbook, keep using `--1kgp-mode cached`

## Raw Palmetto-Backed 1KGP

The raw `1kgp` path is different.

This path is for the heavier `VCF -> VEP -> parsed table` style workflow and is currently implemented as a Palmetto-backed flow for `CDKL5`.

Important rule:

- before using the Palmetto-backed raw `1kgp` flow, establish the Palmetto bridge

## Palmetto Bridge

Run this exact block first:

```bash
rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock $USER@slogin.palmetto.clemson.edu "hostname && whoami"
```

If this bridge is not active, the Palmetto-backed `1kgp` manager commands should not be expected to work.

## Public Palmetto Manager

TGVR exposes one public manager for this flow:

- [manage_1kgp_palmetto.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/manage_1kgp_palmetto.py)

Available actions:

- `setup`
- `submit`
- `status`
- `log`
- `fetch`

## Full Palmetto Sequence

After the bridge is active, use this sequence:

```bash
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 setup
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 submit
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 status
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 log
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 fetch
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
  `workflow/tgvr/00_data/cdkl5/01_variant_curation/1kgp/palmetto_runs/grch38_allvar_noid_20260329_132508`
- rebuilt cached workbook:
  `workflow/tgvr/00_data/cdkl5/01_variant_curation/1kgp/1kgp_cdkl5_grch38.xlsx`

## Current Verified CDKL5 Counts

The memory file records the current verified `CDKL5` counts for the rebuilt cached workbook as:

- `4480` raw rows
- `19` missense rows
- `13` gene-only rows
- `12` unique rows
- `12` protein-change rows

The latest verified staged rerun after the fetch also produced:

- final manual-inclusive master total: `162`
- final kinase-domain total for `1-302`: `114`

Current verified final full-length label counts:

- `Benign` 22
- `Benign/Likely benign` 15
- `Conflicting classifications of pathogenicity` 7
- `Likely benign` 10
- `Likely pathogenic` 26
- `Pathogenic` 10
- `Pathogenic/Likely pathogenic` 24
- `Uncertain significance` 48

## Practical Recommendation

Use this decision rule:

- if you just want to run TGVR locally, use `--1kgp-mode cached`
- if you need a Colab-safe `1kgp` path without Palmetto, use the Colab-independent API-backed workbook build
- if you need to regenerate the heavy raw `1kgp` source path, bridge Palmetto first and use `manage_1kgp_palmetto.py`

## Canonical End-To-End Command Block

To avoid ambiguity later, keep the following exact block as the canonical begin-to-finish `CDKL5` variant-curation sequence:

```bash
cd ~/cdkl5-variants
conda activate tgvr

rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock $USER@slogin.palmetto.clemson.edu "hostname && whoami"

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
