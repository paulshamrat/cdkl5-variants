# TGVR Usage

TGVR is currently driven through explicit stage scripts under [workflow/tgvr/scripts](/home/paul/cdkl5-variants/workflow/tgvr/scripts).

## Current Public Entry Points

- [run_variant_curation.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_variant_curation.py)
- [run_folding.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_folding.py)

At the moment, there is no public `run_binding.py` script in the tree, so the binding stage should be treated as planned/documented work rather than a runnable public CLI.

Most stage commands take `GENE UNIPROT_ID`, for example:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage all
```

## 01 Variant Curation

This section is the main begin-to-finish workflow for `CDKL5` variant curation through one public script:

- [run_variant_curation.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_variant_curation.py)

Typical `CDKL5` variant-curation layout:

```text
workflow/tgvr
├── 00_data
│   └── cdkl5
│       ├── 01_variant_curation
│       │   └── manual
│       │       └── curated_variants.xlsx
│       └── reference
└── outputs
    └── cdkl5
        └── 01_variant_curation
            ├── 1kgp
            ├── clinvar
            ├── gnomad
            ├── logs
            └── master
```

Rule of thumb:

- `00_data/` keeps stable workflow inputs you want to preserve between runs
- `outputs/` keeps downloaded caches, generated stage files, fetched Palmetto runs, logs, and final results

### Optional Prerequisite: Manual Curated File

<details>
<summary><strong>Before you run with manual variants</strong></summary>

<p><code>master --use-manual</code> reads a curated file from <code>workflow/tgvr/00_data/cdkl5/01_variant_curation/manual/</code>.</p>

<p>Current workflow seed: <code>curated_variants.xlsx</code></p>

<p>Required columns: <code>Mutation</code>, <code>Consequence</code>, <code>Source</code></p>

<p>To replace that file, run one of these first:</p>

<pre><code class="language-bash">python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --manual-file /path/to/curated_variants.csv
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --manual-file /path/to/curated_variants.xlsx</code></pre>

</details>

### Recommended End-To-End Run

Use this command block as the main public TGVR path. It keeps the workflow fully local after download, uses the public VCF-backed `1kgp` route, and avoids the extra Palmetto bridge/setup burden:

```bash
cd ~/cdkl5-variants
conda activate tgvr

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage clinvar
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode vcf
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage gnomad
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --use-manual
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-full
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-range --range 1-302 --label "Kinase domain"
```

### 1KGP Mode Choices

Choose one of these three `1kgp` commands depending on the workflow you want.

<details open>
<summary><strong>Recommended Public Path: VCF</strong></summary>

<p>Use:</p>

<pre><code class="language-bash">python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode vcf</code></pre>

<p>Use this when:</p>

<ul>
  <li>you want the main public TGVR workflow</li>
  <li>you want a Palmetto-free run</li>
  <li>you want the comparison-friendly route highlighted in the docs</li>
</ul>

<p>What it does:</p>

<ul>
  <li>downloads or reuses the public <code>phase3</code> chrX crossmap VCF</li>
  <li>stores the source under <code>workflow/tgvr/outputs/&lt;gene&gt;/01_variant_curation/1kgp/source/</code></li>
  <li>subsets the <code>CDKL5</code> region locally with <code>bcftools</code></li>
  <li>annotates the subset through Ensembl</li>
  <li>writes standard TGVR <code>1kgp</code> stage outputs into <code>workflow/tgvr/outputs/&lt;gene&gt;/01_variant_curation/1kgp/</code></li>
</ul>

</details>

<details>
<summary><strong>Fast Local Rerun: Cached Workbook</strong></summary>

<p>Use:</p>

<pre><code class="language-bash">python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached</code></pre>

<p>Use this when:</p>

<ul>
  <li>you already trust the local cached workbook</li>
  <li>you want the fastest rerun</li>
  <li>you do not need to rebuild <code>1kgp</code></li>
</ul>

</details>

<details>
<summary><strong>Legacy Raw Backend: Palmetto</strong></summary>

<p>Before using Palmetto mode, establish the bridge first:</p>

<pre><code class="language-bash">rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock $USER@slogin.palmetto.clemson.edu "hostname && whoami"</code></pre>

<p>Then use:</p>

<pre><code class="language-bash">python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action setup
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action submit
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action status
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action log
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action fetch</code></pre>

<p>Use this when:</p>

<ul>
  <li>you explicitly want the legacy raw backend</li>
  <li>you need comparison against the legacy reconstruction path</li>
  <li>you are willing to do the extra bridge/setup flow</li>
</ul>

</details>

The `all` shortcut is still available, but the explicit staged form is preferred because the master depends on a correct `1kgp` stage.

Important note:

- `--stage all` currently means the non-manual path
- if you want manual curated variants included, rerun the `master` step explicitly with `--use-manual`
- the manual curated file is not part of the `1kgp` stage itself; it is read only by the `master --use-manual` step from `workflow/tgvr/00_data/<gene>/01_variant_curation/manual/curated_variants.csv` (or an installed `.xlsx` equivalent)

Runtime caches now live under `outputs`, for example:

- `workflow/tgvr/outputs/<gene>/01_variant_curation/clinvar/cache/clinvar_variant_summary.txt.gz`
- `workflow/tgvr/outputs/<gene>/01_variant_curation/gnomad/cache/raw_response.json`
- `workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/cache/1kgp_<gene>_grch38.xlsx`
- `workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/palmetto_runs/`

### Current CDKL5 Comparison

Current `1kgp` route comparison:

| Route | 1KGP raw | 1KGP missense | 1KGP gene-only | 1KGP unique | 1KGP prot | Final total | Shared vs legacy final | Extra vs legacy final | Missing vs legacy final | Kinase total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Legacy | 4480 | 19 | 13 | 12 | 12 | 156 | 156 | 0 | 0 | 112 |
| Cached/Palmetto | 4480 | 19 | 13 | 12 | 12 | 162 | 156 | 6 | 0 | 114 |
| Public VCF | 4480 | 22 | 22 | 21 | 21 | 163 | 156 | 7 | 0 | 114 |

Current full-length germline classification comparison:

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

Current kinase-domain (`1-302`) germline classification comparison:

| Germline classification | Legacy | Cached/Palmetto | Public VCF |
|---|---:|---:|---:|
| Benign | 4 | 5 | 5 |
| Benign/Likely benign | 12 | 12 | 12 |
| Conflicting classifications of pathogenicity | 10 | 5 | 5 |
| Likely benign | 1 | 1 | 1 |
| Likely pathogenic | 22 | 26 | 26 |
| Pathogenic | 9 | 9 | 9 |
| Pathogenic/Likely pathogenic | 24 | 24 | 24 |
| Uncertain significance | 30 | 32 | 32 |
| Total kinase-domain variants | 112 | 114 | 114 |

Variant-level comparison:

| Comparison | Shared | Left-only | Right-only |
|---|---:|---:|---:|
| Legacy vs Cached/Palmetto | 156 | 0 | 6 |
| Legacy vs Public VCF | 156 | 0 | 7 |
| Cached/Palmetto vs Public VCF | 162 | 0 | 1 |

Cached/Palmetto-only variants over legacy:

- `G546R`
- `M949I`
- `R59Q`
- `S493C`
- `T288R`
- `Y941C`

Public VCF-only variants over legacy:

- `G546R`
- `M949I`
- `R59Q`
- `R946G`
- `S493C`
- `T288R`
- `Y941C`

Only Public VCF adds beyond Cached/Palmetto:

- `R946G`

Detailed stage notes:

- [01 Variant Curation: 1KGP](/home/paul/cdkl5-variants/workflow/docs/tgvr/stages/01_variant_curation_1kgp.md)

## 02 Folding

Show stage status:

```bash
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage status
```

Prepare DDGun input:

```bash
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage prepare-ddgun --method ddgun_seq --range 1-302 --label "Kinase domain"
```
