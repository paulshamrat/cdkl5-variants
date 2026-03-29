import json
import re
import sys
from pathlib import Path

import pandas as pd


def assess_gene(gene_name):
    gene_upper = gene_name.upper()
    findings = {
        "gene": gene_upper,
        "sources": {},
        "takeaways": [],
    }

    clinvar_path = Path("00_data/clinvar_result.txt")
    if clinvar_path.exists():
        clinvar_df = pd.read_csv(clinvar_path, sep="\t", low_memory=False)
        gene_hits = 0
        if "Gene(s)" in clinvar_df.columns:
            gene_hits = int(clinvar_df["Gene(s)"].astype(str).str.contains(gene_upper, na=False).sum())
        findings["sources"]["clinvar"] = {
            "path": str(clinvar_path),
            "gene_filter_column": "Gene(s)" if "Gene(s)" in clinvar_df.columns else None,
            "gene_hits": gene_hits,
            "usable_from_gene_name_only": "Gene(s)" in clinvar_df.columns,
        }

    gnomad_path = Path("00_data/gnomAD.csv")
    if gnomad_path.exists():
        gnomad_df = pd.read_csv(gnomad_path, nrows=5)
        candidate_columns = [
            col
            for col in gnomad_df.columns
            if any(token in col.lower() for token in ["gene", "symbol", "transcript", "consequence"])
        ]
        direct_gene_column = next(
            (col for col in gnomad_df.columns if any(token in col.lower() for token in ["gene", "symbol"])),
            None,
        )
        findings["sources"]["gnomad"] = {
            "path": str(gnomad_path),
            "candidate_columns": candidate_columns,
            "direct_gene_column": direct_gene_column,
            "usable_from_gene_name_only": direct_gene_column is not None,
        }

    kgp_path = Path("00_data/1kgp_cdkl5_grch38.xlsx")
    if kgp_path.exists():
        kgp_df = pd.read_excel(kgp_path, sheet_name="Sheet1")
        gene_hits = 0
        if "Gene" in kgp_df.columns:
            gene_hits = int(kgp_df["Gene"].astype(str).str.contains(gene_upper, na=False).sum())
        findings["sources"]["kgp"] = {
            "path": str(kgp_path),
            "gene_filter_column": "Gene" if "Gene" in kgp_df.columns else None,
            "gene_hits": gene_hits,
            "usable_from_gene_name_only": "Gene" in kgp_df.columns,
        }

    hector_path = Path("00_data/hector2017.xlsx")
    if hector_path.exists():
        hector_df = pd.read_excel(hector_path, sheet_name="Sheet1")
        direct_gene_column = next(
            (col for col in hector_df.columns if any(token in col.lower() for token in ["gene", "symbol"])),
            None,
        )
        findings["sources"]["hector2017"] = {
            "path": str(hector_path),
            "columns": hector_df.columns.tolist(),
            "direct_gene_column": direct_gene_column,
            "usable_from_gene_name_only": direct_gene_column is not None,
        }

    findings["takeaways"] = [
        "ClinVar can be bootstrapped from only the gene name if the export contains a gene column.",
        "The local 1000 Genomes file can also be bootstrapped from only the gene name because it has a Gene column.",
        "The local gnomAD file does not expose a direct gene-symbol column, so gene name alone is not enough without an upstream gene-specific export or transcript mapping step.",
        "The local Hector file is already gene-specific and lacks a gene column, so it behaves like curated prior knowledge rather than a general-purpose searchable source.",
        "A truly gene-agnostic pipeline needs more than the gene symbol alone for robust automation; at minimum it should also resolve transcripts or a UniProt/reference-sequence mapping.",
    ]
    return findings


def save_live_clinvar_stages(gene_name, clinvar_bulk_path):
    gene_upper = gene_name.upper()
    output_dir = Path("workflow/tgvr/outputs") / gene_name.lower() / "01_variant_curation" / "clinvar"
    output_dir.mkdir(parents=True, exist_ok=True)

    usecols = ["GeneSymbol", "Type", "Name", "PhenotypeList", "ClinicalSignificance"]
    chunks = []
    for chunk in pd.read_csv(
        clinvar_bulk_path,
        sep="\t",
        compression="gzip",
        usecols=usecols,
        low_memory=False,
        chunksize=200000,
    ):
        sub = chunk[chunk["GeneSymbol"].astype(str).str.upper() == gene_upper]
        if not sub.empty:
            chunks.append(sub)

    gene_df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=usecols)
    snv_df = gene_df[gene_df["Type"].astype(str).str.contains("single nucleotide variant", case=False, na=False)].copy()
    missense_df = snv_df[
        snv_df["Name"].astype(str).str.contains(r"\(p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)", regex=True, na=False)
    ].copy()
    condition_df = missense_df[
        missense_df["PhenotypeList"].astype(str).str.contains(gene_upper, case=False, na=False)
    ].copy()
    final_df = condition_df.drop_duplicates().copy()

    stage_files = {
        "gene_rows": output_dir / "01_gene_rows.csv",
        "snv_rows": output_dir / "02_snv_rows.csv",
        "missense_rows": output_dir / "03_missense_rows.csv",
        "condition_rows": output_dir / "04_condition_rows.csv",
        "final_rows": output_dir / "05_final_deduplicated_rows.csv",
    }

    gene_df.to_csv(stage_files["gene_rows"], index=False)
    snv_df.to_csv(stage_files["snv_rows"], index=False)
    missense_df.to_csv(stage_files["missense_rows"], index=False)
    condition_df.to_csv(stage_files["condition_rows"], index=False)
    final_df.to_csv(stage_files["final_rows"], index=False)

    summary = {
        "gene": gene_upper,
        "source": str(clinvar_bulk_path),
        "filters": {
            "gene_symbol": gene_upper,
            "variant_type_contains": "single nucleotide variant",
            "missense_pattern": r"\(p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)",
            "phenotype_contains": gene_upper,
            "deduplicate": True,
        },
        "counts": {
            "gene_rows": int(len(gene_df)),
            "snv_rows": int(len(snv_df)),
            "missense_rows": int(len(missense_df)),
            "condition_rows": int(len(condition_df)),
            "final_rows": int(len(final_df)),
        },
        "files": {key: str(path) for key, path in stage_files.items()},
    }

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def write_markdown_report(findings, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# Local Gene-Name Bootstrap Findings: {findings['gene']}",
        "",
        "## Question",
        f"What can be discovered from the local data sources if the only provided input is the gene name `{findings['gene']}`?",
        "",
        "## Source Assessment",
    ]

    for source_name, source_info in findings["sources"].items():
        lines.extend(
            [
                f"### {source_name}",
                f"- path: `{source_info.get('path')}`",
                f"- usable from gene name only: `{source_info.get('usable_from_gene_name_only')}`",
            ]
        )
        if "gene_filter_column" in source_info:
            lines.append(f"- gene filter column: `{source_info.get('gene_filter_column')}`")
        if "gene_hits" in source_info:
            lines.append(f"- matching rows for gene: `{source_info.get('gene_hits')}`")
        if "candidate_columns" in source_info:
            lines.append(f"- candidate columns: `{source_info.get('candidate_columns')}`")
        if "columns" in source_info:
            lines.append(f"- columns: `{source_info.get('columns')}`")
        lines.append("")

    lines.extend(["## Takeaways"])
    lines.extend([f"- {item}" for item in findings["takeaways"]])
    lines.extend(
        [
            "",
            "## Practical Conclusion",
            f"For `{findings['gene']}`, a gene-name-only bootstrap can partially work on the local data, but it is not sufficient for a robust general pipeline across all sources. Additional identifiers or source-specific resolution steps are still needed.",
        ]
    )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python workflow/tgvr/assess_local_gene_bootstrap.py <GENE>")

    gene_name = sys.argv[1]
    findings = assess_gene(gene_name)
    report_path = (
        Path("workflow/tgvr/reports")
        / f"{gene_name.lower()}_local_findings.md"
    )
    write_markdown_report(findings, report_path)

    clinvar_bulk_path = Path("/tmp/clinvar_variant_summary.txt.gz")
    if clinvar_bulk_path.exists():
        live_summary = save_live_clinvar_stages(gene_name, clinvar_bulk_path)
        print(
            "Saved live ClinVar stage outputs:",
            live_summary["counts"],
        )

    print(f"Assessment completed for gene: {findings['gene']}")
    print(f"Report written to: {report_path}")


if __name__ == "__main__":
    main()
