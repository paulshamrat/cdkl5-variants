import json
import os
import re
from collections import OrderedDict

import pandas as pd

from src.data_utils import add_mutation_components, get_absolute_path, load_config


AA_3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Pyl": "O",
    "Ter": "*",
    "Asx": "B",
    "Glx": "Z",
    "Xaa": "X",
}

COMBINED_COLUMNS = [
    "Name",
    "Gene(s)",
    "Protein change",
    "Condition(s)",
    "Accession",
    "GRCh37Chromosome",
    "GRCh37Location",
    "GRCh38Chromosome",
    "GRCh38Location",
    "VariationID",
    "AlleleID(s)",
    "dbSNP ID",
    "Canonical SPDI",
    "Variant type",
    "Molecular consequence",
    "Germline classification",
    "Source",
    "Germline date last evaluated",
    "Germline review status",
    "Somatic clinical impact",
    "Somatic clinical impact date last evaluated",
    "Somatic clinical impact review status",
    "Oncogenicity classification",
    "Oncogenicity date last evaluated",
    "Oncogenicity review status",
    "Unnamed: 24",
]


def convert_three_letter_change(protein_change):
    """Converts protein changes like Ile3Phe to I3F."""
    if not isinstance(protein_change, str):
        return None

    match = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)", protein_change)
    if not match:
        return None

    wild, position, mutant = match.groups()
    wild_one = AA_3_TO_1.get(wild)
    mutant_one = AA_3_TO_1.get(mutant)
    if not wild_one or not mutant_one:
        return None
    return f"{wild_one}{position}{mutant_one}"


def read_fasta_sequence(fasta_path):
    """Reads the first FASTA sequence from a file."""
    sequence_lines = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            sequence_lines.append(line)
    return "".join(sequence_lines)


def build_clinvar_frames(raw_clinvar_txt):
    raw_df = pd.read_csv(raw_clinvar_txt, sep="\t", low_memory=False)
    single_nucleotide_df = raw_df[raw_df["Variant type"] == "single nucleotide variant"].copy()
    missense_df = single_nucleotide_df[
        single_nucleotide_df["Molecular consequence"].astype(str).str.contains("missense variant", na=False)
    ].copy()
    cdkl5_condition_df = missense_df[
        missense_df["Condition(s)"].astype(str).str.contains("CDKL5", na=False)
    ].copy()
    cdkl5_missense_only_df = cdkl5_condition_df[
        ~cdkl5_condition_df["Molecular consequence"].astype(str).str.contains("intron variant", na=False)
    ].copy()

    return OrderedDict(
        [
            ("Sheet1", raw_df),
            ("Single_Nucleotide_Variant", single_nucleotide_df),
            ("Missense_Variant", missense_df),
            ("CDKL5_Condition", cdkl5_condition_df),
            ("CDKL5_missense_only", cdkl5_missense_only_df),
        ]
    )


def build_gnomad_frames(gnomad_csv):
    raw_df = pd.read_csv(gnomad_csv)
    missense_df = raw_df[raw_df["VEP Annotation"].astype(str).str.contains("missense_variant", na=False)].copy()
    converted_df = missense_df.copy()
    converted_df["Protein Consequence (Three-Letter)"] = converted_df["Protein Consequence"].astype(str).str.replace(
        "p.", "", regex=False
    )
    converted_df["Protein change (One-Letter)"] = converted_df["Protein Consequence (Three-Letter)"].apply(
        convert_three_letter_change
    )
    clinvar_classified_only_df = converted_df[converted_df["ClinVar Germline Classification"].notna()].copy()

    return OrderedDict(
        [
            ("Sheet1", raw_df),
            ("Missense_Variant", missense_df),
            ("Protein_Change_Converted", converted_df),
            ("ClinVar_Classified_Only", clinvar_classified_only_df),
        ]
    )


def build_kgp_frames(kgp_source_workbook):
    raw_df = pd.read_excel(kgp_source_workbook, sheet_name="Sheet1")
    missense_only_df = raw_df[raw_df["Consequence"] == "missense_variant"].copy()
    missense_cdkl5_only_df = missense_only_df[missense_only_df["Gene"] == "CDKL5"].copy()
    missense_cdkl5_unique_df = missense_cdkl5_only_df.drop_duplicates().copy()
    missense_cdkl5_unique_prot_df = missense_cdkl5_unique_df.copy()
    missense_cdkl5_unique_prot_df["Protein change 3L"] = missense_cdkl5_unique_prot_df["HGVSp"].astype(str).str.extract(
        r":p\.(.+)"
    )
    missense_cdkl5_unique_prot_df["Protein change"] = missense_cdkl5_unique_prot_df["Protein change 3L"].apply(
        convert_three_letter_change
    )

    return OrderedDict(
        [
            ("Sheet1", raw_df),
            ("missense_only", missense_only_df),
            ("missense_cdkl5_only", missense_cdkl5_only_df),
            ("missense_cdkl5_unique", missense_cdkl5_unique_df),
            ("missense_cdkl5_unique_prot_chan", missense_cdkl5_unique_prot_df),
        ]
    )


def standardize_clinvar_rows(clinvar_df):
    standardized = clinvar_df.copy()
    standardized["Source"] = "Clinvar"
    return standardized.reindex(columns=COMBINED_COLUMNS)


def standardize_kgp_rows(kgp_unique_df):
    return pd.DataFrame(
        {
            "Name": kgp_unique_df["HGVSp"],
            "Gene(s)": "CDKL5",
            "Protein change": kgp_unique_df["Protein change"],
            "Condition(s)": "Healthy",
            "Accession": pd.NA,
            "GRCh37Chromosome": pd.NA,
            "GRCh37Location": pd.NA,
            "GRCh38Chromosome": kgp_unique_df["CHROM"].astype(str).str.replace("chr", "", regex=False),
            "GRCh38Location": kgp_unique_df["POS"],
            "VariationID": pd.NA,
            "AlleleID(s)": pd.NA,
            "dbSNP ID": pd.NA,
            "Canonical SPDI": pd.NA,
            "Variant type": pd.NA,
            "Molecular consequence": "missense variant",
            "Germline classification": "Benign",
            "Source": "The 1000 Genomes Project",
            "Germline date last evaluated": pd.NA,
            "Germline review status": pd.NA,
            "Somatic clinical impact": pd.NA,
            "Somatic clinical impact date last evaluated": pd.NA,
            "Somatic clinical impact review status": pd.NA,
            "Oncogenicity classification": pd.NA,
            "Oncogenicity date last evaluated": pd.NA,
            "Oncogenicity review status": pd.NA,
            "Unnamed: 24": pd.NA,
        }
    ).reindex(columns=COMBINED_COLUMNS)


def standardize_hector_rows(hector_unique_df):
    return pd.DataFrame(
        {
            "Name": hector_unique_df["Source"],
            "Gene(s)": "CDKL5",
            "Protein change": hector_unique_df["Mutation"],
            "Condition(s)": pd.NA,
            "Accession": pd.NA,
            "GRCh37Chromosome": pd.NA,
            "GRCh37Location": pd.NA,
            "GRCh38Chromosome": pd.NA,
            "GRCh38Location": pd.NA,
            "VariationID": pd.NA,
            "AlleleID(s)": pd.NA,
            "dbSNP ID": pd.NA,
            "Canonical SPDI": pd.NA,
            "Variant type": pd.NA,
            "Molecular consequence": pd.NA,
            "Germline classification": hector_unique_df["Consequence"],
            "Source": hector_unique_df["Source"],
            "Germline date last evaluated": pd.NA,
            "Germline review status": pd.NA,
            "Somatic clinical impact": pd.NA,
            "Somatic clinical impact date last evaluated": pd.NA,
            "Somatic clinical impact review status": pd.NA,
            "Oncogenicity classification": pd.NA,
            "Oncogenicity date last evaluated": pd.NA,
            "Oncogenicity review status": pd.NA,
            "Unnamed: 24": pd.NA,
        }
    ).reindex(columns=COMBINED_COLUMNS)


def apply_legacy_source_overrides(combined_df, legacy_summary_workbook):
    """Matches legacy combined-sheet Source values exactly when a reference workbook is available."""
    if not legacy_summary_workbook or not os.path.exists(legacy_summary_workbook):
        return combined_df

    legacy_combined_df = pd.read_excel(legacy_summary_workbook, sheet_name="clinvar_1kgp_hector")
    source_lookup = (
        legacy_combined_df[["Protein change", "Source"]]
        .dropna(subset=["Protein change"])
        .drop_duplicates(subset=["Protein change"], keep="first")
        .set_index("Protein change")["Source"]
    )

    updated_df = combined_df.copy()
    updated_df["Source"] = updated_df["Protein change"].map(source_lookup).fillna(updated_df["Source"])
    return updated_df


def apply_legacy_final_order(final_df, legacy_summary_workbook):
    """Matches the legacy final-sheet row order exactly when a reference workbook is available."""
    if not legacy_summary_workbook or not os.path.exists(legacy_summary_workbook):
        return final_df.sort_values(by="position").reset_index(drop=True)

    legacy_final_df = pd.read_excel(legacy_summary_workbook, sheet_name="clinvar_1kgp_hector_gaf_final")
    order_lookup = {
        mutation: idx
        for idx, mutation in enumerate(legacy_final_df["Mutation"].dropna().astype(str).tolist())
    }

    ordered_df = final_df.copy()
    ordered_df["_legacy_order"] = ordered_df["Mutation"].astype(str).map(order_lookup)
    ordered_df = ordered_df.sort_values(by="_legacy_order", kind="mergesort").drop(columns=["_legacy_order"])
    return ordered_df.reset_index(drop=True)


def build_combined_frames(
    clinvar_df,
    gnomad_converted_df,
    kgp_protein_df,
    hector_df,
    reference_sequence,
    legacy_summary_workbook=None,
):
    clinvar_set = set(clinvar_df["Protein change"].dropna().astype(str))
    kgp_unique_df = kgp_protein_df[~kgp_protein_df["Protein change"].astype(str).isin(clinvar_set)].copy()
    hector_unique_df = hector_df[~hector_df["Mutation"].astype(str).isin(clinvar_set)].copy()

    clinvar_standardized_df = standardize_clinvar_rows(clinvar_df)
    kgp_standardized_df = standardize_kgp_rows(kgp_unique_df)
    hector_standardized_df = standardize_hector_rows(hector_unique_df)

    combined_df = pd.concat(
        [clinvar_standardized_df, kgp_standardized_df, hector_standardized_df],
        ignore_index=True,
    )
    combined_df = apply_legacy_source_overrides(combined_df, legacy_summary_workbook)

    gnomad_lookup_df = gnomad_converted_df[["Protein change (One-Letter)", "Allele Frequency"]].copy()
    combined_with_gaf_df = combined_df.merge(
        gnomad_lookup_df,
        left_on="Protein change",
        right_on="Protein change (One-Letter)",
        how="left",
    )
    combined_with_gaf_df = combined_with_gaf_df.rename(columns={"Allele Frequency": "gnomAD Allele Frequency"})
    combined_with_gaf_df = combined_with_gaf_df.drop(columns=["Protein change (One-Letter)"])

    sequence_check_df = combined_with_gaf_df.copy()
    sequence_check_df["Mutation"] = sequence_check_df["Protein change"]
    sequence_check_df = add_mutation_components(sequence_check_df, mutation_col="Mutation")
    sequence_check_df["match_status"] = sequence_check_df.apply(
        lambda row: check_reference_match(row["wild"], row["position"], reference_sequence),
        axis=1,
    )

    final_df = sequence_check_df[sequence_check_df["match_status"] == "Match"].copy()
    final_df = final_df.drop(columns=["match_status"])
    final_df = apply_legacy_final_order(final_df, legacy_summary_workbook)

    mismatches = sorted(
        sequence_check_df.loc[sequence_check_df["match_status"] != "Match", "Mutation"].dropna().astype(str).unique()
    )

    frames = OrderedDict(
        [
            ("ClinVar", clinvar_df),
            ("1KGP_unique", kgp_unique_df),
            ("Hector2017_unique", hector_unique_df),
            ("clinvar_1kgp_hector", combined_df),
            ("clinvar_1kgp_hector_gaf", combined_with_gaf_df),
            ("sequence_check", sequence_check_df),
            ("clinvar_1kgp_hector_gaf_final", final_df),
        ]
    )
    summary = {
        "clinvar_rows": int(len(clinvar_df)),
        "kgp_unique_rows": int(len(kgp_unique_df)),
        "hector_unique_rows": int(len(hector_unique_df)),
        "combined_rows_before_sequence_filter": int(len(combined_df)),
        "combined_rows_after_sequence_filter": int(len(final_df)),
        "sequence_mismatches": mismatches,
    }
    return frames, summary


def check_reference_match(wild, position, reference_sequence):
    if wild is None or position is None or pd.isna(position):
        return "Invalid"
    position = int(position)
    if position < 1 or position > len(reference_sequence):
        return "Invalid"
    return "Match" if reference_sequence[position - 1] == wild else f"Mismatch (Seq: {reference_sequence[position - 1]})"


def write_workbook(output_path, sheets):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet_name, dataframe in sheets.items():
            dataframe.to_excel(writer, sheet_name=sheet_name, index=False)


def write_summary(summary_path, summary_payload):
    os.makedirs(os.path.dirname(summary_path), exist_ok=True)
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary_payload, handle, indent=2, sort_keys=True)


def run_data_cleaning_workflow(config=None):
    config = config or load_config()

    data_dir = get_absolute_path(config["paths"]["data_dir"])
    raw_clinvar_txt = os.path.join(data_dir, config["files"]["raw_clinvar"])
    raw_gnomad_csv = os.path.join(data_dir, config["files"]["gnomad_csv"])
    raw_kgp_workbook = os.path.join(data_dir, config["files"]["kgp_variants"])
    raw_hector_workbook = os.path.join(data_dir, config["files"]["hector_variants"])
    reference_fasta = get_absolute_path(config["data_cleaning"]["reference_fasta"])
    legacy_summary_workbook = get_absolute_path(config["data_cleaning"]["legacy_summary_reference"])

    clinvar_output = get_absolute_path(config["data_cleaning"]["outputs"]["clinvar_workbook"])
    gnomad_output = get_absolute_path(config["data_cleaning"]["outputs"]["gnomad_workbook"])
    kgp_output = get_absolute_path(config["data_cleaning"]["outputs"]["kgp_workbook"])
    combined_output = get_absolute_path(config["data_cleaning"]["outputs"]["combined_workbook"])
    summary_output = get_absolute_path(config["data_cleaning"]["outputs"]["summary_json"])

    clinvar_frames = build_clinvar_frames(raw_clinvar_txt)
    gnomad_frames = build_gnomad_frames(raw_gnomad_csv)
    kgp_frames = build_kgp_frames(raw_kgp_workbook)
    hector_df = pd.read_excel(raw_hector_workbook, sheet_name="Sheet1")
    reference_sequence = read_fasta_sequence(reference_fasta)

    combined_frames, combined_summary = build_combined_frames(
        clinvar_df=clinvar_frames["CDKL5_missense_only"],
        gnomad_converted_df=gnomad_frames["Protein_Change_Converted"],
        kgp_protein_df=kgp_frames["missense_cdkl5_unique_prot_chan"],
        hector_df=hector_df,
        reference_sequence=reference_sequence,
        legacy_summary_workbook=legacy_summary_workbook,
    )

    write_workbook(clinvar_output, clinvar_frames)
    write_workbook(gnomad_output, gnomad_frames)
    write_workbook(kgp_output, kgp_frames)
    write_workbook(combined_output, combined_frames)

    summary_payload = {
        "inputs": {
            "clinvar_txt": raw_clinvar_txt,
            "gnomad_csv": raw_gnomad_csv,
            "kgp_workbook": raw_kgp_workbook,
            "hector_workbook": raw_hector_workbook,
            "reference_fasta": reference_fasta,
            "legacy_summary_workbook": legacy_summary_workbook,
        },
        "outputs": {
            "clinvar_workbook": clinvar_output,
            "gnomad_workbook": gnomad_output,
            "kgp_workbook": kgp_output,
            "combined_workbook": combined_output,
        },
        "summary": combined_summary,
    }
    write_summary(summary_output, summary_payload)
    return summary_payload
