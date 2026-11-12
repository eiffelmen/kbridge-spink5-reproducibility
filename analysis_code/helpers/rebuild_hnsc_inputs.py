#!/usr/bin/env python3
"""Build incremental and sensitivity inputs from the pinned public HNSCC analyses."""

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy

import numpy as np
import pandas as pd


TCGA_MODEL = "model_data_M2_age_clinical_stage_sex_smoking_site.tsv"
TCGA_PATIENT_JOIN = "analysis_primary_tumours.tsv"
GSE_MODEL = "primary_clinical_hpv_model_input.csv"
TCGA_METADATA_COLUMNS = (
    "sampleID",
    "SPINK5_log2",
    "time",
    "event",
    "HPV_p16_positive",
    "neoplasm_histologic_grade",
    "radiation_therapy",
    "targeted_molecular_therapy",
)
GSE_METADATA_COLUMNS = ("geo", "tumor_type", "os_event", "treatment")


def project_columns(source, columns, delimiter):
    """Preserve source cell text and the original reduced files' CRLF separators."""
    with source.open(encoding="utf-8", newline="") as stream:
        rows = csv.DictReader(stream, delimiter=delimiter)
        missing = set(columns).difference(rows.fieldnames or ())
        if missing:
            raise ValueError(f"Missing columns in {source}: {', '.join(sorted(missing))}")
        output = io.StringIO(newline="")
        writer = csv.writer(output, delimiter=delimiter)
        writer.writerow(columns)
        writer.writerows([row[column] for column in columns] for row in rows)
    return output.getvalue().encode("utf-8")


def join_tcga_ids(tcga_model, tcga_patient_join):
    """Rebuild and verify the existing design matrix before attaching sample IDs."""
    clinical = pd.read_csv(tcga_patient_join, sep="\t")
    required = (
        "time", "event", "SPINK5_per_SD", "age10", "clinical_stage_num",
        "male", "ever_smoked", "site_group",
    )
    complete = clinical.loc[clinical.time.gt(0)].dropna(subset=required)
    design = complete[["time", "event", "SPINK5_per_SD", "age10", "male", "ever_smoked"]].copy()
    for stage in (2, 3, 4):
        design[f"stage_{stage}_vs_1"] = complete.clinical_stage_num.eq(stage).astype(float)
    for site in ("larynx", "oropharynx", "hypopharynx"):
        design[f"site_{site}_vs_oral_cavity"] = complete.site_group.eq(site).astype(float)

    model = pd.read_csv(tcga_model, sep="\t")
    if list(design.columns) != list(model.columns):
        raise ValueError(f"Unexpected TCGA model columns in {tcga_model}")
    np.testing.assert_allclose(
        design.to_numpy(dtype=float), model.to_numpy(dtype=float), rtol=1e-10, atol=1e-10,
    )
    model.insert(0, "patient_id", complete.sampleID.to_numpy())
    return model.rename(columns={"SPINK5_per_SD": "SPINK5"})


def verify_cohort(data, cohort, size, deaths):
    if len(data) != size or data.event.sum() != deaths or not data.patient_id.is_unique:
        raise ValueError(f"{cohort} input does not match the published sample/event counts and unique IDs")


def extract_markers(signature_source, expression_source):
    populations = ("B lineage", "Fibroblasts")
    with signature_source.open(encoding="utf-8", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        genes = {row["HUGO symbols"] for row in rows if row["Cell population"] in populations}
    genes.add("SPINK5")
    with gzip.open(expression_source, mode="rt", encoding="utf-8") as stream:
        header = next(stream)
        rows = [line for line in stream if line.split("\t", 1)[0] in genes]
    observed = {line.split("\t", 1)[0] for line in rows}
    if observed != genes - {"IGKC"}:
        raise ValueError(f"Marker coverage does not match the published Xena input: {expression_source}")
    return (header + "".join(rows)).encode("utf-8")


def file_inputs(root, xena, puram, signatures):
    upstream_tcga = root / "public_hnsc" / "results"
    upstream_gse = root / "external_hnsc"
    incremental = root / "incremental_prognosis"
    sensitivity = root / "public_hnsc_sensitivity"
    m2 = upstream_tcga / TCGA_MODEL
    patient_join = upstream_tcga / TCGA_PATIENT_JOIN
    external = upstream_gse / GSE_MODEL

    sources = {
        "tcga_model": m2,
        "tcga_patient_join": patient_join,
        "external_model": external,
    }
    manifest = json.loads((incremental / "source_manifest.json").read_text(encoding="utf-8"))
    expected = {entry["key"]: entry["sha256"] for entry in manifest["sources"]}
    artifacts = {}
    for key, source in sources.items():
        contents = source.read_bytes()
        observed = hashlib.sha256(contents).hexdigest()
        if observed != expected[key]:
            raise ValueError(f"SHA-256 mismatch for upstream file {source}: {observed}")
        artifacts[incremental / "inputs" / source.name] = contents

    tcga = join_tcga_ids(m2, patient_join)
    external_model = pd.read_csv(external).rename(
        columns={"geo": "patient_id", "os_native": "time", "death": "event", "SPINK5_z": "SPINK5"}
    )
    verify_cohort(tcga, "TCGA-HNSC", 491, 210)
    verify_cohort(external_model, "GSE65858", 233, 73)
    tcga_verified = tcga.to_csv(sep="\t", index=False).encode("utf-8")
    gse_verified = external_model.to_csv(sep="\t", index=False).encode("utf-8")
    artifacts[incremental / "inputs" / "TCGA_verified_analysis_input.tsv"] = tcga_verified
    artifacts[incremental / "inputs" / "GSE65858_verified_analysis_input.tsv"] = gse_verified

    si = sensitivity / "inputs"
    artifacts[si / "tcga_main.tsv"] = tcga_verified
    artifacts[si / "gse65858_main.csv"] = artifacts[incremental / "inputs" / external.name]
    artifacts[si / "tcga_clinical.tsv"] = project_columns(patient_join, TCGA_METADATA_COLUMNS, "\t")
    artifacts[si / "gse65858_metadata.csv"] = project_columns(
        upstream_gse / "sample_metadata.csv", GSE_METADATA_COLUMNS, ","
    )
    artifacts[si / "mcp_counter_genes.tsv"] = signatures.read_bytes()
    artifacts[si / "tcga_marker_expression.tsv"] = extract_markers(signatures, xena)

    # Reuse the deposited-matrix parser without running the downstream analysis.
    extract_puram = runpy.run_path(str(sensitivity / "run_cellular_source.py"))["extract"]
    artifacts[si / "puram_gene_extract.tsv"] = extract_puram(puram).to_csv(sep="\t").encode("utf-8")
    return artifacts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="Release root containing public_hnsc, external_hnsc, and analysis directories.")
    parser.add_argument("--xena-expression", type=Path, help="Override the default fetched HiSeqV2.gz.")
    parser.add_argument("--puram-matrix", type=Path, help="Override the default fetched GSE103322 matrix.")
    parser.add_argument("--mcp-signatures", type=Path, help="Override the default fetched signature table.")
    args = parser.parse_args()
    root = args.root.resolve()
    sensitivity = root / "public_hnsc_sensitivity"
    artifacts = file_inputs(
        root,
        args.xena_expression or sensitivity / "source_inputs" / "HiSeqV2.gz",
        args.puram_matrix or sensitivity / "source_inputs" / "GSE103322_HNSCC_all_data.txt.gz",
        args.mcp_signatures or sensitivity / "inputs" / "mcp_counter_genes.tsv",
    )
    # Check all existing files before writing any, so a stale mix cannot pass silently.
    for target, contents in artifacts.items():
        if target.exists() and target.read_bytes() != contents:
            raise ValueError(f"Existing file differs from reconstructed input: {target}")
    for target, contents in artifacts.items():
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(contents)
        print(f"{hashlib.sha256(contents).hexdigest()}  {target.relative_to(root)}")


if __name__ == "__main__":
    main()
