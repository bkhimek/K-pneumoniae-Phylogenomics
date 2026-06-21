#!/usr/bin/env python3
"""
Merge QUAST, AMRFinder, and Kleborate outputs into a single summary TSV.
Usage: python merge_summary.py --results <results_dir> --out summary.tsv
"""

import argparse
import os
import pandas as pd
from pathlib import Path


def load_quast(results_dir: Path) -> pd.DataFrame:
    frames = []
    for f in sorted((results_dir / "quast").glob("*_transposed_report.tsv")):
        df = pd.read_csv(f, sep="\t")          # already transposed — no .T needed
        sample_id = f.name.replace("_transposed_report.tsv", "")
        df.insert(0, "sample_id", sample_id)
        df = df.drop(columns=["Assembly"], errors="ignore")
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No QUAST reports found in {results_dir}/quast/")
    merged = pd.concat(frames, ignore_index=True)
    keep = ["sample_id", "# contigs", "Total length", "N50", "L50", "# N's per 100 kbp"]
    keep = [c for c in keep if c in merged.columns]
    return merged[keep].rename(columns={
        "# contigs": "contigs",
        "Total length": "total_length_bp",
        "# N's per 100 kbp": "Ns_per_100kbp"
    })

def load_kleborate(results_dir: Path) -> pd.DataFrame:
    frames = []
    for d in sorted((results_dir / "kleborate").iterdir()):
        if not d.is_dir():
            continue
        tsvs = list(d.glob("klebsiella_pneumo_complex_output.txt"))
        if not tsvs:
            continue
        df = pd.read_csv(tsvs[0], sep="\t")
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No Kleborate reports found in {results_dir}/kleborate/")
    merged = pd.concat(frames, ignore_index=True)

    # Rename v3 hierarchical columns to short names
    rename = {
        "strain":                                           "sample_id",
        "enterobacterales__species__species":               "species",
        "klebsiella_pneumo_complex__mlst__ST":              "ST",
        "klebsiella__ybst__Yersiniabactin":                 "Ybt",
        "klebsiella__cbst__Colibactin":                     "Clb",
        "klebsiella__abst__Aerobactin":                     "Iuc",
        "klebsiella__smst__Salmochelin":                    "Iro",
        "klebsiella__rmst__RmpADC":                         "RmpADC",
        "klebsiella__rmst__rmpA2":                          "rmpA2",
        "klebsiella_pneumo_complex__amr__virulence_score":  "virulence_score",
        "klebsiella_pneumo_complex__amr__resistance_score": "resistance_score",
        "klebsiella_pneumo_complex__amr__Bla_carb":         "Bla_Carb",
        "klebsiella_pneumo_complex__amr__Bla_ESBL":         "Bla_ESBL",
        "klebsiella_pneumo_complex__amr__Bla_inhR":         "Bla_inhR",
    }
    # Only rename columns that exist
    rename = {k: v for k, v in rename.items() if k in merged.columns}
    merged = merged.rename(columns=rename)
    merged["sample_id"] = merged["sample_id"].str.replace(r'\..*$', '', regex=True)
    keep = [v for v in rename.values() if v in merged.columns]
    return merged[keep]

def load_amrfinder(results_dir: Path) -> pd.DataFrame:
    frames = []
    for f in sorted((results_dir / "amrfinder").glob("*_amrfinder.tsv")):
        df = pd.read_csv(f, sep="\t")
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No AMRFinder reports found in {results_dir}/amrfinder/")
#    all_hits = pd.concat(frames, ignore_index=True)

    # Pivot: one row per sample, columns = gene names, value = element subtype
    if "Name" not in all_hits.columns:
        return pd.DataFrame()

    pivot = (
        all_hits.groupby(["Name", "Element symbol"])["Subtype"]
        .first()
        .unstack("Element symbol")
        .reset_index()
        .rename(columns={"Name": "sample_id"})
    )
    pivot.columns.name = None
    return pivot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, help="Path to results/ directory")
    parser.add_argument("--out", default="summary.tsv", help="Output file path")
    args = parser.parse_args()

    results_dir = Path(args.results)

    print("Loading QUAST …")
    quast_df = load_quast(results_dir)

    print("Loading Kleborate …")
    kleb_df = load_kleborate(results_dir)

    print("Loading AMRFinder …")
    amr_df = load_amrfinder(results_dir)

    print("Merging …")
    summary = quast_df.merge(kleb_df, on="sample_id", how="left")
    if not amr_df.empty:
        summary = summary.merge(amr_df, on="sample_id", how="left")

    summary.to_csv(args.out, sep="\t", index=False)
    print(f"Done → {args.out}  ({len(summary)} samples, {len(summary.columns)} columns)")


if __name__ == "__main__":
    main()
