#!/usr/bin/env python3
"""
dataset_diagnostics.py — Build file lists and produce a CSV summary of all
datasets (event counts, cross sections, expected yields at target luminosity).

Usage
-----
# Build file lists from DAS (requires valid VOMS proxy):
python scripts/dataset_diagnostics.py --build-filelists \\
    --mc-yaml datasets/mc_samples.yaml \\
    --data-yaml datasets/data_samples.yaml \\
    --output-dir datasets/filelists/

# Print summary table (no DAS query, uses cached filelists):
python scripts/dataset_diagnostics.py --summary \\
    --output-csv datasets/dataset_summary.csv

# Both at once:
python scripts/dataset_diagnostics.py --build-filelists --summary

Output
------
datasets/filelists/mc_samples_UL18.json    ← pocket-coffea filelist
datasets/filelists/data_samples_UL18.json
datasets/dataset_summary.csv               ← one row per sample
datasets/dataset_summary.parquet           ← same, for pandas

The CSV has columns:
    sample, process_group, era, is_data, n_files, n_events,
    cross_section_pb, int_lumi_per_pb, expected_yield,
    generator, das_path
"""

import os
import sys
import json
import subprocess
import argparse
from pathlib import Path

import yaml
import pandas as pd
from rich.console import Console
from rich.table import Table

console = Console()

# Default luminosity for yield calculation (UL18 full dataset)
DEFAULT_LUMI = 59722.768  # pb^-1


# ─────────────────────────────────────────────────────────────────────────────
# DAS query helpers
# ─────────────────────────────────────────────────────────────────────────────

def das_query(dataset_path: str, query_type: str = "file") -> list[str]:
    """
    Run a DAS query and return the results as a list of strings.
    Requires `dasgoclient` to be on PATH (available inside CMSSW).
    """
    query = f"{query_type} dataset={dataset_path}"
    cmd = ["dasgoclient", "-query", query]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True,
                                timeout=120)
        return [line.strip() for line in result.stdout.splitlines() if line.strip()]
    except FileNotFoundError:
        console.print("[red]dasgoclient not found. Are you inside a CMSSW environment?[/red]")
        raise
    except subprocess.CalledProcessError as e:
        console.print(f"[red]DAS query failed for {dataset_path}:\n{e.stderr}[/red]")
        return []


def get_file_list(dataset_path: str, xrootd_prefix: str = "root://cmsxrootd.fnal.gov/") -> list[str]:
    """Return a list of xrootd file paths for a DAS dataset path."""
    files = das_query(dataset_path, query_type="file")
    return [xrootd_prefix + f for f in files]


def get_event_count(dataset_path: str) -> int:
    """Return the total number of events in a DAS dataset."""
    result = das_query(dataset_path, query_type="summary")
    for line in result:
        try:
            data = json.loads(line)
            return int(data.get("nevents", -1))
        except json.JSONDecodeError:
            continue
    return -1


# ─────────────────────────────────────────────────────────────────────────────
# Build pocket-coffea filelists
# ─────────────────────────────────────────────────────────────────────────────

def build_filelist_json(samples: dict, output_path: Path, is_data: bool = False,
                        xrootd_prefix: str = "root://cmsxrootd.fnal.gov/"):
    """
    Build a pocket-coffea filelist JSON from a samples dict.

    pocket-coffea expects:
    {
        "sample_name": {
            "files": ["root://...", ...],
            "metadata": {...}
        },
        ...
    }
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out = {}

    for name, info in samples.items():
        console.print(f"  Querying DAS for [cyan]{name}[/cyan]...")
        files = get_file_list(info["das_path"], xrootd_prefix)
        n_events = get_event_count(info["das_path"])
        info["n_effective_events"] = n_events

        metadata = {
            "das_path":   info["das_path"],
            "era":        info.get("era", "unknown"),
            "is_data":    is_data,
            "n_events":   n_events,
        }
        if not is_data:
            metadata["cross_section"] = info["cross_section"]
            metadata["process_group"] = info.get("process_group", "unknown")
            metadata["generator"]     = info.get("generator", "unknown")
        else:
            metadata["int_lumi_per_pb"] = info.get("integrated_luminosity", -1)
            metadata["golden_json"]     = info.get("golden_json", "")

        out[name] = {"files": files, "metadata": metadata}
        console.print(f"    → {len(files)} files, {n_events:,} events")

    with open(output_path, "w") as f:
        json.dump(out, f, indent=2)

    console.print(f"[green]Wrote {output_path}[/green]")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Summary table
# ─────────────────────────────────────────────────────────────────────────────

def build_summary_df(mc_yaml_path: str, data_yaml_path: str,
                     filelist_dir: Path, lumi: float = DEFAULT_LUMI) -> pd.DataFrame:
    """
    Build a DataFrame summarising all samples.
    Reads n_events from the cached filelist JSONs if available,
    otherwise from the YAML (which may be -1 if not yet queried).
    """
    rows = []

    def load_filelist_cache(json_path: Path) -> dict:
        if json_path.exists():
            with open(json_path) as f:
                return json.load(f)
        return {}

    # Load cached file lists
    mc_cache   = load_filelist_cache(filelist_dir / "mc_samples_UL18.json")
    data_cache = load_filelist_cache(filelist_dir / "data_samples_UL18.json")

    # MC
    with open(mc_yaml_path) as f:
        mc_samples = yaml.safe_load(f)

    for name, info in mc_samples.items():
        cached = mc_cache.get(name, {})
        meta   = cached.get("metadata", {})
        n_events = meta.get("n_events", info.get("n_effective_events", -1))
        n_files  = len(cached.get("files", []))
        xsec     = info["cross_section"]
        expected_yield = xsec * lumi if n_events > 0 else -1

        rows.append({
            "sample":         name,
            "process_group":  info.get("process_group", ""),
            "era":            info.get("era", ""),
            "is_data":        False,
            "n_files":        n_files,
            "n_events":       n_events,
            "cross_section_pb": xsec,
            "int_lumi_per_pb":  None,
            "expected_yield": expected_yield,
            "generator":      info.get("generator", ""),
            "das_path":       info["das_path"],
        })

    # Data
    with open(data_yaml_path) as f:
        data_samples = yaml.safe_load(f)

    for name, info in data_samples.items():
        cached   = data_cache.get(name, {})
        meta     = cached.get("metadata", {})
        n_events = meta.get("n_events", -1)
        n_files  = len(cached.get("files", []))

        rows.append({
            "sample":           name,
            "process_group":    "data",
            "era":              info.get("era", ""),
            "is_data":          True,
            "n_files":          n_files,
            "n_events":         n_events,
            "cross_section_pb": None,
            "int_lumi_per_pb":  info.get("integrated_luminosity", -1),
            "expected_yield":   None,
            "generator":        "collision",
            "das_path":         info["das_path"],
        })

    return pd.DataFrame(rows)


def print_summary_table(df: pd.DataFrame):
    """Print a rich table to the terminal."""
    table = Table(title=f"Dataset Summary  (lumi = {DEFAULT_LUMI:.1f} pb⁻¹)",
                  show_lines=True)
    table.add_column("Sample",         style="cyan",  no_wrap=False, max_width=45)
    table.add_column("Group",          style="white")
    table.add_column("Era",            style="white")
    table.add_column("Files",          justify="right")
    table.add_column("Events",         justify="right")
    table.add_column("σ (pb)",         justify="right")
    table.add_column("Expected yield", justify="right")

    for _, row in df.iterrows():
        xsec = f"{row['cross_section_pb']:.3g}" if row["cross_section_pb"] else "—"
        yld  = f"{row['expected_yield']:.3g}"   if row["expected_yield"]   else "—"
        nev  = f"{row['n_events']:,}" if row["n_events"] > 0 else "unknown"
        nf   = str(row["n_files"]) if row["n_files"] > 0 else "—"

        table.add_row(
            row["sample"], row["process_group"], row["era"],
            nf, nev, xsec, yld,
        )

    console.print(table)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build pocket-coffea filelists and/or print dataset summary."
    )
    parser.add_argument("--mc-yaml",   default="datasets/mc_samples.yaml")
    parser.add_argument("--data-yaml", default="datasets/data_samples.yaml")
    parser.add_argument("--output-dir", default="datasets/filelists/")
    parser.add_argument("--output-csv", default="datasets/dataset_summary.csv")
    parser.add_argument("--lumi", type=float, default=DEFAULT_LUMI,
                        help="Target luminosity in pb^-1 for yield calculations")
    parser.add_argument("--xrootd-prefix", default="root://cmsxrootd.fnal.gov/",
                        help="xrootd redirector prefix (use cmsxrootd.cern.ch at CERN)")
    parser.add_argument("--build-filelists", action="store_true",
                        help="Query DAS and build filelist JSONs")
    parser.add_argument("--summary", action="store_true",
                        help="Print summary table and write CSV/parquet")
    return parser.parse_args()


def main():
    args = parse_args()
    filelist_dir = Path(args.output_dir)

    if not args.build_filelists and not args.summary:
        console.print("[yellow]Nothing to do. Pass --build-filelists and/or --summary.[/yellow]")
        sys.exit(0)

    if args.build_filelists:
        console.print("\n[bold]── Building MC filelists ──[/bold]")
        with open(args.mc_yaml) as f:
            mc_samples = yaml.safe_load(f)
        build_filelist_json(mc_samples, filelist_dir / "mc_samples_UL18.json",
                            is_data=False, xrootd_prefix=args.xrootd_prefix)

        console.print("\n[bold]── Building data filelists ──[/bold]")
        with open(args.data_yaml) as f:
            data_samples = yaml.safe_load(f)
        build_filelist_json(data_samples, filelist_dir / "data_samples_UL18.json",
                            is_data=True, xrootd_prefix=args.xrootd_prefix)

    if args.summary:
        console.print("\n[bold]── Dataset Summary ──[/bold]")
        df = build_summary_df(args.mc_yaml, args.data_yaml, filelist_dir, args.lumi)
        print_summary_table(df)

        # Write CSV and parquet
        df.to_csv(args.output_csv, index=False)
        parquet_path = args.output_csv.replace(".csv", ".parquet")
        df.to_parquet(parquet_path, index=False)
        console.print(f"\nWrote [green]{args.output_csv}[/green] and [green]{parquet_path}[/green]")

        # Quick stats
        mc_rows   = df[~df.is_data]
        data_rows = df[df.is_data]
        console.print(
            f"\nTotal MC samples:   {len(mc_rows)}"
            f"\nTotal data samples: {len(data_rows)}"
            f"\nTotal data lumi:    {data_rows['int_lumi_per_pb'].sum():.1f} pb⁻¹"
        )


if __name__ == "__main__":
    main()
