#!/usr/bin/env python3
"""
make_plots.py — Read a .coffea output file and produce stacked plots.

Usage
-----
python scripts/make_plots.py \\
    --input outputs/RunII_full.coffea \\
    --metadata datasets/samples_metadata.yaml \\
    --lumi 59722.768 \\
    --output plots/

Produces one PDF per histogram variable in the config, with:
  - Stacked MC histogram (with stat uncertainty band)
  - Data points with error bars
  - Data/MC ratio panel
  - CMS preliminary label
"""

import os
import argparse
from pathlib import Path

import yaml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mplhep as hep
import coffea.util
import hist

plt.style.use(hep.style.CMS)


# ─────────────────────────────────────────────────────────────────────────────

def load_output(path: str) -> dict:
    return coffea.util.load(path)


def load_metadata(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


# ─────────────────────────────────────────────────────────────────────────────

def stack_plot(ax_main, ax_ratio, var_name: str, output: dict,
               metadata: dict, lumi: float, category: str = "SR"):
    """Draw a single stacked MC + data comparison for `var_name`."""

    groups    = metadata["process_groups"]
    sig_name  = metadata.get("signal_sample", None)
    data_grp  = groups.get("data", {})

    # ── Collect MC histograms by process group ────────────────────────────────
    mc_hists   = {}   # group_name → summed hist
    mc_colors  = {}
    mc_labels  = {}
    mc_orders  = {}

    for gname, ginfo in groups.items():
        if gname == "data":
            continue
        h_sum = None
        for sname in ginfo.get("samples", []):
            key = (sname, category, "nominal")
            try:
                h = output["histograms"][var_name][key]
            except KeyError:
                continue
            h_scaled = h * lumi  # multiply by xs*lumi weight already in hist?
            h_sum = h_scaled if h_sum is None else h_sum + h_scaled

        if h_sum is not None:
            mc_hists[gname]  = h_sum
            mc_colors[gname] = ginfo.get("color", "#999999")
            mc_labels[gname] = ginfo.get("label", gname)
            mc_orders[gname] = ginfo.get("stack_order", 99)

    # Sort by stack order
    sorted_groups = sorted(mc_hists, key=lambda g: mc_orders[g], reverse=True)

    # ── Collect data histogram ────────────────────────────────────────────────
    h_data = None
    for sname in data_grp.get("samples", []):
        key = (sname, category, "nominal")
        try:
            h = output["histograms"][var_name][key]
            h_data = h if h_data is None else h_data + h
        except KeyError:
            continue

    if not mc_hists and h_data is None:
        return  # nothing to plot

    # ── Draw ──────────────────────────────────────────────────────────────────
    bottoms = None
    patches = []
    mc_total = None

    for gname in sorted_groups:
        h = mc_hists[gname]
        edges  = h.axes[0].edges
        values = h.values()

        if bottoms is None:
            bottoms = np.zeros_like(values)

        ax_main.bar(
            edges[:-1], values, width=np.diff(edges),
            bottom=bottoms, align="edge",
            color=mc_colors[gname], label=mc_labels[gname],
            linewidth=0.5, edgecolor="black",
        )
        bottoms += values
        mc_total = bottoms.copy()

        patch = mpatches.Patch(color=mc_colors[gname], label=mc_labels[gname])
        patches.append(patch)

    # MC stat uncertainty band
    if mc_total is not None:
        mc_err = np.sqrt(mc_total)  # Poisson approx; use proper variances if available
        ax_main.fill_between(
            np.repeat(edges, 2)[1:-1],
            np.repeat(mc_total - mc_err, 2),
            np.repeat(mc_total + mc_err, 2),
            alpha=0.3, color="grey", label="MC stat. unc.",
            step="mid",
        )

    # Data points
    if h_data is not None:
        centers = 0.5 * (edges[:-1] + edges[1:])
        dvals   = h_data.values()
        derr    = np.sqrt(dvals)
        ax_main.errorbar(
            centers, dvals, yerr=derr,
            fmt="o", color="black", markersize=4, linewidth=1,
            label="Data",
        )

    # Ratio panel
    if h_data is not None and mc_total is not None:
        ratio = np.where(mc_total > 0, h_data.values() / mc_total, np.nan)
        rerr  = np.where(mc_total > 0,
                         np.sqrt(h_data.values()) / mc_total, np.nan)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax_ratio.errorbar(centers, ratio, yerr=rerr,
                          fmt="o", color="black", markersize=3, linewidth=0.8)
        ax_ratio.axhline(1.0, color="grey", linestyle="--", linewidth=0.8)
        mc_ratio_err = np.where(mc_total > 0, mc_err / mc_total, 0.0)
        ax_ratio.fill_between(
            np.repeat(edges, 2)[1:-1],
            np.repeat(1.0 - mc_ratio_err, 2),
            np.repeat(1.0 + mc_ratio_err, 2),
            alpha=0.3, color="grey", step="mid",
        )
        ax_ratio.set_ylim(0.5, 1.5)
        ax_ratio.set_ylabel("Data/MC")


# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",    required=True)
    parser.add_argument("--metadata", default="datasets/samples_metadata.yaml")
    parser.add_argument("--lumi",     type=float, default=59722.768)
    parser.add_argument("--output",   default="plots/")
    parser.add_argument("--category", default="SR",
                        help="Analysis category / region to plot")
    parser.add_argument("--format",   default="pdf",
                        choices=["pdf", "png", "svg"])
    args = parser.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    output   = load_output(args.input)
    metadata = load_metadata(args.metadata)

    histogram_names = list(output.get("histograms", {}).keys())
    print(f"Found {len(histogram_names)} histograms: {histogram_names}")

    for var_name in histogram_names:
        fig, (ax_main, ax_ratio) = plt.subplots(
            2, 1, figsize=(8, 7),
            gridspec_kw={"height_ratios": [3, 1]},
            sharex=True,
        )
        fig.subplots_adjust(hspace=0.05)

        stack_plot(ax_main, ax_ratio, var_name, output, metadata,
                   args.lumi, category=args.category)

        hep.cms.label(
            "Preliminary",
            data=True,
            lumi=args.lumi / 1000.0,   # fb^-1
            ax=ax_main,
        )
        ax_main.set_yscale("log")
        ax_main.legend(fontsize=8, ncol=2)
        ax_main.set_ylabel("Events / bin")
        ax_ratio.set_xlabel(var_name.replace("_", " "))

        out_path = out_dir / f"{var_name}.{args.format}"
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        print(f"  Wrote {out_path}")

    print(f"\nAll plots saved to {out_dir}/")


if __name__ == "__main__":
    main()
