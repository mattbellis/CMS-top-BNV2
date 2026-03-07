#!/usr/bin/env python3
"""
run_analysis.py — Unified runner for the ttbar semileptonic analysis.

Examples
--------
# Tiny test (1 file, 2 chunks, no parallelism):
python scripts/run_analysis.py \\
    --config config/analysis_config.py \\
    --limit-files 1 --limit-chunks 2 \\
    --executor iterative \\
    --output outputs/test.coffea

# Local parallel (uses all cores):
python scripts/run_analysis.py \\
    --config config/analysis_config.py \\
    --limit-files 3 \\
    --executor futures --workers 4 \\
    --output outputs/local.coffea

# Dask on LPC:
python scripts/run_analysis.py \\
    --config config/analysis_config.py \\
    --executor dask --dask-scheduler lpc \\
    --output outputs/dask.coffea

# Dask on lxplus:
python scripts/run_analysis.py \\
    --config config/analysis_config.py \\
    --executor dask --dask-scheduler lxplus \\
    --output outputs/lxplus.coffea

# Only run specific samples (useful for debugging one process):
python scripts/run_analysis.py \\
    --config config/analysis_config.py \\
    --samples TTToSemiLeptonic SingleMuon_UL18A \\
    --executor iterative --limit-files 1 \\
    --output outputs/debug.coffea

# Run without systematics (much faster, for cutflow checks):
python scripts/run_analysis.py \\
    --config config/analysis_config.py \\
    --no-systematics \\
    --executor futures \\
    --output outputs/nosys.coffea
"""

import os
import sys
import argparse
import importlib.util
from pathlib import Path
import time

from pocket_coffea.utils.runner import run_analysis
from pocket_coffea.utils.output import save_output
from rich.console import Console

console = Console()


def load_config(config_path: str):
    """Dynamically import the Configurator from the config module."""
    spec = importlib.util.spec_from_file_location("analysis_config", config_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "cfg"):
        raise AttributeError(f"{config_path} must define a `cfg` Configurator object.")
    return module.cfg


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the ttbar semileptonic analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Config
    parser.add_argument("--config", required=True,
                        help="Path to analysis_config.py")
    # Output
    parser.add_argument("--output", required=True,
                        help="Output .coffea file path")
    # Executor
    parser.add_argument("--executor",
                        choices=["iterative", "futures", "dask", "condor"],
                        default="iterative")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of local workers (futures executor)")
    parser.add_argument("--dask-scheduler",
                        choices=["lpc", "lxplus", "local", "slurm"],
                        default="local",
                        help="Dask cluster type (dask executor only)")
    # Limiting (for tests)
    parser.add_argument("--limit-files",  type=int, default=None,
                        help="Max files per sample (None = all)")
    parser.add_argument("--limit-chunks", type=int, default=None,
                        help="Max chunks per file (None = all)")
    parser.add_argument("--chunksize", type=int, default=100_000,
                        help="Events per chunk")
    # Sample filter
    parser.add_argument("--samples", nargs="*", default=None,
                        help="Restrict to these sample names")
    # Systematics
    parser.add_argument("--no-systematics", action="store_true",
                        help="Disable all systematic variations (nominal only)")
    return parser.parse_args()


def setup_dask_cluster(scheduler: str):
    """Return a dask Client configured for the requested cluster type."""
    from dask.distributed import Client

    if scheduler == "local":
        return Client()

    elif scheduler == "lpc":
        from lpcjobqueue import LPCCondorCluster
        cluster = LPCCondorCluster(
            cores=1, memory="4GB",
            ship_env=True,
            transfer_input_files=["config/", "workflows/", "corrections/"],
        )
        cluster.adapt(minimum=10, maximum=500)
        return Client(cluster)

    elif scheduler == "lxplus":
        from dask_jobqueue import HTCondorCluster
        cluster = HTCondorCluster(
            cores=1, memory="4GB", disk="2GB",
            job_extra_directives={
                "+JobFlavour": '"longnano"',
                "x509userproxy": os.environ.get("X509_USER_PROXY", ""),
            },
        )
        cluster.adapt(minimum=10, maximum=500)
        return Client(cluster)

    elif scheduler == "slurm":
        from dask_jobqueue import SLURMCluster
        cluster = SLURMCluster(
            cores=1, memory="4GB",
            walltime="04:00:00",
        )
        cluster.adapt(minimum=5, maximum=200)
        return Client(cluster)

    else:
        raise ValueError(f"Unknown scheduler: {scheduler}")


def main():
    args = parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # ── Load config ───────────────────────────────────────────────────────────
    console.print(f"\n[bold]Loading config from[/bold] {args.config}")
    cfg = load_config(args.config)

    # Apply overrides from CLI
    if args.limit_files:
        cfg.workflow_options["limit_files"]  = args.limit_files
    if args.limit_chunks:
        cfg.workflow_options["limit_chunks"] = args.limit_chunks
    if args.chunksize:
        cfg.workflow_options["chunksize"]    = args.chunksize
    if args.samples:
        cfg.datasets["filter"]["samples"] = args.samples
    if args.no_systematics:
        cfg.systematic_cfg = {"shape": []}

    console.print(cfg)

    # ── Run ───────────────────────────────────────────────────────────────────
    t0 = time.time()

    if args.executor == "iterative":
        console.print("\n[bold yellow]Running iteratively (single thread, good for debugging).[/bold yellow]")
        from coffea.processor import iterative_executor, Runner
        executor = iterative_executor()
        runner   = Runner(executor=executor, schema=cfg.schema,
                          chunksize=cfg.workflow_options["chunksize"],
                          maxchunks=cfg.workflow_options.get("limit_chunks"))
        output = run_analysis(cfg, runner)

    elif args.executor == "futures":
        console.print(f"\n[bold yellow]Running with {args.workers} local workers.[/bold yellow]")
        from coffea.processor import futures_executor, Runner
        executor = futures_executor(workers=args.workers)
        runner   = Runner(executor=executor, schema=cfg.schema,
                          chunksize=cfg.workflow_options["chunksize"],
                          maxchunks=cfg.workflow_options.get("limit_chunks"))
        output = run_analysis(cfg, runner)

    elif args.executor in ("dask", "condor"):
        console.print(f"\n[bold yellow]Setting up Dask cluster ({args.dask_scheduler})…[/bold yellow]")
        client = setup_dask_cluster(args.dask_scheduler)
        console.print(f"Dask dashboard: [link]{client.dashboard_link}[/link]")
        from coffea.processor import dask_executor, Runner
        executor = dask_executor(client=client)
        runner   = Runner(executor=executor, schema=cfg.schema,
                          chunksize=cfg.workflow_options["chunksize"],
                          maxchunks=cfg.workflow_options.get("limit_chunks"))
        output = run_analysis(cfg, runner)
        client.close()

    else:
        console.print(f"[red]Unknown executor: {args.executor}[/red]")
        sys.exit(1)

    # ── Save ──────────────────────────────────────────────────────────────────
    save_output(output, args.output)
    elapsed = time.time() - t0
    console.print(
        f"\n[green bold]Done in {elapsed:.1f}s.[/green bold]"
        f"  Output: [cyan]{args.output}[/cyan]"
    )


if __name__ == "__main__":
    main()
