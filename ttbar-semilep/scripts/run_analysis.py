#!/usr/bin/env python3
"""
run_analysis.py — Thin wrapper around the `pocket-coffea run` CLI.
"""
import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--executor", default="iterative")
    parser.add_argument("--workers", "-s", type=int, default=4)
    parser.add_argument("--limit-files",  "-lf", type=int, default=None)
    parser.add_argument("--limit-chunks", "-lc", type=int, default=None)
    parser.add_argument("--chunksize", type=int, default=None)
    parser.add_argument("--process-separately", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    Path(args.output).mkdir(parents=True, exist_ok=True)

    cmd = ["pocket-coffea", "run",
           "--cfg", args.config,
           "-o",    args.output,
           "--executor", args.executor]

    if args.executor == "futures":
        cmd += ["-s", str(args.workers)]
    if args.limit_files  is not None: cmd += ["-lf", str(args.limit_files)]
    if args.limit_chunks is not None: cmd += ["-lc", str(args.limit_chunks)]
    if args.chunksize    is not None: cmd += ["--chunksize", str(args.chunksize)]
    if args.process_separately:       cmd += ["--process-separately"]

    print("Running:", " ".join(cmd))
    sys.exit(subprocess.run(cmd).returncode)


if __name__ == "__main__":
    main()
