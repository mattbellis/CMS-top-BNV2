# ttbar Semileptonic Analysis — CMS Run II/III

A modular, scalable analysis framework for reconstructing top quark pairs in the
semileptonic final state, built on [pocket-coffea](https://github.com/PocketCoffea/PocketCoffea).

---

## Repository Layout

```
ttbar-semilep/
├── README.md
├── environment.yml              # Conda environment spec
│
├── datasets/
│   ├── mc_samples.yaml          # MC process definitions + cross sections
│   ├── data_samples.yaml        # Collision data definitions
│   └── samples_metadata.yaml    # Human-readable labels, colors, groupings
│
├── config/
│   ├── analysis_config.py       # Main pocket-coffea Configurator
│   └── parameters/
│       ├── cuts.py              # Selection functions (objects + events)
│       └── weights.py           # Custom weight definitions
│
├── workflows/
│   ├── ttbar_semilep_processor.py   # Main processor (inherits BaseProcessorABC)
│   └── reconstruction.py            # Top quark reco helpers
│
├── scripts/
│   ├── dataset_diagnostics.py   # Produce CSV/parquet of dataset metadata
│   ├── run_analysis.py          # Unified runner: local / dask / condor
│   └── make_plots.py            # Quick diagnostic plots from outputs
│
├── corrections/                 # SF JSONs, JEC tarballs, pileup weights
│   └── README.md
│
├── notebooks/
│   └── explore_output.ipynb     # Interactive exploration of outputs
│
└── outputs/                     # gitignored; .coffea / .parquet outputs go here
```

---

## Setup

### 1. Create the environment

```bash
conda env create -f environment.yml
conda activate ttbar-semilep
```

Or on LPC/lxplus where you want the CMS software stack:

```bash
# On LPC
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv   # inside a CMSSW release, or use a standalone coffea env

pip install pocket-coffea --upgrade
```

### 2. Authenticate for xrootd (needed for remote file access)

```bash
voms-proxy-init --voms cms --valid 192:00
```

### 3. Build / update dataset file lists

pocket-coffea uses a JSON file listing all the ROOT files per dataset.
Build it from DAS with the provided script:

```bash
python scripts/dataset_diagnostics.py --build-filelists \
    --mc-yaml datasets/mc_samples.yaml \
    --data-yaml datasets/data_samples.yaml \
    --output-dir datasets/filelists/
```

This will also write `datasets/dataset_summary.csv` with event counts,
cross sections, and expected luminosities — useful for quick sanity checks.

---

## Workflow: from test to full scale

### Step 1 — Tiny local test (1 file per dataset, no scheduler)

```bash
python scripts/run_analysis.py \
    --config config/analysis_config.py \
    --limit-files 1 \
    --limit-chunks 2 \
    --executor iterative \
    --output outputs/test_run.coffea
```

`--limit-files 1 --limit-chunks 2` means you process at most 2 chunks of
~100k events from a single file per dataset. Runs in ~5 minutes on a laptop.

### Step 2 — Small local test with futures executor (use all local cores)

```bash
python scripts/run_analysis.py \
    --config config/analysis_config.py \
    --limit-files 3 \
    --executor futures \
    --workers 4 \
    --output outputs/local_test.coffea
```

### Step 3 — Dask on LPC / lxplus interactive node

```bash
# Start a dask cluster (see scripts/run_analysis.py --help for scheduler options)
python scripts/run_analysis.py \
    --config config/analysis_config.py \
    --executor dask \
    --dask-scheduler lpc \
    --output outputs/dask_test.coffea
```

### Step 4 — Full Run II on condor (LPC)

```bash
python scripts/run_analysis.py \
    --config config/analysis_config.py \
    --executor condor \
    --output outputs/RunII_full.coffea
```

### Step 5 — Inspect outputs on your laptop

Copy the `.coffea` output file to your laptop, then:

```bash
python scripts/make_plots.py --input outputs/RunII_full.coffea --output plots/
# or open the notebook:
jupyter lab notebooks/explore_output.ipynb
```

---

## Systematic variations

Variations are declared in `config/analysis_config.py` under `systematic_cfg`.
By default the framework runs:

- JEC/JER Up/Down
- b-tagging SF shape variations
- Lepton trigger/ID/ISO SF Up/Down
- Pileup weight Up/Down
- Top pT reweighting On/Off
- PDF + QCD scale variations (via LHE weights)

Each variation produces its own histogram set in the output without duplicating
the full event loop — pocket-coffea handles the branching internally.

---

## Adding a new MC sample

1. Add an entry to `datasets/mc_samples.yaml` with the DAS path and cross section.
2. Re-run `dataset_diagnostics.py --build-filelists`.
3. Add a row to `datasets/samples_metadata.yaml` for labels/groupings.
4. Done — the processor picks it up automatically.

---

## Key dependencies

| Package | Purpose |
|---|---|
| `pocket-coffea` | Analysis framework, dataset mgmt, systematics |
| `coffea` | NanoAOD reading, awkward arrays, histograms |
| `awkward` | Jagged array manipulation |
| `hist` | Histogram filling and serialization |
| `correctionlib` | CMS official SFs (b-tag, lepton, pileup) |
| `dask` / `distributed` | Scale-up scheduler |
| `pandas` | Diagnostics, CSV/parquet I/O |
| `matplotlib` / `mplhep` | Plotting with CMS style |
