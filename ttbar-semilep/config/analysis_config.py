"""
analysis_config.py — pocket-coffea Configurator for the ttbar semileptonic analysis.

Run this file directly to validate the configuration:
    python config/analysis_config.py

Pass it to run_analysis.py via --config.
"""

from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.parameters.cuts.preselection_cuts import (
    passthrough,
    goldenJson,
)
from pocket_coffea.parameters.histograms import HistConf, Axis

from workflows.ttbar_semilep_processor import TTbarSemilepProcessor
from config.parameters.cuts import (
    semilep_preselection,
    trigger_selection_data,
    met_filter_selection,
    high_btag_region,
    low_btag_control_region,
)
from config.parameters.weights import TopPtReweighting, LHEScaleVariation

import os

# ─────────────────────────────────────────────────────────────────────────────
# Dataset JSON paths (built by dataset_diagnostics.py)
# ─────────────────────────────────────────────────────────────────────────────

FILELIST_DIR = os.path.join(os.path.dirname(__file__), "../datasets/filelists")

cfg = Configurator(
    # ── Workflow ──────────────────────────────────────────────────────────────
    workflow=TTbarSemilepProcessor,

    # ── Datasets ─────────────────────────────────────────────────────────────
    # pocket-coffea expects a dict mapping sample name → filelist JSON.
    # These JSONs are produced by dataset_diagnostics.py.
    datasets={
        "jsons": [
            os.path.join(FILELIST_DIR, "mc_samples_UL18.json"),
            os.path.join(FILELIST_DIR, "data_samples_UL18.json"),
        ],
        "filter": {
            "samples": None,   # None = use all; or list of names to restrict
        },
    },

    # ── Parameters passed to processor and cuts ───────────────────────────────
    parameters={},

    # ── Weights ───────────────────────────────────────────────────────────────
    weights_cfg={
        "common": {
            "inclusive": [
                "genWeight",
                "lumi",
                "XS",
                "pileup",
                "sf_mu_id",
                "sf_mu_iso",
                "sf_mu_trigger",
                "sf_btag_shape",   # shape-based b-tag SF (DeepJet)
            ],
            "bysample": {},   # add sample-specific weights here if needed
        },
        "custom_weights": [
            TopPtReweighting,
            LHEScaleVariation,
        ],
    },

    # ── Systematic variations ─────────────────────────────────────────────────
    # Shape systematics (produce separate Up/Down histograms)
    systematic_cfg={
        "shape": [
            "JES",     # Jet Energy Scale (full set of JEC sources)
            "JER",     # Jet Energy Resolution
            "top_pt_reweighting",
            "qcd_scale",
            "sf_mu_id",
            "sf_mu_iso",
            "sf_mu_trigger",
            "sf_btag_shape",
            "pileup",
        ],
    },

    # ── Selection / categories ────────────────────────────────────────────────
    # pocket-coffea applies cuts in order; each defines a "category" in the
    # output histograms.  You get cutflow counts automatically.
    cut_cfg={
        "preselections": [
            # MET filters (data + MC)
            Cut(
                name="met_filters",
                function=met_filter_selection,
                samples="all",
            ),
            # Trigger (data only; MC is trigger-unbiased)
            Cut(
                name="trigger",
                function=trigger_selection_data,
                samples=["SingleMuon_UL18A", "SingleMuon_UL18B",
                         "SingleMuon_UL18C", "SingleMuon_UL18D"],
            ),
            # Main semileptonic selection
            Cut(
                name="semilep_preselection",
                function=semilep_preselection,
                samples="all",
            ),
        ],
        "categories": {
            # Signal region
            "SR": [
                Cut(name="high_btag", function=high_btag_region, samples="all"),
            ],
            # Control region for QCD / W+jets
            "CR_1btag": [
                Cut(name="low_btag", function=low_btag_control_region, samples="all"),
            ],
        },
    },

    # ── Output variables (saved as flat columns) ──────────────────────────────
    # Saved to parquet; one row per event passing all cuts.
    columns={
        "common": {
            "inclusive": [
                ColOut("events", ["run", "luminosityBlock", "event"]),
                ColOut("GoodMuons", ["pt", "eta", "phi", "mass"], flatten=True),
                ColOut("GoodJets",  ["pt", "eta", "phi", "mass",
                                     "btagDeepFlavB"], flatten=True),
                ColOut("MET", ["pt", "phi"]),
            ],
        },
    },

    # ── Histograms ────────────────────────────────────────────────────────────
    variables={
        # Muon kinematics
        "muon_pt": HistConf(
            [Axis(coll="GoodMuons", field="pt", bins=50, start=0, stop=300,
                  label=r"Muon $p_T$ [GeV]", pos=0)]
        ),
        "muon_eta": HistConf(
            [Axis(coll="GoodMuons", field="eta", bins=50, start=-2.5, stop=2.5,
                  label=r"Muon $\eta$", pos=0)]
        ),
        # Jet kinematics
        "jet_pt": HistConf(
            [Axis(coll="GoodJets", field="pt", bins=50, start=0, stop=500,
                  label=r"Jet $p_T$ [GeV]")]
        ),
        "jet_eta": HistConf(
            [Axis(coll="GoodJets", field="eta", bins=50, start=-2.5, stop=2.5,
                  label=r"Jet $\eta$")]
        ),
        "jet_multiplicity": HistConf(
            [Axis(coll="events", field="nGoodJets", bins=12, start=0, stop=12,
                  label="Number of jets")]
        ),
        "bjet_multiplicity": HistConf(
            [Axis(coll="events", field="nBJets", bins=8, start=0, stop=8,
                  label="Number of b-jets")]
        ),
        # MET
        "met_pt": HistConf(
            [Axis(coll="MET", field="pt", bins=50, start=0, stop=400,
                  label=r"$p_T^{miss}$ [GeV]")]
        ),
        # Reconstructed top mass (filled in processor after reco)
        "reco_top_mass": HistConf(
            [Axis(coll="events", field="RecoTopMass", bins=60, start=0, stop=400,
                  label=r"Reconstructed top mass [GeV]")]
        ),
        "reco_W_mass": HistConf(
            [Axis(coll="events", field="RecoWMass", bins=60, start=0, stop=200,
                  label=r"Reconstructed $W$ mass [GeV]")]
        ),
        # HT
        "HT": HistConf(
            [Axis(coll="events", field="HT", bins=60, start=0, stop=2000,
                  label=r"$H_T$ [GeV]")]
        ),
    },

    # ── General settings ──────────────────────────────────────────────────────
    workflow_options={
        "chunksize": 100_000,    # events per dask task; tune for memory
        "maxchunks": None,       # None = process everything; set int for tests
    },
)


if __name__ == "__main__":
    # Quick validation — prints sample list and cut summary
    print(cfg)
