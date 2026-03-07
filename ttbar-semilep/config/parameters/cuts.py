"""
cuts.py — Object preselection and event selection for ttbar semileptonic.

All functions follow the pocket-coffea convention:
    cut_function(events, params, year, sample, **kwargs) -> ak.Array(bool)

Object collections are attached to events in the processor's
apply_object_preselection method before cuts are called.
"""

import awkward as ak
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Object selection parameters
# These are referenced in analysis_config.py and can be overridden there.
# ─────────────────────────────────────────────────────────────────────────────

OBJECT_PARAMS = {
    "muon": {
        "pt_min": 26.0,       # GeV; above SingleMuon HLT threshold
        "eta_max": 2.4,
        "iso_max": 0.15,       # PFRelIso04_all (tight working point)
        "id": "tightId",
    },
    "electron": {
        "pt_min": 15.0,
        "eta_max": 2.4,
        "veto_id": "cutBased",  # cutBased >= 1 (veto WP) to reject extra leptons
    },
    "jet": {
        "pt_min": 30.0,
        "eta_max": 2.4,
        "jet_id": 6,            # tightLepVeto
        "btag_wp": 0.2783,      # DeepJet medium WP, UL18
    },
    "fatjet": {
        "pt_min": 200.0,
        "eta_max": 2.4,
        "msoftdrop_min": 50.0,
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Object preselection helpers
# These are called by the processor to build collections that persist
# across the event loop as events.GoodMuons, events.GoodJets, etc.
# ─────────────────────────────────────────────────────────────────────────────

def select_good_muons(events):
    """Tight muons for the signal lepton requirement."""
    mu = events.Muon
    p = OBJECT_PARAMS["muon"]
    mask = (
        (mu.pt > p["pt_min"])
        & (np.abs(mu.eta) < p["eta_max"])
        & (mu.tightId)
        & (mu.pfRelIso04_all < p["iso_max"])
    )
    return mu[mask]


def select_veto_muons(events):
    """Loose muons — events with extras are vetoed."""
    mu = events.Muon
    mask = (
        (mu.pt > 15.0)
        & (np.abs(mu.eta) < 2.4)
        & (mu.looseId)
        & (mu.pfRelIso04_all < 0.25)
    )
    return mu[mask]


def select_veto_electrons(events):
    """Veto-WP electrons — events with any are vetoed."""
    el = events.Electron
    mask = (
        (el.pt > 15.0)
        & (np.abs(el.eta) < 2.5)
        & (el.cutBased >= 1)   # veto WP
    )
    return el[mask]


def select_good_jets(events):
    """AK4 jets passing kinematic and ID requirements, cleaned vs good muons."""
    j = events.Jet
    p = OBJECT_PARAMS["jet"]

    # Basic kinematic + ID
    base = (
        (j.pt > p["pt_min"])
        & (np.abs(j.eta) < p["eta_max"])
        & (j.jetId == p["jet_id"])
    )

    # DR cleaning vs tight muons (already attached to events by the processor)
    if hasattr(events, "GoodMuons"):
        dr_ok = ak.all(
            j.metric_table(events.GoodMuons) > 0.4, axis=2
        )
        base = base & dr_ok

    return j[base]


def select_btagged_jets(good_jets, wp: float = None):
    """b-tagged subset of good jets using DeepJet score."""
    wp = wp or OBJECT_PARAMS["jet"]["btag_wp"]
    return good_jets[good_jets.btagDeepFlavB > wp]


# ─────────────────────────────────────────────────────────────────────────────
# Event-level selection functions
# Each returns a boolean array of length n_events.
# ─────────────────────────────────────────────────────────────────────────────

def semilep_preselection(events, params, year, sample, **kwargs):
    """
    Baseline semileptonic ttbar preselection:
      - Exactly 1 tight muon (signal lepton)
      - 0 veto electrons
      - Veto extra loose muons
      - >= 4 good AK4 jets
      - >= 2 b-tagged jets
      - MET > 20 GeV
    """
    n_good_mu   = ak.num(events.GoodMuons)
    n_veto_mu   = ak.num(events.VetoMuons)
    n_veto_el   = ak.num(events.VetoElectrons)
    n_good_jets = ak.num(events.GoodJets)
    n_bjets     = ak.num(events.BJets)
    met         = events.MET.pt

    return (
        (n_good_mu == 1)
        & (n_veto_mu == 1)      # tight is subset of veto; extra loose → reject
        & (n_veto_el == 0)
        & (n_good_jets >= 4)
        & (n_bjets >= 2)
        & (met > 20.0)
    )


def trigger_selection_data(events, params, year, sample, **kwargs):
    """
    Apply single-muon HLT trigger (data only).
    Triggers are era-dependent; this covers UL18 as default.
    """
    trigger_map = {
        "UL16APV": ["HLT_IsoMu24", "HLT_IsoTkMu24"],
        "UL16":    ["HLT_IsoMu24", "HLT_IsoTkMu24"],
        "UL17":    ["HLT_IsoMu27"],
        "UL18":    ["HLT_IsoMu24"],
    }
    triggers = trigger_map.get(year, ["HLT_IsoMu24"])
    fired = ak.zeros_like(events.run, dtype=bool)
    for trig in triggers:
        if hasattr(events.HLT, trig):
            fired = fired | getattr(events.HLT, trig)
    return fired


def met_filter_selection(events, params, year, sample, **kwargs):
    """Standard CMS MET filters for NanoAOD."""
    flags = events.Flag
    return (
        flags.goodVertices
        & flags.globalSuperTightHalo2016Filter
        & flags.HBHENoiseFilter
        & flags.HBHENoiseIsoFilter
        & flags.EcalDeadCellTriggerPrimitiveFilter
        & flags.BadPFMuonFilter
        & flags.BadPFMuonDzFilter
        & flags.eeBadScFilter
        & flags.ecalBadCalibFilter
    )


def high_btag_region(events, params, year, sample, **kwargs):
    """Signal region: >= 2 b-tags (already in preselection, but useful as
    an explicit named cut for the cutflow)."""
    return ak.num(events.BJets) >= 2


def low_btag_control_region(events, params, year, sample, **kwargs):
    """Control region: exactly 1 b-tag (for QCD / W+jets studies)."""
    return ak.num(events.BJets) == 1
