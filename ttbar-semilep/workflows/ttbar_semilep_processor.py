"""
ttbar_semilep_processor.py — Main processor for the ttbar semileptonic analysis.

Inherits from pocket-coffea's BaseProcessorABC which provides:
  - Automatic cutflow accumulation
  - Weight / SF application (via weights_cfg in analysis_config.py)
  - Systematic variation branching
  - Histogram filling infrastructure
  - Column (parquet) output

You only need to override the object preselection and optionally
add extra event-level quantities before histograms are filled.
"""

import numpy as np
import awkward as ak

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.lib.objects import (
    jet_correction,         # applies JEC/JER and returns corrected jets
    lepton_selection,
    get_dilepton,
)

from config.parameters.cuts import (
    select_good_muons,
    select_veto_muons,
    select_veto_electrons,
    select_good_jets,
    select_btagged_jets,
    OBJECT_PARAMS,
)
from workflows.reconstruction import reconstruct_tops


class TTbarSemilepProcessor(BaseProcessorABC):
    """
    Processor for ttbar → (W→μν)(W→qq̄)bb.

    The BaseProcessorABC handles:
      - Reading NanoAOD branches
      - Applying golden JSON (data)
      - Looping over systematic variations
      - Filling histograms defined in analysis_config.py
      - Saving columnar outputs

    We override three key methods:
      1. apply_object_preselection  — build GoodMuons, GoodJets, BJets, …
      2. count_objects              — attach scalar counts for histogram axes
      3. define_extra_variables     — attach reco quantities (top mass, HT, …)
    """

    # ── 1. Object preselection ────────────────────────────────────────────────

    def apply_object_preselection(self, variation):
        """
        Called once per systematic variation (jet variations rerun jet
        selection; lepton variations rerun lepton SFs but not collection).

        Attaches collections to self.events:
            GoodMuons, VetoMuons, VetoElectrons,
            GoodJets, BJets
        """
        events = self.events

        # Apply JEC/JER corrections (pocket-coffea handles variation branching)
        if variation in ("nominal", "JES_up", "JES_down", "JER_up", "JER_down"):
            events["Jet"] = jet_correction(
                events, variation, self._year,
                self.params.get("jet_corrections", {}),
            )

        # Muons
        events["GoodMuons"]      = select_good_muons(events)
        events["VetoMuons"]      = select_veto_muons(events)
        events["VetoElectrons"]  = select_veto_electrons(events)

        # Jets (DR-cleaned vs GoodMuons, so muons must come first)
        events["GoodJets"] = select_good_jets(events)
        events["BJets"]    = select_btagged_jets(
            events.GoodJets,
            wp=OBJECT_PARAMS["jet"]["btag_wp"],
        )

    # ── 2. Count objects ──────────────────────────────────────────────────────

    def count_objects(self, variation):
        """
        Attach scalar-per-event counts as new fields on events.
        These are used as histogram axes (see analysis_config.py variables).
        """
        events = self.events
        events["nGoodMuons"]  = ak.num(events.GoodMuons)
        events["nGoodJets"]   = ak.num(events.GoodJets)
        events["nBJets"]      = ak.num(events.BJets)
        events["nLightJets"]  = events["nGoodJets"] - events["nBJets"]

    # ── 3. Extra variables (reco, HT, …) ─────────────────────────────────────

    def define_extra_variables(self):
        """
        Called after event selection is applied but before histograms are
        filled.  Attach per-event derived quantities.

        Only runs on events that pass the full selection (self.events is
        already filtered at this point).
        """
        events = self.events

        # HT = scalar sum of jet pTs
        events["HT"] = ak.sum(events.GoodJets.pt, axis=1)

        # Top quark reconstruction
        reco = reconstruct_tops(events)
        for key, val in reco.items():
            events[key] = val

    # ── Optional: extra diagnostics saved per-chunk ───────────────────────────

    def process_extra(self, events):
        """
        Hook called at the end of each chunk.  Use this to accumulate
        anything not handled by the standard histogram / column infrastructure.

        Here we store a small summary dict that dataset_diagnostics.py
        will collect into the CSV.
        """
        return {
            "n_events_processed": len(events),
            "n_events_passing":   int(ak.sum(self._selection_mask)),
        }
