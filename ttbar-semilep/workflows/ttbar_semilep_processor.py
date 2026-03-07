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
        events = self.events

        events["MuonGood"]     = select_good_muons(events)
        events["MuonVeto"]     = select_veto_muons(events)
        events["ElectronVeto"] = select_veto_electrons(events)

        # Jets (DR-cleaned vs GoodMuons, so muons must come first)
        events["JetGood"]      = select_good_jets(events)
        events["BJetGood"]         = select_btagged_jets(
            events.JetGood,
            wp=OBJECT_PARAMS["jet"]["btag_wp"],
        )


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
        events["HT"] = ak.sum(events.JetGood.pt, axis=1)

        # Top quark reconstruction
        reco = reconstruct_tops(events)
        for key, val in reco.items():
            events[key] = val

    ###############################################################

    def count_objects(self, variation):
        events = self.events
        events["nMuonGood"] = ak.num(events.MuonGood)
        events["nJetGood"]  = ak.num(events.JetGood)
        events["nBJetGood"] = ak.num(events.BJetGood)
        events["nLightJet"] = events["nJetGood"] - events["nBJetGood"]

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
