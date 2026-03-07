"""
weights.py — Custom weight/SF definitions beyond pocket-coffea defaults.

pocket-coffea already handles:
  - Pileup reweighting (from correctionlib)
  - b-tagging SFs (DeepJet, shape method)
  - Muon trigger / ID / ISO SFs
  - JEC/JER

Add custom weights here by subclassing WeightWrapper.
Each gets registered in analysis_config.py under the weights_cfg block.
"""

import numpy as np
import awkward as ak
from pocket_coffea.lib.weights.weights import WeightWrapper


class TopPtReweighting(WeightWrapper):
    """
    Top pT reweighting to correct NLO generator mis-modelling.
    Uses the CMS-recommended function:
        SF(top) = exp(a + b * pT_top)
    where a = 0.0615, b = -0.0005  (Run II tuning)

    Applied as sqrt(SF_top * SF_antitop) at the event level.
    Only meaningful for TTbar samples — returns 1.0 for everything else.
    """

    name = "top_pt_reweighting"
    has_variations = True  # generates Up/Down = 2x and 1.0 variations

    def __init__(self, params, metadata):
        super().__init__(params, metadata)
        self._a = 0.0615
        self._b = -0.0005

    def compute(self, events, size, shape_variation):
        sample = self.metadata["sample"]

        # Only apply to ttbar samples
        if "TT" not in sample:
            ones = np.ones(size)
            return ones, ones, ones   # nominal, up, down

        # Find gen-level top and antitop
        gen = events.GenPart
        is_top     = (np.abs(gen.pdgId) == 6) & gen.hasFlags(["isLastCopy"])
        top_pts    = gen[is_top & (gen.pdgId == 6)].pt
        antitop_pts = gen[is_top & (gen.pdgId == -6)].pt

        # Take first (should be exactly one each in ttbar events)
        pt_top     = ak.firsts(top_pts,    axis=1)
        pt_antitop = ak.firsts(antitop_pts, axis=1)

        # Fill missing with 0 (edge case: outside acceptance)
        pt_top     = ak.fill_none(pt_top,     0.0)
        pt_antitop = ak.fill_none(pt_antitop, 0.0)

        def sf(pt):
            return np.exp(self._a + self._b * np.clip(ak.to_numpy(pt), 0, 500))

        sf_top     = sf(pt_top)
        sf_antitop = sf(pt_antitop)
        weight_nom = np.sqrt(sf_top * sf_antitop)

        # Up: apply twice; Down: do not apply (weight = 1)
        weight_up   = weight_nom ** 2
        weight_down = np.ones(size)

        return weight_nom, weight_up, weight_down


class LHEScaleVariation(WeightWrapper):
    """
    QCD scale (muR, muF) variations from LHE weights.
    Indices follow the standard CMS 9-point grid:
        0: muR=0.5 muF=0.5
        1: muR=0.5 muF=1
        4: nominal (muR=1 muF=1)
        8: muR=2   muF=2
    We take the envelope of the 6 non-diagonal variations as Up/Down.
    """

    name = "qcd_scale"
    has_variations = True

    # Indices to use (exclude diagonal extremes 0 and 8 per CMS convention)
    INDICES = [1, 2, 3, 5, 6, 7]
    NOMINAL_IDX = 4

    def compute(self, events, size, shape_variation):
        if not hasattr(events, "LHEScaleWeight"):
            ones = np.ones(size)
            return ones, ones, ones

        lhe = ak.to_numpy(events.LHEScaleWeight)
        if lhe.shape[1] < 9:
            ones = np.ones(size)
            return ones, ones, ones

        nominal = lhe[:, self.NOMINAL_IDX]
        variations = lhe[:, self.INDICES]

        # Normalise to nominal
        ratios = variations / nominal[:, None]
        weight_up   = np.max(ratios,  axis=1)
        weight_down = np.min(ratios,  axis=1)

        return np.ones(size), weight_up, weight_down
