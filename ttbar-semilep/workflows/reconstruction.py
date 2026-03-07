"""
reconstruction.py — Top quark pair reconstruction in the semileptonic channel.

Strategy (simple chi-squared kinematic fit):

1. Neutrino pz from W mass constraint (quadratic — take real solution
   closest to muon pz, or smaller |pz| if complex).
2. Hadronic top: choose 3-jet combination from untagged jets + 1 b-jet
   that minimises |M_jjb - m_top|.
3. Leptonic top: mu + nu_reco + remaining b-jet.

This is intentionally simple and fast.  A full KinFitter is better for
precision but overkill for a first pass.
"""

import numpy as np
import awkward as ak
from coffea.nanoevents.methods import vector

# Attach vector behavior
ak.behavior.update(vector.behavior)

M_W   = 80.377   # GeV
M_TOP = 172.69   # GeV


def neutrino_pz(muon, met):
    """
    Reconstruct neutrino pz from the W mass constraint.
    Returns (pz_1, pz_2) — caller picks the physical solution.
    Both are ak.Arrays of shape (n_events,).
    """
    mu_pt  = muon.pt
    mu_eta = muon.eta
    mu_phi = muon.phi
    mu_e   = muon.energy

    nu_px = met.pt * np.cos(met.phi)
    nu_py = met.pt * np.sin(met.phi)

    # Scalar products
    mu_px = mu_pt * np.cos(mu_phi)
    mu_py = mu_pt * np.sin(mu_phi)
    mu_pz = mu_pt * np.sinh(mu_eta)

    A = M_W ** 2 / 2.0 + mu_px * nu_px + mu_py * nu_py
    a_coeff = mu_e ** 2 - mu_pz ** 2
    b_coeff = -2.0 * A * mu_pz
    c_coeff = mu_e ** 2 * met.pt ** 2 - A ** 2

    discriminant = b_coeff ** 2 - 4.0 * a_coeff * c_coeff

    # For complex solutions: set discriminant = 0 (real part only)
    disc_safe = ak.where(discriminant < 0, 0.0, discriminant)

    pz1 = (-b_coeff + np.sqrt(disc_safe)) / (2.0 * a_coeff)
    pz2 = (-b_coeff - np.sqrt(disc_safe)) / (2.0 * a_coeff)

    return pz1, pz2


def build_neutrino(met, pz):
    """Build a 4-vector for the neutrino given MET (px, py) and pz."""
    px = met.pt * np.cos(met.phi)
    py = met.pt * np.sin(met.phi)
    e  = np.sqrt(px**2 + py**2 + pz**2)
    return ak.zip(
        {"pt": met.pt, "eta": np.arcsinh(pz / met.pt),
         "phi": met.phi, "energy": e},
        with_name="Momentum4D",
        behavior=vector.behavior,
    )


def reconstruct_tops(events):
    """
    Full semileptonic ttbar reconstruction.

    Requires events to have:
        events.GoodMuons   (exactly 1 entry per event)
        events.GoodJets    (>= 4)
        events.BJets       (>= 2; subset of GoodJets)
        events.MET

    Returns a dict of scalar ak.Arrays (one value per event):
        reco_top_mass, reco_antitop_mass, reco_W_mass,
        reco_top_pt, chi2
    """
    muon = ak.firsts(events.GoodMuons)
    met  = events.MET

    # ── Neutrino pz ──────────────────────────────────────────────────────────
    pz1, pz2 = neutrino_pz(muon, met)
    # Pick solution with smaller |pz| (standard choice)
    pz_nu = ak.where(np.abs(pz1) < np.abs(pz2), pz1, pz2)
    nu = build_neutrino(met, pz_nu)

    # Leptonic W
    lep_W = muon + nu

    # ── Hadronic top reco ─────────────────────────────────────────────────────
    # Use b-jets as b-candidates.  Take the 2 highest-pT b-jets.
    bjets = events.BJets[:, :2]   # shape (n_events, <=2)

    # Light jets = good jets not in the top-2 b-jets
    # (pocket-coffea doesn't provide a "not-b" collection, so we
    #  reconstruct it from indices)
    b_idx   = ak.local_index(bjets,   axis=1)
    all_idx = ak.local_index(events.GoodJets, axis=1)

    # Simple proxy: use the first 2 non-bjet jets as light jets
    # A proper combinatoric search is done below when >=3 light jets exist.
    light_jets = events.GoodJets[events.GoodJets.btagDeepFlavB
                                  < 0.2783][:, :4]

    # For each event, try all pairs of light jets and pick the pair whose
    # invariant mass is closest to M_W, then add the hadronic b-jet.
    # We use the highest-pT b-jet for the hadronic side as a first pass.
    b_had  = bjets[:, 0:1]   # (n_events, 1)
    b_lep  = bjets[:, 1:2]   # (n_events, 1)

    # Build all (i,j) pairs of light jets for each event
    lj_i, lj_j = ak.unzip(ak.argcartesian([light_jets, light_jets], axis=1))
    unique = lj_i < lj_j
    pairs_i = light_jets[lj_i[unique]]
    pairs_j = light_jets[lj_j[unique]]
    dijet_mass = (pairs_i + pairs_j).mass

    # Pick best pair (closest to M_W)
    best_pair_idx = ak.argmin(np.abs(dijet_mass - M_W), axis=1, keepdims=True)
    best_lj_i = pairs_i[best_pair_idx]
    best_lj_j = pairs_j[best_pair_idx]

    reco_W_had = ak.firsts(best_lj_i + best_lj_j)
    reco_top_had = reco_W_had + ak.firsts(b_had)

    # Leptonic top
    reco_top_lep = lep_W + ak.firsts(b_lep)

    # chi2 for quality metric
    chi2 = (
        (ak.fill_none(reco_W_had.mass,  0.0) - M_W)   ** 2 / (10.0 ** 2) +
        (ak.fill_none(reco_top_had.mass, 0.0) - M_TOP) ** 2 / (14.0 ** 2)
    )

    return {
        "RecoTopMass":    ak.fill_none(reco_top_had.mass, -1.0),
        "RecoAntiTopMass": ak.fill_none(reco_top_lep.mass, -1.0),
        "RecoWMass":      ak.fill_none(reco_W_had.mass,   -1.0),
        "RecoTopPt":      ak.fill_none(reco_top_had.pt,   -1.0),
        "RecoChi2":       chi2,
    }
