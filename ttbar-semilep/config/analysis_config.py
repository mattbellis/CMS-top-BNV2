import os
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import cloudpickle

from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_functions import (
    get_nPVgood, goldenJson, eventFlags,
    get_nObj_min, get_nBtagMin, get_HLTsel,
)
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *
from pocket_coffea.parameters import defaults
from pocket_coffea.lib.hist_manager import Axis, HistConf
from pocket_coffea.lib.columns_manager import ColOut

import workflows.ttbar_semilep_processor as workflow
from workflows.ttbar_semilep_processor import TTbarSemilepProcessor

cloudpickle.register_pickle_by_value(workflow)

localdir = os.path.dirname(os.path.abspath(__file__))
default_parameters = defaults.get_default_parameters()
parameters = default_parameters

FILELIST_DIR = os.path.join(localdir, "../datasets/filelists")

cfg = Configurator(
    workflow=TTbarSemilepProcessor,
    parameters=parameters,
    datasets={
        "jsons": [
            os.path.join(FILELIST_DIR, "mc_samples_UL18.json"),
            os.path.join(FILELIST_DIR, "data_samples_UL18.json"),
        ],
        #"filter": {
        #    "samples": [],
        #    "samples_exclude": [],
        #    
        #},
    },
    skim=[
        get_nPVgood(1),
        eventFlags,
        goldenJson,
    ],
    preselections=[passthrough],
    categories={
        "SR_2b": [
            get_nObj_min(4, coll="JetGood"),
            get_nBtagMin(2, coll="BJetGood"),
        ],
        "CR_1b": [
            get_nObj_min(4, coll="JetGood"),
            get_nBtagMin(1, coll="BJetGood"),
        ],
    },
    weights={
        "common": {
            "inclusive": ["genWeight", "lumi", "XS", "pileup"],
            "bycategory": {},
        },
        "bycategory": {},
        "bysample": {},
    },
    variations={
        "weights": {
            "common": {
                "inclusive": ["pileup"],
                "bycategory": {},
            },
            "bysample": {},
            "bycategory": {},
        },
    },
    variables={
        "muon_pt": HistConf(
            [Axis(coll="MuonGood", field="pt", bins=50, start=0, stop=300,
                  label=r"Muon $p_T$ [GeV]")]
        ),
        "muon_eta": HistConf(
            [Axis(coll="MuonGood", field="eta", bins=50, start=-2.5, stop=2.5,
                  label=r"Muon $\eta$")]
        ),
        "jet_pt": HistConf(
            [Axis(coll="JetGood", field="pt", bins=50, start=0, stop=500,
                  label=r"Jet $p_T$ [GeV]")]
        ),
        "jet_eta": HistConf(
            [Axis(coll="JetGood", field="eta", bins=50, start=-2.5, stop=2.5,
                  label=r"Jet $\eta$")]
        ),
        "met_pt": HistConf(
            [Axis(coll="MET", field="pt", bins=50, start=0, stop=400,
                  label=r"$p_T^{miss}$ [GeV]")]
        ),
    },
    columns={},
    workflow_options={},
)
