# corrections/

This directory holds calibration files used by the processor.
They are **not committed to git** (too large); fetch them from CVMFS or
the CMS GitLab mirror and place them here.

## Required files

### Golden JSON (data)
```
Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
```
Source: `cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/`

### Pileup weights
```
puWeights_UL18.json.gz
```
Source: `cms-nanoAOD.web.cern.ch/cms-nanoaod-integration/corrections/`

### b-tagging (DeepJet)
```
DeepJet_106XUL18SF_WPonly_V2p1.csv   ← for cut-based
DeepJet_106XUL18SF_shape_v3.csv      ← for shape-based (recommended)
```
Source: `twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL18`

### Muon SFs (trigger / ID / ISO)
```
Muon_UL2018_ID.json.gz
Muon_UL2018_ISO.json.gz
Muon_UL2018_Trigger.json.gz
```
Source: `cms-nanoAOD.web.cern.ch/cms-nanoaod-integration/corrections/`
(All in correctionlib JSON format, compatible with `correctionlib>=2.4`)

### JEC tarballs
```
Summer19UL18_V5_MC_AK4PFchs.tar.gz    ← MC JEC
Summer19UL18_RunA_V5_DATA_AK4PFchs.tar.gz
...etc...
```
Source: `cms-jerc.web.cern.ch/cms-jerc/data/` or CVMFS

## Quick fetch script (run once on each machine)

```bash
# From within CMSSW or any machine with xrootd access
cd corrections/

# Golden JSON
wget https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt

# Correctionlib JSON bundles (includes pileup, muon SFs, b-tagging)
wget https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/MUO_2018_UL_muon.json
```

On CVMFS (available at lxplus and LPC):
```bash
ls /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/
```
