from scale_factors_config import *
from lbl_params import *
import numpy as np
from numpy import pi

nEvents = -1
printEveryNevents = 1000

base_path = "/data/dust/user/jniedzie/monophoton/"

trigger = "doubleEG2"
# trigger = "singleEG5"

sample = "collisionData"
# sample = "lbl"
# sample = "cep"
# sample = "qed"
# sample = "qed_MG_ee_a"
# sample = "emptyBeams"
# sample = "qed_sc_noPhotos"

skim = f"skimmed_baselineSelections_{trigger}"

inputFilePath = f"{base_path}/{sample}/{skim}/ntuple_0.root"
histogramsOutputFilePath = f"../{skim}_{sample}_histograms.root"

defaultHistParams = (
  # collection      variable          bins    xmin     xmax     dir
  ("convertedPhoton", "et", 200, 0, 100, ""),
  ("convertedPhoton", "eta", 100, -2.2, 2.2, ""),
  ("convertedPhoton", "phi", 100, -3.14, 3.14, ""),
  ("unconvertedPhoton", "et", 200, 0, 100, ""),
  ("unconvertedPhoton", "eta", 100, -2.2, 2.2, ""),
  ("unconvertedPhoton", "phi", 100, -3.14, 3.14, ""),
  ("conversionElectron", "pt", 200, 0, 50, ""),
  ("conversionElectron", "eta", 100, -2.2, 2.2, ""),
  ("conversionElectron", "phi", 100, -3.14, 3.14, ""),
  ("conversionElectron", "nMissHits", 20, 0, 20, ""),
  ("conversionElectron", "hOverE", 100, 0, 0.5, ""),
  ("conversionElectron", "deltaEtaAtVertex", 100, 0, 0.5, ""),
  ("Event", "runNumber", 1000, 320000, 335000, ""),
)

histParams = (
  # photons

  # gen-level
  ("genPhoton", "et", 200, 0, 10, ""),
  ("genPhoton", "energy", 200, 0, 10, ""),

  # event
  ("event", "deltaEt", 100, 0, 1, ""),
  ("event", "ZDCenergyPlus", 10000, 0, 20000, ""),
  ("event", "ZDCenergyMinus", 10000, 0, 20000, ""),
  ("monophoton", "egamma_deltaEta", 1000, -10, 10, ""),
  ("monophoton", "egamma_deltaPhi", 1000, -10, 10, ""),
  ("monophoton", "egamma_deltaR", 1000, -10, 10, ""),
  ("monophoton", "egamma_deltaEta_gt50GeV", 1000, -10, 10, ""),
  ("monophoton", "egamma_deltaPhi_gt50GeV", 1000, -10, 10, ""),
  ("monophoton", "egamma_deltaR_gt50GeV", 1000, -10, 10, ""),
)

histParams2D = ((
  "egamma_et_vs_goodPhoton_et", 1000, 0, 1000, 1000, 0, 1000, ""
), )

for prefix in ["", "Barrel_", "EndCap_"]:
  histParams += (
    ("goodPhoton", f"{prefix}et", 2000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}logEt", 200, -1, 3, ""),
    ("goodPhoton", f"{prefix}seedTime", 1000, -50, 50, ""),
    ("goodPhoton", f"{prefix}topOverCentral", 1000, 0, 10, ""),
    ("goodPhoton", f"{prefix}bottomOverCentral", 1000, 0, 10, ""),
    ("goodPhoton", f"{prefix}leftOverCentral", 1000, 0, 10, ""),
    ("goodPhoton", f"{prefix}rightOverCentral", 1000, 0, 10, ""),
    ("goodPhoton", f"{prefix}minOverCentral", 1000, 0, 1.0, ""),
    ("goodPhoton", f"{prefix}verticalOverCentral", 10000, 0, 10, ""),
    ("goodPhoton", f"{prefix}horizontalOverCentral", 10000, 0, 10, ""),
    ("goodPhoton", f"{prefix}horizontalImbalance", 100, -2, 2, ""),
    ("goodPhoton", f"{prefix}verticalImbalance", 100, -2, 2, ""),
    
    ("goodPhoton", f"{prefix}SCEnergy", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}SCEt", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}SCEta", 100, -5, 5, ""),
    ("goodPhoton", f"{prefix}SCEtaWidth", 1000, 0, 0.2, ""),
    ("goodPhoton", f"{prefix}SCPhi", 100, -5, 5, ""),
    ("goodPhoton", f"{prefix}SCPhiWidth", 1000, 0, 0.2, ""),
    ("goodPhoton", f"{prefix}energy", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}energyBottom", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}energyLeft", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}energyRight", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}energyTop", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}eta", 100, -5, 5, ""),
    ("goodPhoton", f"{prefix}hOverE", 1000, 0, 1.0, ""),
    ("goodPhoton", f"{prefix}hasConversionTracks", 2, 0, 1, ""),
    ("goodPhoton", f"{prefix}maxEnergyCrystal", 1000, 0, 1000, ""),
    ("goodPhoton", f"{prefix}phi", 100, -5, 5, ""),
    ("goodPhoton", f"{prefix}sigmaEta2012", 100, 0, 0.1, ""),
    ("goodPhoton", f"{prefix}sigmaIEtaIEta2012", 100, 0, 0.1, ""),
  )

  histParams2D += (
    (f"goodPhoton_{prefix}absEta_vs_et", 100, 0, 3.0, 1000, 0, 1000, ""),
    (f"goodPhoton_{prefix}eta_vs_et", 100, -3.0, 3.0, 1000, 0, 1000, ""),
    (f"goodPhoton_{prefix}eta_vs_phi", 1000, -3.0, 3.0, 1000, -4.0, 4.0, ""),
    (
      f"goodPhoton_{prefix}eta_vs_phi_vs_et", 1000, -3.0, 3.0, 1000, -4.0, 4.0,
      ""
    ),
    (
      f"goodPhoton_{prefix}eta_vs_phi_gt50GeV", 100, -3.0, 3.0, 100, -4.0, 4.0,
      ""
    ),
  )

eventsTreeNames = [
  "Events",
]
specialBranchSizes = {}

extraEventCollections = {
  "convertedPhoton": {
    "inputCollections": ("photon", ),
    "hasConversionTracks": True,
  },
  "unconvertedPhoton": {
    "inputCollections": ("photon", ),
    "hasConversionTracks": False,
  },
  "conversionElectron": {
    "inputCollections": ("electron", ),
    "conversionVeto": False,
  },
}
