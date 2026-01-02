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
  ("goodPhoton", "SCEnergy", 1000, 0, 1000, ""),
  ("goodPhoton", "SCEt", 1000, 0, 1000, ""),
  ("goodPhoton", "SCEta", 100, -5, 5, ""),
  ("goodPhoton", "SCEtaWidth", 100, 0, 0.1, ""),
  ("goodPhoton", "SCPhi", 100, -5, 5, ""),
  ("goodPhoton", "SCPhiWidth", 1000, 0, 0.2, ""),
  ("goodPhoton", "energy", 1000, 0, 1000, ""),
  ("goodPhoton", "energyBottom", 1000, 0, 1000, ""),
  ("goodPhoton", "energyLeft", 1000, 0, 1000, ""),
  ("goodPhoton", "energyRight", 1000, 0, 1000, ""),
  ("goodPhoton", "energyTop", 1000, 0, 1000, ""),
  # ("goodPhoton", "et", 1000, 0, 1000, ""),  # removing from here for blinding (custom filling)
  ("goodPhoton", "eta", 100, -5, 5, ""),
  ("goodPhoton", "hOverE", 1000, 0, 1.0, ""),
  ("goodPhoton", "hasConversionTracks", 2, 0, 1, ""),
  ("goodPhoton", "maxEnergyCrystal", 1000, 0, 1000, ""),
  ("goodPhoton", "phi", 100, -5, 5, ""),
  # ("goodPhoton", "seedTime", 100, -5, 5, ""),  # removing from here for blinding (custom filling)
  ("goodPhoton", "sigmaEta2012", 100, 0, 0.1, ""),
  ("goodPhoton", "sigmaIEtaIEta2012", 100, 0, 0.1, ""),

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
  ("unfoldingPhoton", "pt", 5, 0, 1, ""),
  ("unfoldingPhoton", "mass", 5, 5, 25, ""),
  ("unfoldingPhoton", "absRap", 2, 0, 2.2, ""),
  ("unfoldingPhoton", "absRap3", 3, 0, 2.2, ""),
  ("unfoldingPhoton", "rap3", 3, -2.2, 2.2, ""),
  ("unfoldingPhoton", "rap4", 4, -2.2, 2.2, ""),
  ("unfoldingPhoton", "costhetastar2", 2, 0, 1, ""),
  ("unfoldingPhoton", "costhetastar3", 3, 0, 1, ""),
  ("unfoldingPhoton", "costhetastar4", 4, 0, 1, ""),
  ("goodPhotonSR", "et", 5, 2, 8, ""),
  ("goodPhotonSR", "eta", 5, -2.2, 2.2, ""),
  ("goodPhotonSR", "phi", 6, -3.14, 3.14, ""),
  ("goodPhoton", "et", 2000, 0, 1000, ""),
  ("goodPhoton", "logEt", 200, -1, 3, ""),
  ("goodPhoton", "seedTime", 1000, -50, 50, ""),
  ("goodPhoton", "topOverCentral", 1000, 0, 10, ""),
  ("goodPhoton", "bottomOverCentral", 1000, 0, 10, ""),
  ("goodPhoton", "leftOverCentral", 1000, 0, 10, ""),
  ("goodPhoton", "rightOverCentral", 1000, 0, 10, ""),
  ("goodPhoton", "minOverCentral", 1000, 0, 1.0, ""),
  ("goodPhoton", "verticalOverCentral", 10000, 0, 10, ""),
  ("goodPhoton", "horizontalOverCentral", 10000, 0, 10, ""),
  ("goodPhoton", "horizontalImbalance", 100, -2, 2, ""),
  ("goodPhoton", "verticalImbalance", 100, -2, 2, ""),

  # gen-level
  ("genPhoton", "et", 200, 0, 10, ""),
  ("genPhoton", "energy", 200, 0, 10, ""),
  ("leadingGenPhoton", "energy", 200, 0, 10, ""),
  ("leadingGenPhotonBarrel", "energy", 200, 0, 10, ""),
  ("leadingGenPhotonBarrelEndcap", "energy", 200, 0, 10, ""),

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

histParams2D = (
  ("goodPhoton_absEta_vs_et", 100, 0, 3.0, 1000, 0, 1000, ""),
  ("goodPhoton_eta_vs_et", 100, -3.0, 3.0, 1000, 0, 1000, ""),
  ("goodPhoton_eta_vs_phi", 1000, -3.0, 3.0, 1000, -4.0, 4.0, ""),
  ("goodPhoton_eta_vs_phi_vs_et", 1000, -3.0, 3.0, 1000, -4.0, 4.0, ""),
  ("goodPhoton_eta_vs_phi_gt50GeV", 100, -3.0, 3.0, 100, -4.0, 4.0, ""),
  ("egamma_et_vs_goodPhoton_et", 1000, 0, 1000, 1000, 0, 1000, ""),
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
