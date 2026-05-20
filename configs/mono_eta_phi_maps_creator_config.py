from lbl_params import *

nEvents = -1

base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/"

inputFilePath = f"{base_path}/collisionData/initial_singleEG5/ntuple_0.root"
histogramsOutputFilePath = "../eta_phi_maps.root"

nBins = 1000

histParams2D = (
    ("photon_eta_phi", nBins, -4, 4, nBins, -4, 4, ""),
    ("goodPhoton_eta_phi", nBins, -4, 4, nBins, -4, 4, ""),
    ("caloTower_eta_phi", nBins, -4, 4, nBins, -4, 4, ""),
    ("goodCaloTower_eta_phi", nBins, -4, 4, nBins, -4, 4, ""),
)

eventsTreeNames = ["Events"]
specialBranchSizes = {}

