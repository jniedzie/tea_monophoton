#include "LbLSelections.hpp"

#include "ConfigManager.hpp"
#include "Logger.hpp"
#include "UserExtensionsHelpers.hpp"

using namespace std;

LbLSelections::LbLSelections() {
  auto& config = ConfigManager::GetInstance();
  config.GetMap("eventCuts", eventCuts);

  lblObjectsManager = make_shared<LbLObjectsManager>();
}

bool LbLSelections::PassesNeutralExclusivity(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  auto towers = event->GetCollection("CaloTower");
  int nPassingTowers = 0;

  for (auto physicsObject : *towers) {
    auto tower = asCaloTower(physicsObject);
    string detectorType = tower->GetDetectorType();

    // We skip these checks to make NEE super-clean
    // if (detectorType == "HadronicCrack" || detectorType == "ElectromagneticCrack") continue;
    // if (tower->IsDead()) continue;

    if (detectorType == "HCAL") {
      // add hotspots (make maps, see if there are any)
      if (tower->IsHadronicEnergyAboveNoiseThreshold()) nPassingTowers++;
      if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;
    }
    else if (detectorType == "ECAL") {
      // add hotspots (?)
      
      // Skipping for super-clean NEE (to be replaced with eta-dependent thresholds)
      // if (tower->IsEtaAboveLimit()) continue;
      if (tower->OverlapsWithOtherObjects(event->GetCollection("goodPhoton"))) continue;
      if (tower->IsElectromagneticEnergyAboveNoiseThreshold()) nPassingTowers++;
      if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;
    }
    else if (detectorType == "HF") {
      // add hotspots (?)
      
      if (tower->IsHadronicEnergyAboveNoiseThreshold()) nPassingTowers++;
      if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;
    }
    else {
      fatal() << "Unknown detector type for calo tower: " << detectorType << endl;
      exit(0);
    }
  }

  if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;

  if (cutFlowManager) cutFlowManager->UpdateCutFlow("neutralExclusivity");

  return true;
}

bool LbLSelections::PassesPhotonSelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  auto goodPhotons = event->GetCollection("goodPhoton");
  int nPhotons = goodPhotons->size();
  if (nPhotons < eventCuts.at("min_Nphotons") || nPhotons > eventCuts.at("max_Nphotons")) return false;

  cutFlowManager->UpdateCutFlow("nPhotons");
  return true;
}

bool LbLSelections::PassesChargedExclusivity(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  
  try {
    int nAdditionalMuonSegmentsCSC = event->Get("nAdditionalMuonCSCsegments");
    if (nAdditionalMuonSegmentsCSC > eventCuts.at("max_NmuonSegmentsCSC")) return false;
    if (cutFlowManager) cutFlowManager->UpdateCutFlow("nMuonSegmentsCSC");
  }
  catch (Exception& e) {
    error() << "nAdditionalMuonCSCsegments not found in event. Skipping this cut." << endl;
  }
  
  int nElectrons = event->GetCollection("goodElectron")->size();
  if (nElectrons > eventCuts.at("max_Nelectrons") || nElectrons < eventCuts.at("min_Nelectrons")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nElectrons");

  auto tracks = event->GetCollection("goodTrack");
  auto electrons = event->GetCollection("goodElectron");
  int nNonOverlappingTracks = 0;
  for (auto physicsObject : *tracks) {
    auto track = asTrack(physicsObject);
    if (!track->OverlapsWithOtherObjects(electrons)) nNonOverlappingTracks++;
  }
  if (nNonOverlappingTracks > eventCuts.at("max_Ntracks")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nTracks");

  int nMuons = event->GetCollection("goodMuon")->size();
  if (nMuons > eventCuts.at("max_Nmuons")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nMuons");

  int nStandaloneMuons = event->GetCollection("standaloneMuon")->size();
  if (nStandaloneMuons > eventCuts.at("max_NstandaloneMuons")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nStandaloneMuons");

  return true;
}

bool LbLSelections::PassesBeamHaloFilters(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  map<int, string> options = {
    {0, "none"},
    {1, "isBeamHaloLoose"},
    {2, "isBeamHaloTight"},
    {3, "isBeamHaloGlobalTight2016"},
    {4, "isbeamHaloGlobalSuperTight2016"},
  };
  int option = eventCuts.at("beamHaloFilter");
  
  if (options.find(option) == options.end()) {
    fatal() << "Invalid beam halo filter option: " << option << ". Valid options are: " << endl;
    for (const auto& opt : options) {
      fatal() << opt.first << ": " << opt.second << endl;
    }
    exit(0);
  }

  string beamHaloOption = options.at(option);
  if (beamHaloOption == "none") {
    cutFlowManager->UpdateCutFlow("beamHaloFilters");
    return true;
  }
  
  if (event->GetAs<bool>(beamHaloOption)) return false;
  cutFlowManager->UpdateCutFlow("beamHaloFilters");
  
  return true;
}

bool LbLSelections::PassesZDC(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  if (eventCuts.at("ZDC_cut") == 0) {
    cutFlowManager->UpdateCutFlow("ZDC");
    return true;
  }

  // If ZDC cuts are applied, we need to skip the bad runs
  int badRunsMin = 326571;
  int badRunsMax = 326676;
  int runNumber = event->GetAs<int>("runNumber");
  if (runNumber >= badRunsMin && runNumber <= badRunsMax) return false;

  shared_ptr<PhysicsObjects> zdcEnergies;

  try {
    zdcEnergies = event->GetCollection("ZDC");
  } catch (Exception& e) {
    warn() << "No ZDC collection found in event. Will skip ZDC cuts." << endl;
    cutFlowManager->UpdateCutFlow("ZDC");
    return true;
  }

  float totalEnergyPlus = 0;
  float totalEnergyMinus = 0;

  for (auto physicsObject : *zdcEnergies) {
    auto zdcEnergy = asZDCEnergy(physicsObject);

    if (zdcEnergy->GetSide() > 0) {
      totalEnergyPlus += zdcEnergy->GetEnergy();
    } else if (zdcEnergy->GetSide() < 0) {
      totalEnergyMinus += zdcEnergy->GetEnergy();
    }
  }

  // <1n: 1600, <2n: 4000, <3n: 7000 GeV, <4n: 10000 GeV
  if (eventCuts.at("ZDC_cut") == 1) {
    if (totalEnergyPlus > 7000 && totalEnergyMinus > 7000) return false;
  } else if (eventCuts.at("ZDC_cut") == 2) {
    if (totalEnergyPlus > 1600 || totalEnergyMinus > 1600) return false;  // 0n0n
  } else if (eventCuts.at("ZDC_cut") == 3) {
    if (totalEnergyPlus > 4000 || totalEnergyMinus > 4000) return false;  // 0n0n + 1n0n + 0n1n + 1n1n
  } else if (eventCuts.at("ZDC_cut") == 4) {
    auto photons = event->GetCollection("goodPhoton");
    auto photon = asPhoton(photons->at(0));
    float photonEta = photon->GetEta();
    float photonEt = photon->GetEt();
    float zdcEnergy = (photonEta < 0) ? totalEnergyPlus : totalEnergyMinus;

    // if (photonEt > 30 && (zdcEnergy >= 10000)) {

    //   ofstream photonFile("/afs/desy.de/user/j/jniedzie/tea_monophoton/debugging/info" + to_string(event->GetAs<int>("eventNumber")) +
    //   ".txt"); if (photonFile.is_open()) {
    //     // print all photon and zdc info using their operator<< overloads
    //     photonFile << "Event: " << event->GetAs<int>("eventNumber") << ", Lumi section: " << event->GetAs<int>("lumiSection")
    //                << ", Run: " << event->GetAs<int>("runNumber") << endl;
    //     photonFile << "Photon information:\n" << *photon << endl;

    //     // photonFile << "Photon Et: " << photonEt << "\teta: " << photonEta << "\tphi: " << photon->GetPhi() << endl;
    //     photonFile << "\tZDC+: " << totalEnergyPlus << "\tZDC-: " << totalEnergyMinus << "\tzdcEnergy: " << zdcEnergy << endl;
    //   }
    // }

    if (zdcEnergy < 10000) return false;  // 0nXn, X>3n (on the side opposite to the photon)
    // if (totalEnergyPlus < 10000 || totalEnergyMinus < 10000) return false;  // 0nXn, X≥3n
  } else {
    warn() << "Unknown ZDC cut type: " << eventCuts.at("ZDC_cut") << ". Will skip ZDC cuts." << endl;
  }

  cutFlowManager->UpdateCutFlow("ZDC");

  return true;
}
