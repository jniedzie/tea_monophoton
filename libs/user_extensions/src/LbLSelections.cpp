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

    if (detectorType == "HCAL") {
      if (tower->IsDead()) continue;
      if (tower->IsInHadronicCrack()) continue;
      if (tower->IsHadronicEnergyAboveNoiseThreshold()) nPassingTowers++;
      if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;
    }
    if (detectorType == "ECAL") {
      if (tower->IsDead()) continue;
      if (tower->IsEtaAboveLimit()) continue;
      if (tower->IsInHEM()) continue;
      if (tower->IsInElectromagneticCrack()) continue;
      if (tower->OverlapsWithOtherObjects(event->GetCollection("goodPhoton"))) continue;
      if (tower->OverlapsWithOtherObjects(event->GetCollection("goodElectron"))) continue;
      if (tower->IsElectromagneticEnergyAboveNoiseThreshold()) nPassingTowers++;
      if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;
    }
    if (detectorType == "HF") {
      if (tower->IsDead()) continue;
      if (tower->OverlapsWithOtherObjects(event->GetCollection("goodPhoton"))) continue;
      if (tower->OverlapsWithOtherObjects(event->GetCollection("goodElectron"))) continue;

      // This for HF is actually equivalent to checking the EM energy:
      if (tower->IsHadronicEnergyAboveNoiseThreshold()) nPassingTowers++;
      if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;
    }
  }
  if (nPassingTowers > eventCuts.at("max_Ntowers")) return false;

  if (cutFlowManager) cutFlowManager->UpdateCutFlow("neutralExclusivity");

  return true;
}

bool LbLSelections::PassesZeroPhotonAndElectronSelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  auto goodPhotons = event->GetCollection("goodPhoton");
  int nPhotons = goodPhotons->size();
  if (nPhotons != 0) return false;

  auto electrons = event->GetCollection("goodElectron");
  int nElectrons = electrons->size();
  if (nElectrons != 0) return false;

  cutFlowManager->UpdateCutFlow("zeroPhotonElectron");
  return true;
}

bool LbLSelections::PassesSinglePhotonSelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  auto goodPhotons = event->GetCollection("goodPhoton");
  int nPhotons = goodPhotons->size();
  if (nPhotons != 1) return false;

  cutFlowManager->UpdateCutFlow("singlePhoton");
  return true;
}

bool LbLSelections::PassesThreePhotonsSelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  auto goodPhotons = event->GetCollection("goodPhoton");
  int nPhotons = goodPhotons->size();
  if (nPhotons != 3) return false;

  cutFlowManager->UpdateCutFlow("threePhotons");
  return true;
}

bool LbLSelections::PassesDiphotonSelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager,
                                            shared_ptr<map<string, int>> cutFlow) {
  auto genMatchedPhotons = asLbLEvent(event)->GetGenMatchedRecoPhotons();
  cutFlow->at("00_initial") += 2;
  cutFlow->at("01_genMatched") += genMatchedPhotons.size();

  if (genMatchedPhotons.size() != 2) {
    warn() << "Number of gen-matched reco photons != 2" << endl;
  } else {
    lblObjectsManager->IsGoodPhoton(asPhoton(genMatchedPhotons.at(0)), cutFlow);
    lblObjectsManager->IsGoodPhoton(asPhoton(genMatchedPhotons.at(1)), cutFlow);
  }

  auto goodPhotons = event->GetCollection("goodPhoton");
  int nPhotons = goodPhotons->size();
  if (nPhotons != 2) return false;

  cutFlowManager->UpdateCutFlow("twoGoodPhotons");

  if (nPhotons < 2) {
    fatal() << "Requested diphoton selections, but <2 good photons found in event." << endl;
    exit(0);
  }
  auto photon1 = goodPhotons->at(0);
  auto photon2 = goodPhotons->at(1);
  TLorentzVector photon1vec, photon2vec;
  photon1vec.SetPtEtaPhiM(photon1->Get("et"), photon1->Get("eta"), photon1->Get("phi"), 0);
  photon2vec.SetPtEtaPhiM(photon2->Get("et"), photon2->Get("eta"), photon2->Get("phi"), 0);

  auto diphoton = photon1vec + photon2vec;
  if (diphoton.M() < eventCuts.at("min_diphotonMass")) return false;
  cutFlowManager->UpdateCutFlow("diphotonMass");
  return true;
}

bool LbLSelections::PassesDielectronSelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager, bool sameCharge) {
  auto electrons = event->GetCollection("goodElectron");

  // check number of electrons
  int nElectrons = electrons->size();
  if (nElectrons != 2) return false;
  cutFlowManager->UpdateCutFlow("twoGoodElectrons");

  // check electron charge
  auto electron1 = electrons->at(0);
  auto electron2 = electrons->at(1);
  if (((int)electron1->Get("charge") * (int)electron2->Get("charge") > 0) && !sameCharge) return false;
  if (((int)electron1->Get("charge") * (int)electron2->Get("charge") < 0) && sameCharge) return false;
  cutFlowManager->UpdateCutFlow("electronCharge");

  // check charge exclusivity
  if (!PassesDielectronChargedExclusivity(event, cutFlowManager)) return false;

  TLorentzVector electron1vec, electron2vec;
  electron1vec.SetPtEtaPhiM(electron1->Get("pt"), electron1->Get("eta"), electron1->Get("phi"), 0);
  electron2vec.SetPtEtaPhiM(electron2->Get("pt"), electron2->Get("eta"), electron2->Get("phi"), 0);

  // check dielectron mass
  auto dielectron = electron1vec + electron2vec;
  if (dielectron.M() < eventCuts.at("min_dielectronMass")) return false;
  cutFlowManager->UpdateCutFlow("dielectronMass");

  // check dielectron pt
  if (dielectron.Pt() > eventCuts.at("max_dielectronPt")) return false;
  cutFlowManager->UpdateCutFlow("dielectronPt");

  return true;
}

bool LbLSelections::PassesDielectronChargedExclusivity(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  auto electrons = event->GetCollection("goodElectron");
  auto tracks = event->GetCollection("goodTrack");

  int nNonOverlappingTracks = 0;
  for (auto physicsObject : *tracks) {
    auto track = asTrack(physicsObject);
    if (!track->OverlapsWithOtherObjects(electrons)) nNonOverlappingTracks++;
  }
  if (nNonOverlappingTracks > eventCuts.at("max_Ntracks")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nTracks");

  return true;
}

bool LbLSelections::PassesChargedExclusivity(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  int nElectrons = event->GetCollection("goodElectron")->size();
  if (nElectrons > eventCuts.at("max_Nelectrons")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nElectrons");

  int nTracks = event->GetCollection("goodTrack")->size();
  if (nTracks > eventCuts.at("max_Ntracks")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nTracks");

  int nMuons = event->GetCollection("goodMuon")->size();
  if (nMuons > eventCuts.at("max_Nmuons")) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("nMuons");

  return true;
}

bool LbLSelections::PassesDiphotonPt(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  if (event->GetCollection("goodPhoton")->size() != 2) return false;

  auto photon1 = event->GetCollection("goodPhoton")->at(0);
  auto photon2 = event->GetCollection("goodPhoton")->at(1);
  TLorentzVector photon1vec, photon2vec;
  photon1vec.SetPtEtaPhiM(photon1->Get("et"), photon1->Get("eta"), photon1->Get("phi"), 0);
  photon2vec.SetPtEtaPhiM(photon2->Get("et"), photon2->Get("eta"), photon2->Get("phi"), 0);

  auto diphoton = photon1vec + photon2vec;
  if (diphoton.Pt() > eventCuts.at("max_diphotonPt")) return false;
  cutFlowManager->UpdateCutFlow("diphotonPt");

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

    //   ofstream photonFile("/afs/desy.de/user/j/jniedzie/tea_monophoton/debugging/info" + to_string(event->GetAs<int>("eventNumber")) + ".txt");
    //   if (photonFile.is_open()) {
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

bool LbLSelections::PassesTracksPlusPhotonsSelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  if (!PassesDiphotonSelection(event, cutFlowManager)) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("twoGoodPhotons");

  auto photons = event->GetCollection("goodPhoton");
  auto electrons = event->GetCollection("goodElectron");
  auto tracks = event->GetCollection("track");
  auto photon1 = asPhoton(photons->at(0));
  auto photon2 = asPhoton(photons->at(1));

  shared_ptr<Track> track1 = nullptr;
  shared_ptr<Track> track2 = nullptr;

  for (auto physicsObject : *tracks) {
    auto track = asTrack(physicsObject);
    if (track->OverlapsWithOtherObjects(photons) && !track->OverlapsWithOtherObjects(electrons)) {
      if (track1 == nullptr)
        track1 = track;
      else if (track2 == nullptr)
        track2 = track;
      else
        return false;
    }
  }

  if (track1 == nullptr || track2 == nullptr) return false;
  if (cutFlowManager) cutFlowManager->UpdateCutFlow("twoGoodTracks");

  return true;
}

bool LbLSelections::PassesAcoplanaritySelection(shared_ptr<Event> event, shared_ptr<CutFlowManager> cutFlowManager) {
  if (event->GetCollection("goodPhoton")->size() != 2) return false;

  auto photon1 = event->GetCollection("goodPhoton")->at(0);
  auto photon2 = event->GetCollection("goodPhoton")->at(1);
  TLorentzVector photon1vec, photon2vec;
  photon1vec.SetPtEtaPhiM(photon1->Get("et"), photon1->Get("eta"), photon1->Get("phi"), 0);
  photon2vec.SetPtEtaPhiM(photon2->Get("et"), photon2->Get("eta"), photon2->Get("phi"), 0);

  double deltaPhi = fabs(photon1vec.DeltaPhi(photon2vec));
  double acoplanarity = 1 - (deltaPhi / TMath::Pi());
  if (acoplanarity > eventCuts.at("max_diphotonAcoplanarity")) return false;
  cutFlowManager->UpdateCutFlow("diphotonAcoplanarity");

  return true;
}