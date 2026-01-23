#include "LbLHistogramsFiller.hpp"

#include "ConfigManager.hpp"
#include "ExtensionsHelpers.hpp"
#include "UserExtensionsHelpers.hpp"

using namespace std;

LbLHistogramsFiller::LbLHistogramsFiller(shared_ptr<HistogramsHandler> histogramsHandler_) : histogramsHandler(histogramsHandler_) {
  // Create a config manager
  auto& config = ConfigManager::GetInstance();

  try {
    config.GetMap("dataBlinding", dataBlinding);
  } catch (const Exception& e) {
    warn() << "No data blinding parameters found. Will not apply any data blinding." << endl;
    dataBlinding["max_et"] = 9999;
  }

  eventProcessor = make_unique<EventProcessor>();
}

LbLHistogramsFiller::~LbLHistogramsFiller() {}

void LbLHistogramsFiller::FillMonoPhotonHistograms(const shared_ptr<Event> event) {
  auto photons = event->GetCollection("goodPhoton");
  auto photon = asPhoton(photons->at(0));

  if (asLbLEvent(event)->IsData() && photon->GetEt() > dataBlinding["max_et"]) return;
  FillMonoPhotonHistograms(photon);

  if (fabs(photon->GetEta()) > 1.2) {
    FillMonoPhotonHistograms(photon, "EndCap_");
  } else {
    FillMonoPhotonHistograms(photon, "Barrel_");
  }
}

void LbLHistogramsFiller::FillMonoPhotonHistograms(const shared_ptr<Photon> photon, string prefix) {
  histogramsHandler->Fill("goodPhoton_" + prefix + "absEta_vs_et", fabs(photon->GetEta()), photon->GetEt());
  histogramsHandler->Fill("goodPhoton_" + prefix + "eta_vs_et", photon->GetEta(), photon->GetEt());
  histogramsHandler->Fill("goodPhoton_" + prefix + "eta_vs_phi", photon->GetEta(), photon->GetPhi());
  histogramsHandler->Fill("goodPhoton_" + prefix + "et", photon->GetEt());
  histogramsHandler->Fill("goodPhoton_" + prefix + "logEt", TMath::Log10(photon->GetEt()));
  histogramsHandler->Fill("goodPhoton_" + prefix + "seedTime", photon->GetSeedTime());
  histogramsHandler->Fill("goodPhoton_" + prefix + "topOverCentral", photon->GetEnergyTop() / photon->GetEnergyCentral());
  histogramsHandler->Fill("goodPhoton_" + prefix + "bottomOverCentral", photon->GetEnergyBottom() / photon->GetEnergyCentral());
  histogramsHandler->Fill("goodPhoton_" + prefix + "leftOverCentral", photon->GetEnergyLeft() / photon->GetEnergyCentral());
  histogramsHandler->Fill("goodPhoton_" + prefix + "rightOverCentral", photon->GetEnergyRight() / photon->GetEnergyCentral());
  histogramsHandler->Fill("goodPhoton_" + prefix + "minOverCentral", photon->GetMinEnergy() / photon->GetEnergyCentral());
  histogramsHandler->Fill("goodPhoton_" + prefix + "verticalOverCentral", photon->GetVerticalOverCentralEnergy());
  histogramsHandler->Fill("goodPhoton_" + prefix + "horizontalOverCentral", photon->GetHorizontalOverCentralEnergy());
  histogramsHandler->Fill("goodPhoton_" + prefix + "horizontalImbalance", photon->GetHorizontalImbalance());
  histogramsHandler->Fill("goodPhoton_" + prefix + "verticalImbalance", photon->GetVerticalImbalance());
  histogramsHandler->Fill("goodPhoton_" + prefix + "swissCross", photon->GetSwissCross());
  histogramsHandler->Fill("goodPhoton_" + prefix + "horizontalImbalance_vs_seedTime", photon->GetHorizontalImbalance(), photon->GetSeedTime());
  histogramsHandler->Fill("goodPhoton_" + prefix + "verticalImbalance_vs_seedTime", photon->GetVerticalImbalance(), photon->GetSeedTime());

  vector<string> defaultBranches = {
      "SCEnergy",         "SCEt",       "SCEta",        "SCEtaWidth",        "SCPhi", "SCPhiWidth", "energy",
      "energyBottom",     "energyLeft", "energyRight",  "energyTop",         "eta",   "hOverE",     "hasConversionTracks",
      "maxEnergyCrystal", "phi",        "sigmaEta2012", "sigmaIEtaIEta2012",
  };

  for (string branch : defaultBranches) {
    histogramsHandler->Fill("goodPhoton_" + prefix + branch, photon->GetAs<float>(branch));
  }
}

void LbLHistogramsFiller::FillEGammaHistograms(const shared_ptr<Event> event) {
  auto eGammaObjects = event->GetCollection("egamma");
  auto photons = event->GetCollection("goodPhoton");
  auto photon = asPhoton(photons->at(0));
  float et = photon->Get("et");
  if (asLbLEvent(event)->IsData() && et > dataBlinding["max_et"]) return;

  shared_ptr<PhysicsObject> closestEgamma = nullptr;
  float minDeltaR = 9999;

  if (!closestEgamma) {
    warn() << "No EGamma objects found in the event." << endl;
    return;
  }

  for (int i = 0; i < eGammaObjects->size(); i++) {
    auto eGamma = eGammaObjects->at(i);

    float eta1 = eGamma->Get("eta");
    float phi1 = eGamma->Get("phi");

    // check deltaR from the good photon
    float deltaR = sqrt(pow(eta1 - photon->GetAs<float>("eta"), 2) + pow(TVector2::Phi_mpi_pi(phi1 - photon->GetAs<float>("phi")), 2));

    if (deltaR < minDeltaR) {
      minDeltaR = deltaR;
      closestEgamma = eGamma;
    }

    for (int j = i + 1; j < eGammaObjects->size(); j++) {
      auto eGamma2 = eGammaObjects->at(j);

      float eta2 = eGamma2->Get("eta");
      float phi2 = eGamma2->Get("phi");

      histogramsHandler->Fill("monophoton_egamma_deltaR", sqrt(pow(eta1 - eta2, 2) + pow(TVector2::Phi_mpi_pi(phi1 - phi2), 2)));
      histogramsHandler->Fill("monophoton_egamma_deltaEta", fabs(eta1 - eta2));
      histogramsHandler->Fill("monophoton_egamma_deltaPhi", fabs(TVector2::Phi_mpi_pi(phi1 - phi2)));
    }
  }

  if (!closestEgamma) return;

  histogramsHandler->Fill("egamma_et_vs_goodPhoton_et", closestEgamma->Get("et"), et);
}

void LbLHistogramsFiller::SaveHighEtPhotonsInfo(const shared_ptr<Event> event) {
  auto photons = event->GetCollection("goodPhoton");
  auto photon = asPhoton(photons->at(0));

  // histogramsHandler->Fill("goodPhoton_eta_vs_phi_gt50GeV", (float)photon->Get("eta"), (float)photon->Get("phi"), et);
  // open a text file to write the photon information
  ofstream photonFile("/afs/desy.de/user/j/jniedzie/tea_lbl/photon_info" + to_string(event->GetAs<int>("eventNumber")) + ".txt");
  if (photonFile.is_open()) {
    photonFile << "Event: " << event->GetAs<int>("eventNumber") << ", Lumi section: " << event->GetAs<int>("lumiSection")
               << ", Run: " << event->GetAs<int>("runNumber") << endl;
    // print event vertex
    auto vertex = event->GetCollection("vertex")->at(0);
    photonFile << "\n\nVertex information:" << endl;
    if (vertex) {
      photonFile << "vertex_x: " << vertex->GetAs<float>("x") << endl;
      photonFile << "vertex_y: " << vertex->GetAs<float>("y") << endl;
      photonFile << "vertex_z: " << vertex->GetAs<float>("z") << endl;
    } else {
      photonFile << "No vertex information available." << endl;
    }

    photonFile << "\n\nPhoton information:" << endl;
    photonFile << "photon_energy: " << photon->GetAs<float>("energy") << endl;
    photonFile << "photon_et: " << photon->GetAs<float>("et") << endl;
    photonFile << "photon_eta: " << photon->GetAs<float>("eta") << endl;
    photonFile << "photon_phi: " << photon->GetAs<float>("phi") << endl;
    photonFile << "photon_hOverE: " << photon->GetAs<float>("hOverE") << endl;
    photonFile << "photon_sigmaEta2012: " << photon->GetAs<float>("sigmaEta2012") << endl;
    photonFile << "photon_sigmaIEtaIEta2012: " << photon->GetAs<float>("sigmaIEtaIEta2012") << endl;
    photonFile << "photon_maxEnergyCrystal: " << photon->GetAs<float>("maxEnergyCrystal") << endl;
    photonFile << "photon_energyTop: " << photon->GetAs<float>("energyTop") << endl;
    photonFile << "photon_energyBottom: " << photon->GetAs<float>("energyBottom") << endl;
    photonFile << "photon_energyLeft: " << photon->GetAs<float>("energyLeft") << endl;
    photonFile << "photon_energyRight: " << photon->GetAs<float>("energyRight") << endl;
    photonFile << "photon_hasConversionTracks: " << photon->GetAs<int>("hasConversionTracks") << endl;
    photonFile << "photon_seedTime: " << photon->GetAs<float>("seedTime") << endl;
    photonFile << "photon_SCEnergy: " << photon->GetAs<float>("SCEnergy") << endl;
    photonFile << "photon_SCEt: " << photon->GetAs<float>("SCEt") << endl;
    photonFile << "photon_SCEta: " << photon->GetAs<float>("SCEta") << endl;
    photonFile << "photon_SCPhi: " << photon->GetAs<float>("SCPhi") << endl;
    photonFile << "photon_SCEtaWidth: " << photon->GetAs<float>("SCEtaWidth") << endl;
    photonFile << "photon_SCPhiWidth: " << photon->GetAs<float>("SCPhiWidth") << endl;
    photonFile << "photon_ecalClusterIsoR2: " << photon->GetAs<float>("ecalClusterIsoR2") << endl;
    photonFile << "photon_ecalClusterIsoR3: " << photon->GetAs<float>("ecalClusterIsoR3") << endl;
    photonFile << "photon_ecalClusterIsoR4: " << photon->GetAs<float>("ecalClusterIsoR4") << endl;
    photonFile << "photon_ecalClusterIsoR5: " << photon->GetAs<float>("ecalClusterIsoR5") << endl;
    photonFile << "photon_hcalRechitIsoR1: " << photon->GetAs<float>("hcalRechitIsoR1") << endl;
    photonFile << "photon_hcalRechitIsoR2: " << photon->GetAs<float>("hcalRechitIsoR2") << endl;
    photonFile << "photon_hcalRechitIsoR3: " << photon->GetAs<float>("hcalRechitIsoR3") << endl;
    photonFile << "photon_hcalRechitIsoR4: " << photon->GetAs<float>("hcalRechitIsoR4") << endl;
    photonFile << "photon_hcalRechitIsoR5: " << photon->GetAs<float>("hcalRechitIsoR5") << endl;

    // find calo tower closest to the photon:
    auto caloTowers = event->GetCollection("CaloTower");
    shared_ptr<CaloTower> closestCaloTower = nullptr;

    float minDeltaR = 9999;
    for (auto caloTowerObj : *caloTowers) {
      auto caloTower = asCaloTower(caloTowerObj);
      float towerEta = caloTower->Get("eta");
      float towerPhi = caloTower->Get("phi");
      float photonEta = photon->Get("eta");
      float photonPhi = photon->Get("phi");

      float deltaR = sqrt(pow(towerEta - photonEta, 2) + pow(TVector2::Phi_mpi_pi(towerPhi - photonPhi), 2));

      if (deltaR < minDeltaR) {
        minDeltaR = deltaR;
        closestCaloTower = caloTower;
      }
    }
    photonFile << "\n\nClosest Calo Tower information:" << endl;
    if (closestCaloTower) {
      photonFile << "closestCaloTower_energy: " << closestCaloTower->GetAs<float>("energy") << endl;
      photonFile << "closestCaloTower_et: " << closestCaloTower->GetAs<float>("et") << endl;
      photonFile << "closestCaloTower_eta: " << closestCaloTower->GetAs<float>("eta") << endl;
      photonFile << "closestCaloTower_phi: " << closestCaloTower->GetAs<float>("phi") << endl;
      photonFile << "closestCaloTower_hadE: " << closestCaloTower->GetAs<float>("hadE") << endl;
      photonFile << "closestCaloTower_emE: " << closestCaloTower->GetAs<float>("emE") << endl;
    } else {
      photonFile << "No closest Calo Tower found." << endl;
    }

    // print all photons (also the ones that don't pass good photon selection)
    auto allPhotons = event->GetCollection("photon");

    photonFile << "\n\nAll photons information:" << endl;
    for (const auto& photonObj : *allPhotons) {
      auto pho = asPhoton(photonObj);
      photonFile << "\nphoton_et: " << pho->GetAs<float>("et") << endl;
      photonFile << "photon_eta: " << pho->GetAs<float>("eta") << endl;
      photonFile << "photon_phi: " << pho->GetAs<float>("phi") << endl;
      photonFile << "photon_energy: " << pho->GetAs<float>("energy") << endl;
      photonFile << "photon_hOverE: " << pho->GetAs<float>("hOverE") << endl;
      photonFile << "photon_sigmaEta2012: " << pho->GetAs<float>("sigmaEta2012") << endl;
      photonFile << "photon_sigmaIEtaIEta2012: " << pho->GetAs<float>("sigmaIEtaIEta2012") << endl;
      photonFile << "photon_maxEnergyCrystal: " << pho->GetAs<float>("maxEnergyCrystal") << endl;
      photonFile << "photon_energyTop: " << pho->GetAs<float>("energyTop") << endl;
      photonFile << "photon_energyBottom: " << pho->GetAs<float>("energyBottom") << endl;
      photonFile << "photon_energyLeft: " << pho->GetAs<float>("energyLeft") << endl;
      photonFile << "photon_energyRight: " << pho->GetAs<float>("energyRight") << endl;
    }

    // print tracks info
    auto tracks = event->GetCollection("track");
    photonFile << "\n\nTracks information:" << endl;
    if (tracks) {
      for (const auto& trackObj : *tracks) {
        auto track = asTrack(trackObj);
        photonFile << "\ntrack_pt: " << track->GetAs<float>("pt") << endl;
        photonFile << "track_eta: " << track->GetAs<float>("eta") << endl;
        photonFile << "track_phi: " << track->GetAs<float>("phi") << endl;
        photonFile << "track_charge: " << track->GetAs<int>("charge") << endl;
        photonFile << "track_dxy: " << track->GetAs<float>("dxy") << endl;
        photonFile << "track_dz: " << track->GetAs<float>("dz") << endl;
      }
    } else {
      photonFile << "No tracks information available." << endl;
    }

    // print electrons info
    auto electrons = event->GetCollection("electron");
    photonFile << "\n\nElectrons information:" << endl;
    if (electrons) {
      for (const auto& electronObj : *electrons) {
        auto electron = asElectron(electronObj);
        photonFile << "\nelectron_pt: " << electron->GetAs<float>("pt") << endl;
        photonFile << "electron_eta: " << electron->GetAs<float>("eta") << endl;
        photonFile << "electron_phi: " << electron->GetAs<float>("phi") << endl;
        photonFile << "electron_charge: " << electron->GetAs<int>("charge") << endl;
      }
    } else {
      photonFile << "No electrons information available." << endl;
    }

    // print muon info
    auto muons = event->GetCollection("muon");
    photonFile << "\n\nMuons information:" << endl;
    if (muons) {
      for (const auto& muonObj : *muons) {
        auto muon = asMuon(muonObj);
        photonFile << "\nmuon_pt: " << muon->GetAs<float>("pt") << endl;
        photonFile << "muon_eta: " << muon->GetAs<float>("eta") << endl;
        photonFile << "muon_phi: " << muon->GetAs<float>("phi") << endl;
        photonFile << "muon_charge: " << muon->GetAs<int>("charge") << endl;
      }
    } else {
      photonFile << "No muons information available." << endl;
    }

    // print EGamma objects
    auto eGammaObjects = event->GetCollection("egamma");
    photonFile << "\n\nEGamma objects information:" << endl;
    if (eGammaObjects) {
      for (const auto& eGammaObj : *eGammaObjects) {
        photonFile << "\neGamma_et: " << eGammaObj->GetAs<float>("et") << endl;
        photonFile << "eGamma_eta: " << eGammaObj->GetAs<float>("eta") << endl;
        photonFile << "eGamma_phi: " << eGammaObj->GetAs<float>("phi") << endl;
      }

      for (int i = 0; i < eGammaObjects->size(); i++) {
        auto eGamma = eGammaObjects->at(i);

        float eta1 = eGamma->Get("eta");
        float phi1 = eGamma->Get("phi");

        for (int j = i + 1; j < eGammaObjects->size(); j++) {
          auto eGamma2 = eGammaObjects->at(j);

          float eta2 = eGamma2->Get("eta");
          float phi2 = eGamma2->Get("phi");

          histogramsHandler->Fill("monophoton_egamma_deltaR_gt50GeV",
                                  sqrt(pow(eta1 - eta2, 2) + pow(TVector2::Phi_mpi_pi(phi1 - phi2), 2)));
          histogramsHandler->Fill("monophoton_egamma_deltaEta_gt50GeV", fabs(eta1 - eta2));
          histogramsHandler->Fill("monophoton_egamma_deltaPhi_gt50GeV", fabs(TVector2::Phi_mpi_pi(phi1 - phi2)));
        }
      }

    } else {
      photonFile << "No EGamma objects information available." << endl;
    }

    photonFile.close();
  } else {
    warn() << "Unable to open photon_info.txt for writing." << endl;
  }
}

void LbLHistogramsFiller::FillGenLevelHistograms(const shared_ptr<Event> event) {
  auto photons = event->GetCollection("genPhoton");
  auto electrons = event->GetCollection("genElectron");

  for (auto physObject : *photons) {
    auto photon = asPhoton(physObject)->GetFourMomentum();
    if (fabs(photon.Eta()) > 5.2) continue;

    bool overlapsWithElectron = false;

    for (auto electron : *electrons) {
      if (photon.DeltaR(asElectron(electron)->GetFourMomentum()) < 0.1) {
        overlapsWithElectron = true;
        break;
      }
    }
    if (overlapsWithElectron) continue;

    histogramsHandler->Fill("genPhoton_et", photon.Pt());
    histogramsHandler->Fill("genPhoton_energy", photon.E());
  }
}

void LbLHistogramsFiller::FillEventLevelHistograms(const shared_ptr<Event> event) {
  try {
    auto zdcEnergies = event->GetCollection("ZDC");

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
    histogramsHandler->Fill("event_ZDCenergyPlus", totalEnergyPlus);
    histogramsHandler->Fill("event_ZDCenergyMinus", totalEnergyMinus);
    histogramsHandler->Fill("event_ZDCenergyPlusLogX", TMath::Log10(totalEnergyPlus));
    histogramsHandler->Fill("event_ZDCenergyMinusLogX", TMath::Log10(totalEnergyMinus));
  } catch (const Exception& e) {
    warn() << "Cannot fill ZDC histograms, since ZDC collection was not found." << endl;
  }
}

void LbLHistogramsFiller::Fill(const shared_ptr<Event> event) {
  auto photons = event->GetCollection("goodPhoton");

  if (photons->size() != 1) {
    fatal() << "Number of good photons != 1. This should never happen." << endl;
    exit(13);
  }

  FillMonoPhotonHistograms(event);
  FillEGammaHistograms(event);
  FillEventLevelHistograms(event);
  FillGenLevelHistograms(event);

  // auto photon = asPhoton(photons->at(0));
  // float et = photon->Get("et");
  // if (et > 50) SaveHighEtPhotonsInfo(event);
}
