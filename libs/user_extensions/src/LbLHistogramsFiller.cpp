#include "LbLHistogramsFiller.hpp"

#include "ConfigManager.hpp"
#include "ExtensionsHelpers.hpp"
#include "UserExtensionsHelpers.hpp"

using namespace std;

LbLHistogramsFiller::LbLHistogramsFiller(shared_ptr<HistogramsHandler> histogramsHandler_) : histogramsHandler(histogramsHandler_) {
  // Create a config manager
  auto &config = ConfigManager::GetInstance();

  config.GetMap("caloEtaEdges", caloEtaEdges);

  try {
    config.GetMap("dataBlinding", dataBlinding);
  } catch (const Exception &e) {
    warn() << "No data blinding parameters found. Will not apply any data blinding." << endl;
    dataBlinding["max_et"] = 9999;
  }

  eventProcessor = make_unique<EventProcessor>();

  acoplanarityFunction = new TF1("fit_fun", "[0] * exp([1]*x) + [2] * exp([3]*x)", 0, 0.1);
  vector<float> params = {102020, -798.6, 693.552, -112.731};
  for (int i = 0; i < 4; i++) {
    acoplanarityFunction->SetParameter(i, params[i]);
  }

  // auto file = TFile::Open("/afs/desy.de/user/j/jniedzie/tea_lbl/merged_afterTrigger_histograms_default.root");
  // genDielectronDeltaPt = (TH1D *)file->Get("genDielectron_deltaPt");
}

LbLHistogramsFiller::~LbLHistogramsFiller() {}

void LbLHistogramsFiller::FillCaloHistograms(const shared_ptr<Event> event) {
  auto towers = event->GetCollection("CaloTower");

  float maxEnergyHE = -1, maxEnergyHFplus = -1, maxEnergyHFminus = -1;
  shared_ptr<CaloTower> leadingTowerHE, leadingTowerHFplus, leadingTowerHFminus;

  // info() << "Event number: " << (int)event->Get("eventNumber") << "\trun: " << (int)event->Get("runNumber") << endl;

  for (auto physicsObject : *towers) {
    auto tower = asCaloTower(physicsObject);

    if (tower->IsDead()) continue;
    if (tower->IsInHadronicCrack()) continue;
    if (tower->IsInHEM()) continue;

    float eta = (float)physicsObject->Get("eta");
    float energy = physicsObject->Get("energy");

    // HE
    if ((fabs(eta) > caloEtaEdges["minHE"] && fabs(eta) < caloEtaEdges["maxHE"])) {
      if (energy > maxEnergyHE) {
        maxEnergyHE = energy;
        leadingTowerHE = tower;
      }
    }

    // HF+
    if (eta > caloEtaEdges["minHF"] && eta < caloEtaEdges["maxHF"]) {
      if (energy > maxEnergyHFplus) {
        maxEnergyHFplus = energy;
        leadingTowerHFplus = tower;
      }
    }

    // HF-
    if (eta > -caloEtaEdges["maxHF"] && eta < -caloEtaEdges["minHF"]) {
      if (energy > maxEnergyHFminus) {
        maxEnergyHFminus = energy;
        leadingTowerHFminus = tower;
      }
    }
  }

  if (maxEnergyHE > 0) {
    histogramsHandler->Fill("goodCaloTowerHE_energyHad", leadingTowerHE->GetAs<float>("hadE"));
    histogramsHandler->Fill("goodCaloTowerHE_energyTransverse", leadingTowerHE->GetAs<float>("et"));
    histogramsHandler->Fill("goodCaloTowerHE_energy", leadingTowerHE->GetAs<float>("energy"));
  }

  if (maxEnergyHFplus > 0) {
    // info() << (int)event->Get("eventNumber") << "\t" << leadingTowerHFplus->GetAs<float>("eta") << "\t" <<
    // leadingTowerHFplus->GetAs<float>("energy") << endl;
    histogramsHandler->Fill("goodCaloTowerHFplus_energyHad", leadingTowerHFplus->GetAs<float>("hadE"));
    histogramsHandler->Fill("goodCaloTowerHFplus_energyTransverse", leadingTowerHFplus->GetAs<float>("et"));
    histogramsHandler->Fill("goodCaloTowerHFplus_energy", leadingTowerHFplus->GetAs<float>("energy"));
  }

  if (maxEnergyHFminus > 0) {
    // info() << "Filling HF- tower with eta: " << leadingTowerHFminus->GetAs<float>("eta") << "\tenergy: " <<
    // leadingTowerHFminus->GetAs<float>("energy") << endl;
    histogramsHandler->Fill("goodCaloTowerHFminus_energyHad", leadingTowerHFminus->GetAs<float>("hadE"));
    histogramsHandler->Fill("goodCaloTowerHFminus_energyTransverse", leadingTowerHFminus->GetAs<float>("et"));
    histogramsHandler->Fill("goodCaloTowerHFminus_energy", leadingTowerHFminus->GetAs<float>("energy"));
  }
}

void LbLHistogramsFiller::FillThreePhotonHistograms(const shared_ptr<Event> event) {
  auto photons = event->GetCollection("goodPhoton");
  if (photons->size() != 3) return;

  auto photon0 = photons->at(0);
  auto photon1 = photons->at(1);
  auto photon2 = photons->at(2);

  TLorentzVector photon0vec, photon1vec, photon2vec;
  photon0vec.SetPtEtaPhiM(photon0->Get("et"), photon0->Get("eta"), photon0->Get("phi"), 0);
  photon1vec.SetPtEtaPhiM(photon1->Get("et"), photon1->Get("eta"), photon1->Get("phi"), 0);
  photon2vec.SetPtEtaPhiM(photon2->Get("et"), photon2->Get("eta"), photon2->Get("phi"), 0);

  // diphoton plots
  histogramsHandler->Fill("diphotonThreePhoton_mass", (photon0vec + photon1vec).M());
  histogramsHandler->Fill("diphotonThreePhoton_mass", (photon0vec + photon2vec).M());
  histogramsHandler->Fill("diphotonThreePhoton_mass", (photon1vec + photon2vec).M());

  float deltaR01 = photon0vec.DeltaR(photon1vec);
  float deltaR02 = photon0vec.DeltaR(photon2vec);
  float deltaR12 = photon1vec.DeltaR(photon2vec);

  if (deltaR01 < deltaR02 && deltaR01 < deltaR12) {
    histogramsHandler->Fill("diphotonThreePhoton_massClosestDeltaR", (photon0vec + photon1vec).M());
  } else if (deltaR02 < deltaR12) {
    histogramsHandler->Fill("diphotonThreePhoton_massClosestDeltaR", (photon0vec + photon2vec).M());
  } else {
    histogramsHandler->Fill("diphotonThreePhoton_massClosestDeltaR", (photon1vec + photon2vec).M());
  }

  // triphoton plots
  auto triphoton = photon0vec + photon1vec + photon2vec;
  histogramsHandler->Fill("triphoton_mass", triphoton.M());
  histogramsHandler->Fill("triphoton_pt", triphoton.Pt());
  histogramsHandler->Fill("triphoton_eta", triphoton.Eta());
  histogramsHandler->Fill("triphoton_phi", triphoton.Phi());
}

void LbLHistogramsFiller::FillMonoPhotonHistograms(const shared_ptr<Event> event) {
  auto photons = event->GetCollection("goodPhoton");
  auto photon = asPhoton(photons->at(0));
  float et = photon->Get("et");

  // on the fly cuts (to be moved to selections)
  float horizontalOverCentralEnergy = photon->GetHorizontalOverCentralEnergy();
  float verticalOverCentralEnergy = photon->GetVerticalOverCentralEnergy();
  float horizontalImbalance = photon->GetHorizontalImbalance();
  float verticalImbalance = photon->GetVerticalImbalance();

  // if (horizontalOverCentralEnergy < 0.005 || verticalOverCentralEnergy < 0.005) return;
  if (horizontalImbalance < -1.0 || verticalImbalance < -1.0) return;

  if (asLbLEvent(event)->IsData() && et > dataBlinding["max_et"]) return;

  histogramsHandler->Fill("goodPhoton_absEta_vs_et", fabs((float)photon->Get("eta")), et);
  histogramsHandler->Fill("goodPhoton_eta_vs_et", (float)photon->Get("eta"), et);
  histogramsHandler->Fill("goodPhoton_eta_vs_phi", (float)photon->Get("eta"), (float)photon->Get("phi"));
  // histogramsHandler->Fill("goodPhoton_eta_vs_phi_vs_et", (float)photon->Get("eta"), (float)photon->Get("phi"), et);

  histogramsHandler->Fill("goodPhoton_et", et);
  histogramsHandler->Fill("goodPhoton_logEt", TMath::Log10(et));
  histogramsHandler->Fill("goodPhoton_seedTime", photon->GetAs<float>("seedTime"));

  float energyTop = photon->GetAs<float>("energyTop");
  float energyBottom = photon->GetAs<float>("energyBottom");
  float energyLeft = photon->GetAs<float>("energyLeft");
  float energyRight = photon->GetAs<float>("energyRight");
  float energyCentral = photon->GetAs<float>("maxEnergyCrystal");

  float minEnergy = min({energyTop, energyBottom, energyLeft, energyRight});

  histogramsHandler->Fill("goodPhoton_topOverCentral", energyTop/energyCentral);
  histogramsHandler->Fill("goodPhoton_bottomOverCentral", energyBottom/energyCentral);
  histogramsHandler->Fill("goodPhoton_leftOverCentral", energyLeft/energyCentral);
  histogramsHandler->Fill("goodPhoton_rightOverCentral", energyRight/energyCentral);
  histogramsHandler->Fill("goodPhoton_minOverCentral", minEnergy/energyCentral);
  histogramsHandler->Fill("goodPhoton_verticalOverCentral", photon->GetVerticalOverCentralEnergy());
  histogramsHandler->Fill("goodPhoton_horizontalOverCentral", photon->GetHorizontalOverCentralEnergy());
  histogramsHandler->Fill("goodPhoton_horizontalImbalance", horizontalImbalance);
  histogramsHandler->Fill("goodPhoton_verticalImbalance", verticalImbalance);

  auto eGammaObjects = event->GetCollection("egamma");

  shared_ptr<PhysicsObject> closestEgamma = nullptr;
  float minDeltaR = 9999;

  if (eGammaObjects) {
    

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
  }

  if (closestEgamma) {
    histogramsHandler->Fill("egamma_et_vs_goodPhoton_et", closestEgamma->Get("et"), et);
  } else {
    warn() << "No EGamma objects found in the event." << endl;
  }

  if (et > 50) {
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
      for (const auto &photonObj : *allPhotons) {
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
        for (const auto &trackObj : *tracks) {
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
        for (const auto &electronObj : *electrons) {
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
        for (const auto &muonObj : *muons) {
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
        for (const auto &eGammaObj : *eGammaObjects) {
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
}

void LbLHistogramsFiller::FillPhotonHistograms(const shared_ptr<Event> event) {
  auto photons = event->GetCollection("goodPhoton");

  if (photons->size() == 1) {
    FillMonoPhotonHistograms(event);
    return;
  }

  if (photons->size() != 2) return;
  auto photon1 = photons->at(0);
  auto photon2 = photons->at(1);
  TLorentzVector photon1vec, photon2vec;
  photon1vec.SetPtEtaPhiM(photon1->Get("et"), photon1->Get("eta"), photon1->Get("phi"), 0);
  photon2vec.SetPtEtaPhiM(photon2->Get("et"), photon2->Get("eta"), photon2->Get("phi"), 0);

  auto diphoton = photon1vec + photon2vec;

  double deltaPhi = photon1vec.DeltaPhi(photon2vec);
  double acoplanarity = 1 - (fabs(deltaPhi) / TMath::Pi());

  histogramsHandler->Fill("diphoton_pt", diphoton.Pt());
  histogramsHandler->Fill("diphoton_massNew", diphoton.M());
  histogramsHandler->Fill("diphoton_mass100", diphoton.M());
  histogramsHandler->Fill("diphoton_mass200", diphoton.M());
  histogramsHandler->Fill("diphoton_mass", diphoton.M());
  histogramsHandler->Fill("diphoton_rapidity", diphoton.Rapidity());
  histogramsHandler->Fill("diphoton_acoplanarity20", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity25", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity30", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity32", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity33", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity34", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity35", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity40", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity45", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity50", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity55", acoplanarity);
  histogramsHandler->Fill("diphoton_acoplanarity60", acoplanarity);

  auto tracks = event->GetCollection("track");
  auto electrons = event->GetCollection("goodElectron");
  shared_ptr<Track> track1 = nullptr;

  for (auto physicsObject : *tracks) {
    auto track = asTrack(physicsObject);
    if (track->OverlapsWithOtherObjects(photons) && !track->OverlapsWithOtherObjects(electrons)) {
      if (track1 == nullptr) {
        track1 = track;
        break;
      }
    }
  }
  if (track1) {
    float trackPt = track1->Get("pt");
    if (trackPt < 0.3) histogramsHandler->Fill("diphoton_acoplanarityTrack0to0p3", acoplanarity);
    if (trackPt > 0.3 && trackPt < 0.65) histogramsHandler->Fill("diphoton_acoplanarityTrack0p3to0p65", acoplanarity);
    if (trackPt > 0.6 && trackPt < 2.0) histogramsHandler->Fill("diphoton_acoplanarityTrack0p65to2p0", acoplanarity);
    if (trackPt > 2.0) histogramsHandler->Fill("diphoton_acoplanarityTrack2p0toInf", acoplanarity);
  }
  histogramsHandler->Fill("diphoton_acoplanarity1040", acoplanarity);

  auto hist = histogramsHandler->GetHistogram1D({"diphoton_acoplanarity1", ""});
  float binWidth1 = hist->GetBinWidth(1);
  float binWidth = hist->GetBinWidth(hist->GetXaxis()->FindFixBin(acoplanarity));
  // histogramsHandler->Fill("diphoton_acoplanarity1", acoplanarity, binWidth1/binWidth);
  // histogramsHandler->Fill("diphoton_acoplanarity1", acoplanarity, 1/binWidth);
  histogramsHandler->Fill("diphoton_acoplanarity1", acoplanarity);

  hist = histogramsHandler->GetHistogram1D({"diphoton_acoplanarity2", ""});
  binWidth1 = hist->GetBinWidth(1);
  binWidth = hist->GetBinWidth(hist->GetXaxis()->FindFixBin(acoplanarity));
  // histogramsHandler->Fill("diphoton_acoplanarity2", acoplanarity, binWidth1/binWidth);
  histogramsHandler->Fill("diphoton_acoplanarity2", acoplanarity);

  hist = histogramsHandler->GetHistogram1D({"diphoton_acoplanarity3", ""});
  binWidth1 = hist->GetBinWidth(1);
  binWidth = hist->GetBinWidth(hist->GetXaxis()->FindFixBin(acoplanarity));
  // histogramsHandler->Fill("diphoton_acoplanarity3", acoplanarity, binWidth1/binWidth);
  histogramsHandler->Fill("diphoton_acoplanarity3", acoplanarity);

  hist = histogramsHandler->GetHistogram1D({"diphoton_acoplanarity4", ""});
  binWidth1 = hist->GetBinWidth(1);
  binWidth = hist->GetBinWidth(hist->GetXaxis()->FindFixBin(acoplanarity));
  // histogramsHandler->Fill("diphoton_acoplanarity4", acoplanarity, binWidth1/binWidth);
  histogramsHandler->Fill("diphoton_acoplanarity4", acoplanarity);

  hist = histogramsHandler->GetHistogram1D({"diphoton_acoplanarity5", ""});
  binWidth1 = hist->GetBinWidth(1);
  binWidth = hist->GetBinWidth(hist->GetXaxis()->FindFixBin(acoplanarity));
  // histogramsHandler->Fill("diphoton_acoplanarity4", acoplanarity, binWidth1/binWidth);
  histogramsHandler->Fill("diphoton_acoplanarity5", acoplanarity);

  histogramsHandler->Fill("diphoton_seedTime", photon1->Get("seedTime"), photon2->Get("seedTime"));

  if (acoplanarity < 0.01) {
    histogramsHandler->Fill("diphotonSR_pt", diphoton.Pt());
    histogramsHandler->Fill("diphotonSR_mass", diphoton.M());
    histogramsHandler->Fill("diphotonSR_massNew", diphoton.M());
    histogramsHandler->Fill("diphotonSR_mass100", diphoton.M());
    histogramsHandler->Fill("diphotonSR_mass200", diphoton.M());
    histogramsHandler->Fill("diphotonSR_rapidity", diphoton.Rapidity());

    histogramsHandler->Fill("goodPhotonSR_et", photon1->Get("et"));
    histogramsHandler->Fill("goodPhotonSR_eta", photon1->Get("eta"));
    histogramsHandler->Fill("goodPhotonSR_phi", photon1->Get("phi"));
    histogramsHandler->Fill("goodPhotonSR_et", photon2->Get("et"));
    histogramsHandler->Fill("goodPhotonSR_eta", photon2->Get("eta"));
    histogramsHandler->Fill("goodPhotonSR_phi", photon2->Get("phi"));

    histogramsHandler->Fill("diphoton_seedTimeSR", photon1->Get("seedTime"), photon2->Get("seedTime"));
    histogramsHandler->Fill("unfoldingPhoton_pt", diphoton.Pt());
    histogramsHandler->Fill("unfoldingPhoton_mass", diphoton.M());
    histogramsHandler->Fill("unfoldingPhoton_absRap", fabs(diphoton.Rapidity()));
    histogramsHandler->Fill("unfoldingPhoton_absRap3", fabs(diphoton.Rapidity()));
    histogramsHandler->Fill("unfoldingPhoton_rap3", diphoton.Rapidity());
    histogramsHandler->Fill("unfoldingPhoton_rap4", diphoton.Rapidity());
  }
}

void LbLHistogramsFiller::FillGenLevelHistograms(const shared_ptr<Event> event) {
  float leadingPhotonEnergy = 99999;
  float leadingPhotonEnergyBarrel = 99999;
  float leadingPhotonEnergyBarrelEndcap = 99999;

  float leadingPhotonEt = 0;
  float leadingPhotonEtBarrel = 0;
  float leadingPhotonEtBarrelEndcap = 0;

  auto photons = event->GetCollection("genPhoton");
  auto electrons = event->GetCollection("genElectron");

  for (auto physObject : *photons) {
    auto photon = asPhoton(physObject)->GetFourMomentum();

    if (fabs(photon.Eta()) > 5.2) continue;
    if (GetDielectronAcoplanarity(event) > 0.01) continue;

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

    if (photon.Pt() > leadingPhotonEt) {
      leadingPhotonEt = photon.Pt();
      leadingPhotonEnergy = photon.E();
    }
    if (photon.Pt() > leadingPhotonEtBarrel && fabs(photon.Eta()) < 1.4442) {
      leadingPhotonEtBarrel = photon.Pt();
      leadingPhotonEnergyBarrel = photon.E();
    }

    if (photon.Pt() > leadingPhotonEtBarrelEndcap && fabs(photon.Eta()) < 3.0) {
      leadingPhotonEtBarrelEndcap = photon.Pt();
      leadingPhotonEnergyBarrelEndcap = photon.E();
    }
  }
  if (GetDielectronAcoplanarity(event) < 0.01) {
    histogramsHandler->Fill("leadingGenPhoton_energy", leadingPhotonEnergy);
    histogramsHandler->Fill("leadingGenPhotonBarrel_energy", leadingPhotonEnergyBarrel);
    histogramsHandler->Fill("leadingGenPhotonBarrelEndcap_energy", leadingPhotonEnergyBarrelEndcap);
  }

  if (electrons->size() == 2) {
    auto electron_1 = asElectron(electrons->at(0));
    auto electron_2 = asElectron(electrons->at(1));

    float deltaPhi = GetPhiModulation(electron_1, electron_2);
    float deltaPt = electron_1->GetFourMomentum().Pt() - electron_2->GetFourMomentum().Pt();

    histogramsHandler->Fill("genDielectron_deltaPhi", deltaPhi);
    histogramsHandler->Fill("genDielectron_deltaPt", deltaPt);
    histogramsHandler->Fill("genDielectron_pt", genDielectron.Pt());

    if (GetDielectronAcoplanarity(event) < 0.01) {
      histogramsHandler->Fill("genDielectronSR_deltaPhi", deltaPhi);
      histogramsHandler->Fill("genDielectronSR_deltaPt", deltaPt);
    }
  }
}

float LbLHistogramsFiller::GetDielectronAcoplanarity(const shared_ptr<Event> event) {
  auto electrons = event->GetCollection("goodElectron");
  if (electrons->size() != 2) return 1000;
  return GetDielectronAcoplanarity(asElectron(electrons->at(0)), asElectron(electrons->at(1)));
}

float LbLHistogramsFiller::GetDielectronAcoplanarity(const shared_ptr<Electron> &electron1, const shared_ptr<Electron> &electron2) {
  auto electron1vec = electron1->GetFourMomentum();
  auto electron2vec = electron2->GetFourMomentum();
  double deltaPhi = electron1vec.DeltaPhi(electron2vec);
  double acoplanarity = 1 - (fabs(deltaPhi) / TMath::Pi());
  return acoplanarity;
}

float LbLHistogramsFiller::GetPhiModulation(const shared_ptr<Electron> &electron1, const shared_ptr<Electron> &electron2) {
  TLorentzVector electron, positron;
  if (electron1->GetCharge() > 0) {
    positron = electron1->GetFourMomentum();
    electron = electron2->GetFourMomentum();
  } else {
    positron = electron2->GetFourMomentum();
    electron = electron1->GetFourMomentum();
  }

  // add a little random gaussian spread to pt of electron and positron
  // electron.SetPtEtaPhiM(electron.Pt() + gRandom->Gaus(0, 0.007*electron.Pt()), electron.Eta(), electron.Phi(), electron.M());
  // positron.SetPtEtaPhiM(positron.Pt() + gRandom->Gaus(0, 0.007*positron.Pt()), positron.Eta(), positron.Phi(), positron.M());

  // draw random sign (+1 or -1)
  // int sign = gRandom->Rndm() > 0.5 ? 1 : -1;
  // positron.SetPtEtaPhiM(electron.Pt() + sign * 0.01712, positron.Eta(), positron.Phi(), positron.M());

  // positron.SetPtEtaPhiM(electron.Pt() + genDielectronDeltaPt->GetRandom(), positron.Eta(), positron.Phi(), positron.M());

  TLorentzVector dielectron = electron + positron;

  float deltaPhi = fabs(dielectron.Phi() - electron.Phi());
  if (deltaPhi > TMath::Pi()) deltaPhi = 2. * TMath::Pi() - deltaPhi;
  return deltaPhi;
}

void LbLHistogramsFiller::FillElectronHistograms(const shared_ptr<Event> event) {
  auto electrons = event->GetCollection("goodElectron");

  if (electrons->size() != 2) return;

  auto electron1 = asElectron(electrons->at(0));
  auto electron2 = asElectron(electrons->at(1));
  float acoplanarity = GetDielectronAcoplanarity(electron1, electron2);
  float deltaPhi = GetPhiModulation(electron1, electron2);
  float deltaPt = fabs(electron1->GetFourMomentum().Pt() - electron2->GetFourMomentum().Pt());

  histogramsHandler->Fill("dielectron_pt", dielectron.Pt());
  histogramsHandler->Fill("dielectron_mass", dielectron.M());
  histogramsHandler->Fill("dielectron_rapidity", dielectron.Rapidity());
  histogramsHandler->Fill("dielectron_acoplanarity", acoplanarity);

  if (acoplanarity < 0.01) {
    histogramsHandler->Fill("dielectronacoLt0p01_mass", dielectron.M());
  }

  if (electron1->GetCharge() > 0 && electron2->GetCharge() > 0) {
    histogramsHandler->Fill("dielectronPP_mass", dielectron.M());
    if (acoplanarity < 0.01) {
      histogramsHandler->Fill("dielectronPPacoLt0p01_mass", dielectron.M());
    }
  } else if (electron1->GetCharge() < 0 && electron2->GetCharge() < 0) {
    histogramsHandler->Fill("dielectronMM_mass", dielectron.M());
    if (acoplanarity < 0.01) {
      histogramsHandler->Fill("dielectronMMacoLt0p01_mass", dielectron.M());
    }
  }

  histogramsHandler->Fill("dielectron_deltaPt", deltaPt);
  histogramsHandler->Fill("electrons_deltaPhi", electron1->GetFourMomentum().Phi() - electron2->GetFourMomentum().Phi());

  float theta = electron1->GetFourMomentum().Phi() - electron2->GetFourMomentum().Phi();
  histogramsHandler->Fill("dielectron_theta", theta);

  float acoWeight = 1 / acoplanarityFunction->Eval(acoplanarity);

  histogramsHandler->Fill("dielectron_deltaPhi", deltaPhi);
  histogramsHandler->Fill("dielectron_deltaPhiAcoWeighted", deltaPhi, acoWeight);
  if ((float)electron1->Get("pt") > 6) {
    histogramsHandler->Fill("dielectron_deltaPhiPtGt6GeV", deltaPhi);
  }

  if (acoplanarity > 0.01) {
    histogramsHandler->Fill("dielectron_deltaPhi0p01", deltaPhi);
  }

  if (acoplanarity > 0.005) {
    histogramsHandler->Fill("dielectron_deltaPhi0p005", deltaPhi);
  }

  if (acoplanarity > 0.001) {
    histogramsHandler->Fill("dielectron_deltaPhi0p001", deltaPhi);
  }

  histogramsHandler->Fill("dielectron_acoplanarityAcoWeighted", acoplanarity, acoWeight);

  if (acoplanarity > 0.01) {
    histogramsHandler->Fill("dielectronSR_pt", dielectron.Pt());
    histogramsHandler->Fill("dielectronSR_mass", dielectron.M());
    histogramsHandler->Fill("dielectronSR_rapidity", dielectron.Rapidity());
    histogramsHandler->Fill("dielectronSR_deltaPt", deltaPt);

    histogramsHandler->Fill("goodElectronSR_pt", electron1->Get("pt"));
    histogramsHandler->Fill("goodElectronSR_eta", electron1->Get("eta"));
    histogramsHandler->Fill("goodElectronSR_phi", electron1->Get("phi"));

    histogramsHandler->Fill("goodElectronSR_pt", electron2->Get("pt"));
    histogramsHandler->Fill("goodElectronSR_eta", electron2->Get("eta"));
    histogramsHandler->Fill("goodElectronSR_phi", electron2->Get("phi"));
  }
}

void LbLHistogramsFiller::FillEventLevelHistograms(const shared_ptr<Event> event) {
  auto lblEvent = asLbLEvent(event);
  float deltaEt = lblEvent->GetDeltaEt();
  histogramsHandler->Fill("event_deltaEt", deltaEt);

  auto photons = event->GetCollection("goodPhoton");

  if (photons->size() == 2) {
    double deltaPhi = asPhoton(photons->at(0))->GetFourMomentum().DeltaPhi(asPhoton(photons->at(1))->GetFourMomentum());
    double acoplanarity = 1 - (fabs(deltaPhi) / TMath::Pi());
    float cosThetaStar = fabs(lblEvent->GetCosThetaStar(false));

    histogramsHandler->Fill("event_cosThetaStar", cosThetaStar);

    if (acoplanarity < 0.01) {
      histogramsHandler->Fill("eventSR3_cosThetaStar", cosThetaStar);
      histogramsHandler->Fill("eventSR4_cosThetaStar", cosThetaStar);
      histogramsHandler->Fill("eventSR5_cosThetaStar", cosThetaStar);
      histogramsHandler->Fill("eventSR10_cosThetaStar", cosThetaStar);

      histogramsHandler->Fill("unfoldingPhoton_costhetastar2", cosThetaStar);
      histogramsHandler->Fill("unfoldingPhoton_costhetastar3", cosThetaStar);
      histogramsHandler->Fill("unfoldingPhoton_costhetastar4", cosThetaStar);
    }
  }

  auto electrons = event->GetCollection("goodElectron");
  if (electrons->size() == 2) {
    double deltaPhi = asElectron(electrons->at(0))->GetFourMomentum().DeltaPhi(asElectron(electrons->at(1))->GetFourMomentum());
    double acoplanarity = 1 - (fabs(deltaPhi) / TMath::Pi());
    float cosThetaStar = fabs(asLbLEvent(event)->GetCosThetaStar(true));

    histogramsHandler->Fill("event_electronsCosThetaStar", cosThetaStar);

    if (acoplanarity < 0.01) {
      histogramsHandler->Fill("eventSR_electronsCosThetaStar", cosThetaStar);
    }
  }

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
  } catch (const Exception &e) {
    warn() << "Cannot fill ZDC histograms, since ZDC collection was not found." << endl;
  }
}

void LbLHistogramsFiller::Fill(const shared_ptr<Event> event) {
  FillPhotonHistograms(event);
  FillEventLevelHistograms(event);
  FillGenLevelHistograms(event);
}
