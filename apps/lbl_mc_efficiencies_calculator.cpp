#include "ConfigManager.hpp"
#include "CutFlowManager.hpp"
#include "EventProcessor.hpp"
#include "EventReader.hpp"
#include "EventWriter.hpp"
#include "ExtensionsHelpers.hpp"
#include "HistogramsFiller.hpp"
#include "HistogramsHandler.hpp"
#include "LbLObjectsManager.hpp"
#include "LbLSelections.hpp"
#include "Logger.hpp"
#include "Profiler.hpp"
#include "UserExtensionsHelpers.hpp"

using namespace std;

bool runForPhotons = false;

void CheckArgs(int argc, char **argv) {
  if (argc != 2 && argc != 4) {
    fatal() << "Usage: " << argv[0] << " config_path" << endl;
    fatal() << "or" << endl;
    fatal() << argv[0] << " config_path input_path output_path" << endl;
    exit(1);
  }
}

float GetAcoplanarity(const TLorentzVector &vec1, const TLorentzVector &vec2) { return (1 - (fabs(vec1.DeltaPhi(vec2)) / TMath::Pi())); }

bool GoodID(const shared_ptr<Photon> photon) {
  if (!photon->PassesShowerShape()) return false;
  if (!photon->PassesHoverE()) return false;
  if (!photon->PassesSeedTimeCuts()) return false;
  if (!photon->PassesSwissCross()) return false;
  return true;
}

bool GoodID(const shared_ptr<Electron> electron) {
  if (!electron->PassesMissingHitsCuts()) return false;
  if (!electron->PassesHoverE()) return false;
  if (!electron->PassesDeltaEtaAtVertex()) return false;
  if (!electron->PassesIsolationCuts()) return false;
  return true;
}

bool GoodReco(const shared_ptr<Electron> electron) {
  if (!electron->PassesPtCuts()) return false;
  if (electron->IsInCrack()) return false;
  if (electron->IsEtaAboveLimit()) return false;
  if (electron->IsInHEM()) return false;
  return true;
}

bool GoodReco(const shared_ptr<Photon> photon) {
  if (!photon->PassesConversionCuts()) return false;
  if (photon->IsEtaAboveLimit()) return false;
  if (photon->IsInCrack()) return false;
  if (photon->IsInHEM()) return false;
  if (!photon->PassesEtCuts()) return false;
  return true;
}

bool InAcceptance(const TLorentzVector &vec1, const TLorentzVector &vec2) {
  auto diphoton = vec1 + vec2;
  if (fabs(vec1.Eta()) > 2.2 || fabs(vec2.Eta()) > 2.2) return false;
  if (vec1.Pt() < 2.0 || vec2.Pt() < 2.0) return false;
  if (diphoton.M() < 5.0) return false;
  if (diphoton.Pt() > 1.0) return false;
  if (GetAcoplanarity(vec1, vec2) > 0.01) return false;

  return true;
}

pair<shared_ptr<PhysicsObject>, shared_ptr<PhysicsObject>> GetMatchedRecoParticles(TLorentzVector &particle1, TLorentzVector &particle2,
                                                                                   shared_ptr<PhysicsObjects> recoParticles,
                                                                                   bool isPhoton = true) {
  shared_ptr<PhysicsObject> recoParticle1 = nullptr;
  shared_ptr<PhysicsObject> recoParticle2 = nullptr;

  float deltaRcut = 0.5;

  for (auto physicsObject : *recoParticles) {
    TLorentzVector recoVec;

    if (isPhoton) {
      if (!GoodReco(asPhoton(physicsObject))) continue;
      recoVec = asPhoton(physicsObject)->GetFourMomentum();
    } else {
      if (!GoodReco(asElectron(physicsObject))) continue;
      recoVec = asElectron(physicsObject)->GetFourMomentum();
    }

    if (particle1.DeltaR(recoVec) < deltaRcut) {
      recoParticle1 = physicsObject;
    }
    if (particle2.DeltaR(recoVec) < deltaRcut) {
      recoParticle2 = physicsObject;
    }
  }
  return {recoParticle1, recoParticle2};
}

int main(int argc, char **argv) {
  CheckArgs(argc, argv);
  ConfigManager::Initialize(argv[1]);
  auto &config = ConfigManager::GetInstance();
  if (argc == 4) {
    config.SetInputPath(argv[2]);
    config.SetHistogramsOutputPath(argv[3]);
  }

  auto eventReader = make_shared<EventReader>();
  auto eventProcessor = make_unique<EventProcessor>();
  auto lblSelections = make_unique<LbLSelections>();
  auto lblObjectsManager = make_unique<LbLObjectsManager>();

  int n_total = eventReader->GetNevents();
  int n_acc = 0;
  int n_reco = 0;
  int n_id = 0;
  int n_good_for_trig = 0;
  int n_trig = 0;
  int n_che = 0;
  int n_nee = 0;

  // auto histDeltaRphotons = new TH1D("histDeltaRphotons", "histDeltaRphotons", 100, 0, 5);
  // auto histDeltaRelectrons = new TH1D("histDeltaRelectrons", "histDeltaRelectrons", 100, 0, 5);

  // for (int iEvent = 0; iEvent < eventReader->GetNevents(); iEvent++) {
  //   auto event = eventReader->GetEvent(iEvent);

  //   auto mcParticles = event->GetCollection("genParticle");
  //   auto genPhotons = GetGenParticles(event, 22);
  //   auto genElectrons = GetGenParticles(event, 11);
  //   auto recoPhotons = event->GetCollection("photon");
  //   auto recoElectrons = event->GetCollection("electron");

  //   if (genPhotons->size() == 2) {
  //     auto photon1 = GetFourMomentum(genPhotons->at(0));
  //     auto photon2 = GetFourMomentum(genPhotons->at(1));
  //     if (InAcceptance(photon1, photon2)) {
  //       for(auto physObject : *recoPhotons){
  //         auto photon = asPhoton(physObject);
  //         if(!GoodReco(photon)) continue;
  //         histDeltaRphotons->Fill(photon1.DeltaR(photon->GetFourMomentum()));
  //         histDeltaRphotons->Fill(photon2.DeltaR(photon->GetFourMomentum()));
  //       }
  //     }
  //   }

  //   if (genElectrons->size() == 2) {
  //     auto electron1 = GetFourMomentum(genElectrons->at(0));
  //     auto electron2 = GetFourMomentum(genElectrons->at(1));
  //     if (InAcceptance(electron1, electron2)) {
  //       for(auto physObject : *recoElectrons){
  //         auto electron = asElectron(physObject);
  //         if(!GoodReco(electron)) continue;
  //         histDeltaRelectrons->Fill(electron1.DeltaR(electron->GetFourMomentum()));
  //         histDeltaRelectrons->Fill(electron2.DeltaR(electron->GetFourMomentum()));
  //       }
  //     }
  //   }
  // }

  // histDeltaRphotons->SaveAs("../histDeltaRphotons.root");
  // histDeltaRelectrons->SaveAs("../histDeltaRelectrons.root");
  // exit(0);

  for (int iEvent = 0; iEvent < eventReader->GetNevents(); iEvent++) {
    auto event = eventReader->GetEvent(iEvent);
    lblObjectsManager->InsertGoodPhotonsCollection(event);
    lblObjectsManager->InsertGoodElectronsCollection(event);
    lblObjectsManager->InsertGoodTracksCollection(event);
    lblObjectsManager->InsertGoodMuonsCollection(event, false);
    lblObjectsManager->InsertGenPhotonsCollection(event);
    lblObjectsManager->InsertGenElectronsCollection(event);

    auto genParticles = event->GetCollection(runForPhotons ? "genPhoton" : "genElectron");
    auto recoParticles = event->GetCollection(runForPhotons ? "photon" : "electron");

    if (genParticles->size() != 2) continue;

    // check acceptance

    TLorentzVector genParticle1, genParticle2;

    if (runForPhotons) {
      genParticle1 = asPhoton(genParticles->at(0))->GetFourMomentum();
      genParticle2 = asPhoton(genParticles->at(1))->GetFourMomentum();
    } else {
      genParticle1 = asElectron(genParticles->at(0))->GetFourMomentum();
      genParticle2 = asElectron(genParticles->at(1))->GetFourMomentum();
    }
    if (!InAcceptance(genParticle1, genParticle2)) continue;
    n_acc++;

    // check reco
    auto [recoObj1, recoObj2] = GetMatchedRecoParticles(genParticle1, genParticle2, recoParticles, runForPhotons);
    if(!recoObj1 || !recoObj2) continue;
    
    TLorentzVector recoParticle1, recoParticle2;
    if (runForPhotons) {
      recoParticle1 = asPhoton(recoParticles->at(0))->GetFourMomentum();
      recoParticle2 = asPhoton(recoParticles->at(1))->GetFourMomentum();
    } else {
      recoParticle1 = asElectron(recoParticles->at(0))->GetFourMomentum();
      recoParticle2 = asElectron(recoParticles->at(1))->GetFourMomentum();
    }
    
    if (!InAcceptance(recoParticle1, recoParticle2)) continue;
    n_reco++;

    // check id
    if(runForPhotons){
      if (!GoodID(asPhoton(recoObj1)) || !GoodID(asPhoton(recoObj2))) continue;
    } else {
      if (!GoodID(asElectron(recoObj1)) || !GoodID(asElectron(recoObj2))) continue;
    }
    n_id++;

    // check trigger
    if (!(int)event->Get("DoubleEG2")) continue;
    n_trig++;

    // check charged exclusivity
    if(runForPhotons){
      if (!lblSelections->PassesChargedExclusivity(event)) continue;
    }
    else{
      if(!lblSelections->PassesDielectronChargedExclusivity(event)) continue;
    }
    n_che++;

    // check neutral exclusivity
    if (!lblSelections->PassesNeutralExclusivity(event)) continue;
    n_nee++;
  }

  info() << "Initial number of events: " << n_total << endl;
  info() << "Number of accepted events: " << n_acc << endl;
  info() << "\nNumber of reconstructed events: " << n_reco << endl;
  info() << "\nNumber of identified events: " << n_id << endl;
  info() << "\nNumber of triggered events: " << n_trig << endl;
  info() << "\nNumber of charged exclusive events: " << n_che << endl;
  info() << "\nNumber of neutral exclusive events: " << n_nee << endl;

  float n_acc_err = sqrt(n_acc);
  float n_reco_err = sqrt(n_reco);
  float n_id_err = sqrt(n_id);
  float n_trig_err = sqrt(n_trig);
  float n_che_err = sqrt(n_che);
  float n_nee_err = sqrt(n_nee);

  float eff_reco = (float)n_reco / n_acc;
  float eff_id = (float)n_id / n_reco;
  float eff_trig = (float)n_trig / n_id;
  float eff_che = (float)n_che / n_trig;
  float eff_nee = (float)n_nee / n_che;
  float total_eff = (float)n_nee / n_acc;

  float eff_reco_err = eff_reco * sqrt(pow(n_reco_err / n_reco, 2) + pow(n_acc_err / n_acc, 2));
  float eff_id_err = eff_id * sqrt(pow(n_id_err / n_id, 2) + pow(n_reco_err / n_reco, 2));
  float eff_trig_err = eff_trig * sqrt(pow(n_trig_err / n_trig, 2) + pow(n_id_err / n_id, 2));
  float eff_che_err = eff_che * sqrt(pow(n_che_err / n_che, 2) + pow(n_trig_err / n_trig, 2));
  float eff_nee_err = eff_nee * sqrt(pow(n_nee_err / n_nee, 2) + pow(n_che_err / n_che, 2));
  float total_eff_err = total_eff * sqrt(pow(n_nee_err / n_nee, 2) + pow(n_acc_err / n_acc, 2));

  info() << "Eff_acc: " << (float)n_acc / n_total << endl;

  info() << "\n\nTwo-photon Eff_reco: " << eff_reco << "+/-" << eff_reco_err << " (single photon: " << sqrt(eff_reco) << ")" << endl;
  info() << "Two-photon Eff_id: " << eff_id << "+/-" << eff_id_err << " (single photon: " << sqrt(eff_id) << ")" << endl;
  info() << "Eff_trig: " << eff_trig << "+/-" << eff_trig_err << endl;
  info() << "Eff_che: " << eff_che << "+/-" << eff_che_err << endl;
  info() << "Eff_nee: " << eff_nee << "+/-" << eff_nee_err << endl;

  info() << "Product of efficiencies (squared reco&ID): " << eff_reco * eff_reco * eff_id * eff_id * eff_che * eff_nee * eff_trig << endl;
  info() << "Product of efficiencies (reco&ID not squared): " << eff_reco * eff_id * eff_che * eff_nee * eff_trig << endl;

  info() << "Total efficiency: " << total_eff << "+/-" << total_eff_err << endl;

  auto &logger = Logger::GetInstance();
  logger.Print();

  return 0;
}