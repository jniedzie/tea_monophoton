#include "LbLObjectsManager.hpp"

#include "ConfigManager.hpp"
#include "Logger.hpp"

using namespace std;

LbLObjectsManager::LbLObjectsManager() {
  auto& config = ConfigManager::GetInstance();
  config.GetMap("detectorParams", detectorParams);
  config.GetMap("caloEtaEdges", caloEtaEdges);
  config.GetVector("knownPids", knownPids);
}

bool LbLObjectsManager::IsGoodPhoton(const shared_ptr<Photon> photon, shared_ptr<map<string, int>> cutFlow) {
  if (cutFlow) cutFlow->at("00_initial")++;
  if (!photon->PassesConversionCuts()) return false;
  if (cutFlow) cutFlow->at("01_conversionCuts")++;
  if (!photon->PassesEtCuts()) return false;
  if (cutFlow) cutFlow->at("02_etCuts")++;
  if (!photon->PassesSwissCross()) return false;
  if (cutFlow) cutFlow->at("03_swissCross")++;
  if (photon->IsEtaAboveLimit()) return false;
  if (cutFlow) cutFlow->at("04_etaCuts")++;
  if (photon->IsInCrack()) return false;
  if (cutFlow) cutFlow->at("05_crackCuts")++;
  if (photon->IsInHotSpot()) return false;
  if (cutFlow) cutFlow->at("06_hotSpotCuts")++;
  if (photon->IsInHEM()) return false;
  if (cutFlow) cutFlow->at("07_HEMCuts")++;
  if (!photon->PassesShowerShape()) return false;
  if (cutFlow) cutFlow->at("08_showerShape")++;
  if (!photon->PassesHoverE()) return false;
  if (cutFlow) cutFlow->at("09_hoverE")++;
  if (!photon->PassesSeedTimeCuts()) return false;
  if (cutFlow) cutFlow->at("10_seedTime")++;

  return true;
}

bool LbLObjectsManager::IsGoodElectron(const shared_ptr<Electron> electron) {
  if (!electron->PassesPtCuts()) return false;
  if (electron->IsInCrack()) return false;
  if (electron->IsEtaAboveLimit()) return false;
  if (electron->IsInHEM()) return false;
  if (!electron->PassesMissingHitsCuts()) return false;
  if (!electron->PassesHoverE()) return false;
  if (!electron->PassesDeltaEtaAtVertex()) return false;
  if (!electron->PassesIsolationCuts()) return false;

  return true;
}

bool LbLObjectsManager::IsGoodTrack(const shared_ptr<Track> track) {
  if (!track->PassesPtCuts()) return false;
  if (track->IsEtaAboveLimit()) return false;
  if (!track->PassesValidHitsCuts()) return false;
  if (!track->PassesChi2Cuts()) return false;
  if (!track->PassesDistanceToPVCuts()) return false;

  return true;
}

bool LbLObjectsManager::IsGoodMuon(const shared_ptr<Muon> muon) {
  if (!muon->PassesPtCuts()) return false;
  if (muon->IsEtaAboveLimit()) return false;

  return true;
}

void LbLObjectsManager::InsertGoodPhotonsCollection(shared_ptr<Event> event, shared_ptr<map<string, int>> cutFlow) {
  auto photons = event->GetCollection("photon");
  auto goodPhotons = make_shared<PhysicsObjects>();

  for (auto physicsObject : *photons) {
    auto photon = asPhoton(physicsObject);
    if (!IsGoodPhoton(photon, cutFlow)) continue;
    goodPhotons->push_back(physicsObject);
  }

  event->AddCollection("goodPhoton", goodPhotons);
}

void LbLObjectsManager::InsertGoodElectronsCollection(shared_ptr<Event> event) {
  auto electrons = event->GetCollection("electron");
  auto goodElectrons = make_shared<PhysicsObjects>();

  for (auto physicsObject : *electrons) {
    auto electron = asElectron(physicsObject);
    if (!IsGoodElectron(electron)) continue;
    goodElectrons->push_back(physicsObject);
  }

  event->AddCollection("goodElectron", goodElectrons);
}

void LbLObjectsManager::InsertGoodTracksCollection(shared_ptr<Event> event) {
  auto tracks = event->GetCollection("track");
  auto goodTracks = make_shared<PhysicsObjects>();

  for (auto physicsObject : *tracks) {
    auto track = asTrack(physicsObject);
    if (!IsGoodTrack(track)) continue;
    goodTracks->push_back(physicsObject);
  }

  event->AddCollection("goodTrack", goodTracks);
}

void LbLObjectsManager::InsertGoodMuonsCollection(shared_ptr<Event> event) {
  shared_ptr<PhysicsObjects> muons;

  try {
    muons = event->GetCollection("muon");
  } catch (const Exception& e) {
    error() << "No muon collection found in event. Will not insert goodMuon collection." << endl;
    return;
  }
  auto goodMuons = make_shared<PhysicsObjects>();

  for (auto physicsObject : *muons) {
    auto muon = asMuon(physicsObject);
    if (!IsGoodMuon(muon)) continue;
    goodMuons->push_back(physicsObject);
  }

  event->AddCollection("goodMuon", goodMuons);
}

void LbLObjectsManager::InsertGenPhotonsCollection(shared_ptr<Event> event) {
  auto genParticles = event->GetCollection("genParticle");
  auto genPhotons = GetGenParticles(event, 22);
  event->AddCollection("genPhoton", genPhotons);
}

void LbLObjectsManager::InsertGenElectronsCollection(shared_ptr<Event> event) {
  auto genParticles = event->GetCollection("genParticle");
  auto genPhotons = GetGenParticles(event, 11);
  event->AddCollection("genElectron", genPhotons);
}

shared_ptr<PhysicsObjects> LbLObjectsManager::GetGenParticles(const shared_ptr<Event> event, int pid) {
  auto mcParticles = event->GetCollection("genParticle");
  auto genParticles = make_shared<PhysicsObjects>();

  for (auto particle : *mcParticles) {
    int particlePid = GetParticlePid(particle);
    if (abs(particlePid) == pid) genParticles->push_back(particle);
  }
  return genParticles;
}

int LbLObjectsManager::GetParticlePid(const shared_ptr<PhysicsObject> particle) {
  int particlePid = particle->Get("pid");
  int convertedPid;

  if (find(knownPids.begin(), knownPids.end(), abs(particlePid)) != knownPids.end()) {
    return particlePid;
  } else {
    // this is a hack needed if pid was stored as float in the tree
  float* floatPtr = reinterpret_cast<float*>(&particlePid);
    float floatValue = *floatPtr;
    convertedPid = round(floatValue);

    if (find(knownPids.begin(), knownPids.end(), abs(convertedPid)) != knownPids.end()) {
      return convertedPid;
    }
  }

  fatal() << "Unknown PID: " << particlePid << ", converted: " << convertedPid
          << ". If one of them makes sense, add it to knownPids in the config. Otherwise, there's a serious problem in decoding it from "
             "the ntuple."
          << endl;
  exit(1);
}