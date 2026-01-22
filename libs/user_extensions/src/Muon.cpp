#include "Muon.hpp"

#include "ConfigManager.hpp"

using namespace std;

Muon::Muon(shared_ptr<PhysicsObject> physicsObject_) : physicsObject(physicsObject_) {
  auto &config = ConfigManager::GetInstance();
  config.GetMap("muonCuts", muonCuts);
  config.GetMap("detectorParams", detectorParams);
  config.GetMap("caloEtaEdges", caloEtaEdges);

  try {
    eta = Get("Eta", false);
  } catch (Exception &e) {
    eta = Get("eta");
  }
  try {
    phi = Get("Phi", false);
  } catch (Exception &e) {
    phi = Get("phi");
  }
  try {
    pt = Get("Pt", false);
  } catch (Exception &e) {
    pt = Get("pt");
  }

  absEta = fabs(eta);
}

bool Muon::PassesPtCuts() { return pt > muonCuts["min_pt"]; }

bool Muon::IsEtaAboveLimit() { return absEta > muonCuts["max_absEta"]; }
