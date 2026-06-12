#include "Electron.hpp"

#include "ConfigManager.hpp"

using namespace std;

Electron::Electron(std::shared_ptr<PhysicsObject> physicsObject_) : physicsObject(physicsObject_) {
  auto& config = ConfigManager::GetInstance();
  config.GetMap("electronCuts", electronCuts);
  config.GetMap("detectorParams", detectorParams);
  config.GetMap("caloEtaEdges", caloEtaEdges);

  eta = Get("eta");
  phi = Get("phi");
  absEta = fabs(eta);

  try {
    etaSC = Get("SCEta", false);
    phiSC = Get("SCPhi", false);
    absEtaSC = fabs(etaSC);
  } catch (const std::exception& e) {
    warn() << "No SC eta/phi found for electron. Using regular eta/phi instead." << endl;
    etaSC = eta;
    phiSC = phi;
    absEtaSC = absEta;
  }

  detRegion = absEta < caloEtaEdges["maxEB"] ? "barrel" : "endcap";
}

bool Electron::PassesPtCuts() { return (float)Get("pt") > electronCuts["min_pt"]; }

bool Electron::IsEtaAboveLimit() { return absEtaSC > electronCuts["max_absEtaSC"]; }

bool Electron::IsInCrack() { return (absEtaSC > detectorParams["crack_start"] && absEtaSC < detectorParams["crack_end"]); }

bool Electron::IsInHEM() {
  return (etaSC > detectorParams["hem_etaStart"] && etaSC < detectorParams["hem_etaEnd"] && phiSC > detectorParams["hem_phiStart"] &&
          phiSC < detectorParams["hem_phiEnd"]);
}

bool Electron::PassesMissingHitsCuts() { return (int)Get("nMissHits") <= electronCuts["max_nMissingHits"]; }

bool Electron::PassesHoverE() { return (float)Get("hOverE") < electronCuts["max_hOverE"]; }

bool Electron::PassesDeltaEtaAtVertex() { return abs((float)Get("deltaEtaAtVertex")) < electronCuts["max_deltaEtaAtVertex"]; }

bool Electron::PassesIsolationCuts() {
  if ((float)Get("PFChIso") > electronCuts["max_PFChIso_" + detRegion]) return false;
  if ((float)Get("PFPhoIso") > electronCuts["max_PFPhoIso_" + detRegion]) return false;
  if ((float)Get("PFNeuIso") > electronCuts["max_PFNeuIso_" + detRegion]) return false;

  return true;
}