#include "CaloTower.hpp"

#include "ConfigManager.hpp"
#include "Logger.hpp"

using namespace std;

CaloTower::CaloTower(std::shared_ptr<PhysicsObject> physicsObject_) : physicsObject(physicsObject_) {
  etaEdges = {0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.870, 0.957, 1.044, 1.131,
              1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650,
              2.853, 3.000, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

  eta = Get("eta");
  phi = Get("phi");
  absEta = fabs(eta);
  iEta = GetIeta();

  auto& config = ConfigManager::GetInstance();
  config.GetMap("deadEtas", deadEtas);
  config.GetMap("caloEtaEdges", caloEtaEdges);
  config.GetMap("caloMatching", caloMatching);
  config.GetMap("detectorParams", detectorParams);
  config.GetMap("caloNoiseThresholds", caloNoiseThresholds);
  config.GetMap("caloNoiseVariables", caloNoiseVariables);

  hadronicSubdetector = GetHadronicSubdetectorName();
  electromagneticSubdetector = GetElectromagneticSubdetectorName();
}

string CaloTower::GetDetectorType() {
  if (hadronicSubdetector == "HFp" || hadronicSubdetector == "HFm") return "HF";
  if (hadronicSubdetector == "HB" || hadronicSubdetector == "HE") return "HCAL";
  if (electromagneticSubdetector == "EB" || electromagneticSubdetector == "EE") return "ECAL";
  if (hadronicSubdetector == "HadronicCrack" || electromagneticSubdetector == "ElectromagneticCrack") return "crack";
  return "unknown";
}

bool CaloTower::IsDead() {
  string hadronicSubdetector = GetHadronicSubdetectorName();
  if (deadEtas.find(hadronicSubdetector) == deadEtas.end()) return false;

  auto deadEtasForDetector = deadEtas[hadronicSubdetector];
  return find(deadEtasForDetector.begin(), deadEtasForDetector.end(), iEta) != deadEtasForDetector.end();
}

bool CaloTower::IsInHadronicCrack() {
  return (absEta >= detectorParams["crackHadron_start"] && absEta <= detectorParams["crackHadron_end"]);
}

bool CaloTower::IsInElectromagneticCrack() { return (absEta >= detectorParams["crack_start"] && absEta <= detectorParams["crack_end"]); }

string CaloTower::GetHadronicSubdetectorName() {
  if (eta > -caloEtaEdges["maxHF"] && eta < -caloEtaEdges["minHF"]) return "HFm";
  if (eta > caloEtaEdges["minHF"] && eta < caloEtaEdges["maxHF"]) return "HFp";
  if (absEta < caloEtaEdges["maxHB"]) return "HB";
  if (absEta > caloEtaEdges["minHE"] && absEta < caloEtaEdges["maxHE"]) return "HE";
  if (IsInHadronicCrack()) return "HadronicCrack";

  error() << "ERROR - could not determine calo tower had sub-det! Eta: " << eta << endl;
  return "";
}

string CaloTower::GetElectromagneticSubdetectorName() {
  if (eta > -caloEtaEdges["maxHF"] && eta < -caloEtaEdges["minHF"]) return "HFm";
  if (eta > caloEtaEdges["minHF"] && eta < caloEtaEdges["maxHF"]) return "HFp";
  if (absEta < caloEtaEdges["maxEB"]) return "EB";
  if (absEta > caloEtaEdges["minEE"] && absEta < caloEtaEdges["maxEE"]) return "EE";
  if (IsInElectromagneticCrack()) return "ElectromagneticCrack";

  error() << "ERROR - could not determine calo tower EM sub-det!" << endl;
  return "";
}

int CaloTower::GetIeta() {
  int iEta = 1;
  while (absEta > etaEdges[iEta] && iEta < etaEdges.size() - 1) iEta++;
  if (eta < 0) iEta = -iEta;
  return iEta;
}

bool CaloTower::IsInHEM() {
  return (eta > detectorParams["hem_etaStart"] && eta < detectorParams["hem_etaEnd"] && phi > detectorParams["hem_phiStart"] &&
          phi < detectorParams["hem_phiEnd"]);
}

bool CaloTower::OverlapsWithOtherObjects(shared_ptr<PhysicsObjects> otherObjects) {
  float phi = Get("phi");

  float maxDeltaEta = caloMatching["maxDeltaEta_" + electromagneticSubdetector];
  float maxDeltaPhi = caloMatching["maxDeltaPhi_" + electromagneticSubdetector];

  for (auto otherObject : *otherObjects) {
    float deltaEta = fabs((float)otherObject->Get("SCEta") - eta);
    float deltaPhi = fabs(TVector2::Phi_mpi_pi((float)otherObject->Get("SCPhi") - phi));
    if (deltaEta < maxDeltaEta && deltaPhi < maxDeltaPhi) return true;
  }
  return false;
}

bool CaloTower::IsEtaAboveLimit() {
  return (absEta > detectorParams["caloTower_etaMax"] && hadronicSubdetector != "HFp" && hadronicSubdetector != "HFm");
}

bool CaloTower::IsHadronicEnergyAboveNoiseThreshold() {
  if (hadronicSubdetector == "HadronicCrack") return false;
  return (float)Get(caloNoiseVariables[hadronicSubdetector]) > caloNoiseThresholds[hadronicSubdetector];
}

bool CaloTower::IsElectromagneticEnergyAboveNoiseThreshold() {
  if (electromagneticSubdetector == "ElectromagneticCrack") return false;
  return (float)Get(caloNoiseVariables[electromagneticSubdetector]) > caloNoiseThresholds[electromagneticSubdetector];
}