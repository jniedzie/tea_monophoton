#ifndef CaloTower_hpp
#define CaloTower_hpp

#include "Helpers.hpp"
#include "PhysicsObject.hpp"

class CaloTower;
typedef Collection<std::shared_ptr<CaloTower>> CaloTowers;

class CaloTower {
 public:
  CaloTower(std::shared_ptr<PhysicsObject> physicsObject_);

  auto Get(std::string branchName) { return physicsObject->Get(branchName); }

  template <typename T>
  T GetAs(std::string branchName) { return physicsObject->GetAs<T>(branchName); }
  std::string GetOriginalCollection() { return physicsObject->GetOriginalCollection(); }
  void Reset() { physicsObject->Reset(); }

  float GetEta() { return eta; }
  float GetPhi() { return phi; }
  
  // Can return: ECAL, HCAL, HF, crack, unknown
  std::string GetDetectorType();

  std::string GetHadronicSubdetectorName();
  std::string GetElectromagneticSubdetectorName();

  bool IsDead();
  bool IsInHadronicCrack();
  bool IsInElectromagneticCrack();
  bool IsInHEM();
  bool IsEtaAboveLimit();
  bool IsHadronicEnergyAboveNoiseThreshold();
  bool IsElectromagneticEnergyAboveNoiseThreshold();

  bool OverlapsWithOtherObjects(std::shared_ptr<PhysicsObjects> otherObjects);

 private:
  std::shared_ptr<PhysicsObject> physicsObject;

  std::string hadronicSubdetector, electromagneticSubdetector;

  float eta, phi, absEta, iEta;
  std::vector<float> etaEdges;

  std::map<std::string, std::vector<int>> deadEtas;
  std::map<std::string, float> caloEtaEdges, caloMatching, detectorParams, caloNoiseThresholds;
  std::map<std::string, std::string> caloNoiseVariables;

  int GetIeta();
};

#endif /* CaloTower_hpp */
