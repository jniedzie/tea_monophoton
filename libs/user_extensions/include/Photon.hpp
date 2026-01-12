#ifndef Photon_hpp
#define Photon_hpp

#include "Helpers.hpp"
#include "PhysicsObject.hpp"

class Photon;
typedef Collection<std::shared_ptr<Photon>> Photons;

class Photon {
 public:
  Photon(std::shared_ptr<PhysicsObject> physicsObject_);

  auto Get(std::string branchName) { return physicsObject->Get(branchName); }
  
  template <typename T>
  T GetAs(std::string branchName) { return physicsObject->GetAs<T>(branchName); }
  std::string GetOriginalCollection() { return physicsObject->GetOriginalCollection(); }
  void Reset() { physicsObject->Reset(); }

  float GetSwissCross();

  bool IsInHEM();
  bool IsEtaAboveLimit();
  bool IsInCrack();
  bool IsInHotSpot();
  bool PassesShowerShape();
  bool PassesHoverE();
  bool PassesSwissCross();
  bool PassesEtCuts();
  bool PassesSeedTimeCuts();
  bool PassesConversionCuts();

  bool OverlapsWithOtherObjects(std::shared_ptr<PhysicsObjects> otherObjects);

  TLorentzVector GetFourMomentum();

  float GetVerticalOverCentralEnergy();
  float GetHorizontalOverCentralEnergy();

  float GetHorizontalImbalance();
  float GetVerticalImbalance();

 private:
  std::shared_ptr<PhysicsObject> physicsObject;

  std::map<std::string, float> photonCuts, detectorParams, caloEtaEdges, electronPhotonMatching;
  float eta, phi, absEta, etaSC, phiSC, absEtaSC;

  std::string detRegion;

  std::map<std::string, std::vector<float>> hotSpots;
};

#endif /* Photon_hpp */
