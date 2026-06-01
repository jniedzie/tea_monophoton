#ifndef Electron_hpp
#define Electron_hpp

#include "Helpers.hpp"
#include "PhysicsObject.hpp"

class Electron;
typedef Collection<std::shared_ptr<Electron>> Electrons;

class Electron {
 public:
  Electron(std::shared_ptr<PhysicsObject> physicsObject_);

  auto Get(std::string branchName) { return physicsObject->Get(branchName); }

  template <typename T>
  T GetAs(std::string branchName) { return physicsObject->GetAs<float>(branchName); }
  std::string GetOriginalCollection() { return physicsObject->GetOriginalCollection(); }
  void Reset() { physicsObject->Reset(); }

  bool PassesPtCuts();
  bool IsEtaAboveLimit();
  bool IsInCrack();
  bool IsInHEM();
  bool PassesMissingHitsCuts();
  bool PassesHoverE();
  bool PassesDeltaEtaAtVertex();
  bool PassesIsolationCuts();

  float GetPt() { return GetAs<float>("pt"); }
  float GetEta() { return eta; }
  float GetPhi() { return phi; }

  int GetCharge() {
    try {
      return Get("charge");
    } catch (const std::exception& e) {
      return (int)Get("pid") > 0 ? 1 : -1;
    }
  }

  TLorentzVector GetFourMomentum() {
    TLorentzVector fourMomentum;
    try {
      fourMomentum.SetPtEtaPhiM(Get("pt"), eta, phi, 0.000511);
    } catch (const std::exception& e) {
      fourMomentum.SetPtEtaPhiM(Get("et"), eta, phi, 0.000511);
    }
    return fourMomentum;
  }

 private:
  std::shared_ptr<PhysicsObject> physicsObject;

  std::map<std::string, float> electronCuts, detectorParams, caloEtaEdges;
  float eta, phi, absEta, etaSC, phiSC, absEtaSC;

  std::string detRegion;
};

#endif /* Electron_hpp */
