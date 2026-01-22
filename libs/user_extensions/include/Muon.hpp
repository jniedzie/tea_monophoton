#ifndef Muon_hpp
#define Muon_hpp

#include "Helpers.hpp"
#include "PhysicsObject.hpp"

class Muon;
typedef Collection<std::shared_ptr<Muon>> Muons;

class Muon {
 public:
  Muon(std::shared_ptr<PhysicsObject> physicsObject_);

  auto Get(std::string branchName, bool verbose=true) { return physicsObject->Get(branchName, verbose); }
  
  template <typename T>
  T GetAs(std::string branchName) { return physicsObject->GetAs<T>(branchName); }
  std::string GetOriginalCollection() { return physicsObject->GetOriginalCollection(); }
  void Reset() { physicsObject->Reset(); }

  bool PassesPtCuts();
  bool IsEtaAboveLimit();

 private:
  std::shared_ptr<PhysicsObject> physicsObject;

  std::map<std::string, float> muonCuts, detectorParams, caloEtaEdges;
  float eta, phi, pt, absEta;
};

#endif /* Muon_hpp */
