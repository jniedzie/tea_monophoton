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
  float GetEt() { return GetAs<float>("et"); }
  float GetEta() { return eta; }
  float GetPhi() { return phi; }
  float GetSeedTime() { return GetAs<float>("seedTime"); }

  float GetEnergyTop() { return energyTop; }
  float GetEnergyBottom() { return energyBottom; }
  float GetEnergyLeft() { return energyLeft; }
  float GetEnergyRight() { return energyRight; }
  float GetEnergyCentral() { return energyCentral; }
  float GetMinEnergy() { return energyMin; }

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

  // implement a function (maybe an operator<<) that would allow to print all info either to cout or to a file stream
  friend std::ostream& operator<<(std::ostream& os, Photon& photon) {
    os << "Photon: " << std::endl;
    os << "\tEt: " << photon.GetEt() << std::endl;
    os << "\tEta: " << photon.GetEta() << std::endl;
    os << "\tPhi: " << photon.GetPhi() << std::endl;
    os << "\tSeed Time: " << photon.GetSeedTime() << std::endl;
    os << "\tEnergy Top: " << photon.GetEnergyTop() << std::endl;
    os << "\tEnergy Bottom: " << photon.GetEnergyBottom() << std::endl;
    os << "\tEnergy Left: " << photon.GetEnergyLeft() << std::endl;
    os << "\tEnergy Right: " << photon.GetEnergyRight() << std::endl;
    os << "\tEnergy Central: " << photon.GetEnergyCentral() << std::endl;
    os << "\tMin Energy: " << photon.GetMinEnergy() << std::endl;
    os << "\tSwiss Cross: " << photon.GetSwissCross() << std::endl;
    os << "\tDetector Region: " << photon.detRegion << std::endl;
    os << "\tVertical Over Central Energy: " << photon.GetVerticalOverCentralEnergy() << std::endl;
    os << "\tHorizontal Over Central Energy: " << photon.GetHorizontalOverCentralEnergy() << std::endl;
    os << "\tHorizontal Imbalance: " << photon.GetHorizontalImbalance() << std::endl;
    os << "\tVertical Imbalance: " << photon.GetVerticalImbalance() << std::endl;
    os << "\tSigma Eta 2012: " << photon.GetAs<float>("sigmaEta2012") << std::endl;
    os << "\tSigma Ieta Ieta 2012: " << photon.GetAs<float>("sigmaIEtaIEta2012") << std::endl;
    os << "\tSCEta Width: " << photon.GetAs<float>("SCEtaWidth") << std::endl;
    os << "\tSCPhi Width: " << photon.GetAs<float>("SCPhiWidth") << std::endl;
    os << "\tH over E: " << photon.GetAs<float>("hOverE") << std::endl;


    return os;
  }

 private:
  std::shared_ptr<PhysicsObject> physicsObject;

  std::map<std::string, float> photonCuts, detectorParams, caloEtaEdges, electronPhotonMatching;
  float eta, phi, absEta, etaSC, phiSC, absEtaSC;
  float energyTop, energyBottom, energyLeft, energyRight, energyCentral, energyMin;

  std::string detRegion;

  std::map<std::string, std::vector<float>> hotSpots;
};

#endif /* Photon_hpp */
