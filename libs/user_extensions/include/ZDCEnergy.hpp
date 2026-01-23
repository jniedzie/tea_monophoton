#ifndef ZDCEnergy_hpp
#define ZDCEnergy_hpp

#include "Helpers.hpp"
#include "PhysicsObject.hpp"

class ZDCEnergy;
typedef Collection<std::shared_ptr<ZDCEnergy>> ZDCEnergys;

class ZDCEnergy {
 public:
  ZDCEnergy(std::shared_ptr<PhysicsObject> physicsObject_) : physicsObject(physicsObject_) {}

  auto Get(std::string branchName) { return physicsObject->Get(branchName); }
  
  template <typename T>
  T GetAs(std::string branchName) { return physicsObject->GetAs<T>(branchName); }
  std::string GetOriginalCollection() { return physicsObject->GetOriginalCollection(); }
  void Reset() { physicsObject->Reset(); }

  float GetEnergy() { return Get("energy");}
  int GetSide() { return Get("zSide");}

  float GetSaturation() { return Get("saturation");}
  int GetSection() { return Get("section");}
  int GetChannel() { return Get("channel");}

  // implement a function (maybe an operator<<) that would allow to print all info either to cout or to a file stream
  friend std::ostream& operator<<(std::ostream& os, ZDCEnergy& zdcEnergy) {
    os << "ZDCEnergy: " << std::endl;
    os << "\tEnergy: " << zdcEnergy.GetEnergy() << std::endl;
    os << "\tSide: " << zdcEnergy.GetSide() << std::endl;
    os << "\tSaturation: " << zdcEnergy.GetSaturation() << std::endl;
    os << "\tSection: " << zdcEnergy.GetSection() << std::endl;
    os << "\tChannel: " << zdcEnergy.GetChannel() << std::endl;
    return os;
  }

 private:
  std::shared_ptr<PhysicsObject> physicsObject;
};

#endif /* ZDCEnergy_hpp */
