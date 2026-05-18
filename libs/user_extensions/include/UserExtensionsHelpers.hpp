#ifndef UserExtensionsHelpers_hpp
#define UserExtensionsHelpers_hpp

#include "ExtensionsHelpers.hpp"
#include "CaloTower.hpp"
#include "Photon.hpp"
#include "Electron.hpp"
#include "Track.hpp"
#include "Muon.hpp"
#include "ZDCEnergy.hpp"
#include "LbLEvent.hpp"
#include "PhysicsObject.hpp"
#include "Event.hpp"



inline std::shared_ptr<CaloTower> asCaloTower(const std::shared_ptr<PhysicsObject> physicsObject) {
  return std::make_shared<CaloTower>(physicsObject);
}

inline std::shared_ptr<Photon> asPhoton(const std::shared_ptr<PhysicsObject> physicsObject) {
  return std::make_shared<Photon>(physicsObject);
}

inline std::shared_ptr<Electron> asElectron(const std::shared_ptr<PhysicsObject> physicsObject) {
  return std::make_shared<Electron>(physicsObject);
}

inline std::shared_ptr<Track> asTrack(const std::shared_ptr<PhysicsObject> physicsObject) {
  return std::make_shared<Track>(physicsObject);
}

inline std::shared_ptr<Muon> asMuon(const std::shared_ptr<PhysicsObject> physicsObject, bool isStandalone) {
  return std::make_shared<Muon>(physicsObject, isStandalone);
}

inline std::shared_ptr<ZDCEnergy> asZDCEnergy(const std::shared_ptr<PhysicsObject> physicsObject) {
  return std::make_shared<ZDCEnergy>(physicsObject);
}

inline std::shared_ptr<LbLEvent> asLbLEvent(const std::shared_ptr<Event> physicsObject) {
  return std::make_shared<LbLEvent>(physicsObject);
}

#endif /* UserExtensionsHelpers_hpp */
