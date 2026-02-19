#pragma once

#include "CutFlowManager.hpp"
#include "Event.hpp"
#include "Helpers.hpp"
#include "UserExtensionsHelpers.hpp"

class LbLObjectsManager {
 public:
  LbLObjectsManager();
  ~LbLObjectsManager() = default;

  void InsertGoodPhotonsCollection(std::shared_ptr<Event> event, std::shared_ptr<std::map<std::string, int>> cutFlow = nullptr);
  void InsertGoodElectronsCollection(std::shared_ptr<Event> event);
  void InsertGoodTracksCollection(std::shared_ptr<Event> event);
  void InsertGoodMuonsCollection(std::shared_ptr<Event> event);

  void InsertGenPhotonsCollection(std::shared_ptr<Event> event);
  void InsertGenElectronsCollection(std::shared_ptr<Event> event);

  bool IsGoodPhoton(const std::shared_ptr<Photon> photon, std::shared_ptr<std::map<std::string, int>> cutFlow = nullptr);

 private:
  bool IsGoodElectron(const std::shared_ptr<Electron> electron);
  bool IsGoodTrack(const std::shared_ptr<Track> track);
  bool IsGoodMuon(const std::shared_ptr<Muon> track);

  std::shared_ptr<PhysicsObjects> GetGenParticles(const std::shared_ptr<Event> event, int pid);
  int GetParticlePid(const std::shared_ptr<PhysicsObject> particle);

  std::map<std::string, float> detectorParams, caloEtaEdges;
  std::vector<int> knownPids;
};
