#ifndef LbLEvent_hpp
#define LbLEvent_hpp

#include <memory>
#include <optional>
#include <string>

#include "Event.hpp"
#include "Helpers.hpp"

class LbLEvent {
 public:
  LbLEvent(std::shared_ptr<Event> event_) : event(event_) {}

  auto Get(std::string branchName) { return event->Get(branchName); }
  
  template <typename T>
  T GetAs(std::string branchName) { return event->GetAs<T>(branchName); }
  std::shared_ptr<PhysicsObjects> GetCollection(std::string name) const { return event->GetCollection(name); }
  void AddExtraCollections() { event->AddExtraCollections(); }

  float GetDeltaEt();
  float GetCosThetaStar(bool doElectrons=false);
  float GetDiphotonAcoplanarity();
  std::vector<std::shared_ptr<PhysicsObject>> GetGenPhotons();
  std::vector<std::shared_ptr<PhysicsObject>> GetGenMatchedRecoPhotons();
  std::optional<bool> HasCollisionInPreviousBXs(int nBXs);

  bool IsData() { return GetAs<int>("runNumber") != 1; }

 private:
  struct MonoCollisionCache;
  static const MonoCollisionCache& GetMonoCollisionCache();
  static std::string GetMonoUtilsDir();
  static MonoCollisionCache LoadMonoCollisionCache();

  std::shared_ptr<Event> event;
};

#endif /* LbLEvent_hpp */
