#include "CutFlowManager.hpp"
#include "Event.hpp"
#include "LbLObjectsManager.hpp"
#include "Helpers.hpp"

class LbLSelections {
 public:
  LbLSelections();
  ~LbLSelections() = default;

  bool PassesNeutralExclusivity(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesPhotonSelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesChargedExclusivity(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesBeamHaloFilters(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesZDC(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);

private:
  std::map<std::string, float> eventCuts;
  std::shared_ptr<LbLObjectsManager> lblObjectsManager;
};