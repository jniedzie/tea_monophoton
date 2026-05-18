#include "CutFlowManager.hpp"
#include "Event.hpp"
#include "LbLObjectsManager.hpp"
#include "Helpers.hpp"

class LbLSelections {
 public:
  LbLSelections();
  ~LbLSelections() = default;

  bool PassesZeroPhotonAndElectronSelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesNeutralExclusivity(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesDiphotonSelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr, std::shared_ptr<std::map<std::string, int>> cutFlow=nullptr);
  bool PassesSinglePhotonSelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesThreePhotonsSelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesDielectronSelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr, bool sameCharge=false);
  bool PassesChargedExclusivity(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesDiphotonPt(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesBeamHaloFilters(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesZDC(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesTracksPlusPhotonsSelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);

  bool PassesDielectronChargedExclusivity(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);
  bool PassesAcoplanaritySelection(std::shared_ptr<Event> event, std::shared_ptr<CutFlowManager> cutFlowManager=nullptr);

private:
  std::map<std::string, float> eventCuts;
  std::shared_ptr<LbLObjectsManager> lblObjectsManager;
};