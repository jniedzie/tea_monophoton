#ifndef LbLHistogramsFiller_hpp
#define LbLHistogramsFiller_hpp

#include "Event.hpp"
#include "EventProcessor.hpp"
#include "Helpers.hpp"
#include "HistogramsHandler.hpp"
#include "UserExtensionsHelpers.hpp"

class LbLHistogramsFiller {
 public:
  LbLHistogramsFiller(std::shared_ptr<HistogramsHandler> histogramsHandler_);
  ~LbLHistogramsFiller();

  void Fill(const std::shared_ptr<Event> event);

 private:
  std::shared_ptr<HistogramsHandler> histogramsHandler;
  std::unique_ptr<EventProcessor> eventProcessor;

  void FillMonoPhotonHistograms(const std::shared_ptr<Event> event);
  void FillMonoPhotonHistograms(const std::shared_ptr<Photon> photon, std::string prefix="");
  void FillEGammaHistograms(const std::shared_ptr<Event> event);
  void FillGenLevelHistograms(const std::shared_ptr<Event> event);
  void FillEventLevelHistograms(const std::shared_ptr<Event> event);
  void SaveHighEtPhotonsInfo(const std::shared_ptr<Event> event);

  std::map<std::string, float> dataBlinding;
};

#endif /* LbLHistogramsFiller_hpp */
