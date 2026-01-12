#include "ConfigManager.hpp"
#include "CutFlowManager.hpp"
#include "EventReader.hpp"
#include "EventWriter.hpp"
#include "ExtensionsHelpers.hpp"
#include "HistogramsHandler.hpp"
#include "Profiler.hpp"
#include "HistogramsFiller.hpp"
#include "ArgsManager.hpp"

// If you also created a histogram filler, you can include it here
// #include "MyHistogramsFiller.hpp"

using namespace std;

int main(int argc, char **argv) {
  
  auto args = make_unique<ArgsManager>(argc, argv);
  if (!args->GetString("config").has_value()) {
    fatal() << "No config file provided" << endl;
    exit(1);
  }

  ConfigManager::Initialize(args->GetString("config").value());
  auto &config = ConfigManager::GetInstance();
  
  if (args->GetString("input_path").has_value()) {
    config.SetInputPath(args->GetString("input_path").value());
  }
  if (args->GetString("output_trees_path").has_value()) {
    config.SetTreesOutputPath(args->GetString("output_trees_path").value());
  }

  auto eventReader = make_shared<EventReader>();
  auto eventWriter = make_shared<EventWriter>(eventReader);
  
  for (int iEvent = 0; iEvent < eventReader->GetNevents(); iEvent++) {
    auto event = eventReader->GetEvent(iEvent);


    // get the photon collection, and keep only the one with highest eta
    auto physicsObjects = event->GetCollection("photon");
    shared_ptr<PhysicsObject> highestEtaPhoton = nullptr;
    float maxAbsEta = -1.0;

    for (auto physicsObject : *physicsObjects) {
      float eta = fabs(physicsObject->GetAs<float>("eta"));
      if (eta > maxAbsEta) {
        maxAbsEta = eta;
        highestEtaPhoton = physicsObject;
      }
    }

    if (highestEtaPhoton != nullptr) {
      auto nowPhotonCollection = make_shared<PhysicsObjects>();
      nowPhotonCollection->push_back(highestEtaPhoton);
      event->ReplaceCollection("Photon", nowPhotonCollection);
      eventWriter->AddCurrentEvent("Events");
    }
  }
  
  eventWriter->Save();

  auto &logger = Logger::GetInstance();
  logger.Print();

  return 0;
}