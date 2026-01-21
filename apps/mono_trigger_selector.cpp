#include "ArgsManager.hpp"
#include "ConfigManager.hpp"
#include "CutFlowManager.hpp"
#include "EventProcessor.hpp"
#include "EventReader.hpp"
#include "EventWriter.hpp"
#include "ExtensionsHelpers.hpp"
#include "HistogramsFiller.hpp"
#include "HistogramsHandler.hpp"
#include "LbLObjectsManager.hpp"
#include "LbLSelections.hpp"
#include "Logger.hpp"
#include "Profiler.hpp"
#include "UserExtensionsHelpers.hpp"

using namespace std;

int main(int argc, char** argv) {
  auto args = make_unique<ArgsManager>(argc, argv);

  if (!args->GetString("config").has_value()) {
    fatal() << "No config file provided" << endl;
    exit(1);
  }

  ConfigManager::Initialize(args->GetString("config").value());
  auto& config = ConfigManager::GetInstance();

  if (args->GetString("input_path").has_value()) {
    config.SetInputPath(args->GetString("input_path").value());
  }

  if (args->GetString("output_trees_path").has_value()) {
    config.SetTreesOutputPath(args->GetString("output_trees_path").value());
  }

  auto eventReader = make_shared<EventReader>();
  auto eventWriter = make_shared<EventWriter>(eventReader);
  auto eventProcessor = make_unique<EventProcessor>();
  
  vector<string> eventsTreeNames;
  config.GetVector("eventsTreeNames", eventsTreeNames);

  info() << "N events: " << eventReader->GetNevents() << endl;

  for (int iEvent = 0; iEvent < eventReader->GetNevents(); iEvent++) {
    auto event = eventReader->GetEvent(iEvent);

    if(!eventProcessor->PassesTriggerCuts(event)) continue;

    for (string eventsTreeName : eventsTreeNames) {
      eventWriter->AddCurrentEvent(eventsTreeName);
    }
  }
  eventWriter->Save();

  auto& logger = Logger::GetInstance();
  logger.Print();

  return 0;
}