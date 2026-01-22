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
  vector<string> requiredArgs = {"config"};
  vector<string> optionalArgs = {"input_path", "output_trees_path"};
  auto args = make_unique<ArgsManager>(argc, argv, requiredArgs, optionalArgs);
  ConfigManager::Initialize(args);
  

  auto eventReader = make_shared<EventReader>();
  auto eventWriter = make_shared<EventWriter>(eventReader);
  auto eventProcessor = make_unique<EventProcessor>();
  auto cutFlowManager = make_shared<CutFlowManager>(eventReader, eventWriter);
  auto lblSelections = make_unique<LbLSelections>();
  auto lblObjectsManager = make_unique<LbLObjectsManager>();

  cutFlowManager->RegisterCut("initial");
  cutFlowManager->RegisterCut("singlePhoton");
  cutFlowManager->RegisterCut("nElectrons");
  cutFlowManager->RegisterCut("nTracks");
  cutFlowManager->RegisterCut("nMuons");
  cutFlowManager->RegisterCut("neutralExclusivity");
  cutFlowManager->RegisterCut("ZDC");

  auto& config = ConfigManager::GetInstance();
  vector<string> eventsTreeNames;
  config.GetVector("eventsTreeNames", eventsTreeNames);

  auto& profiler = Profiler::GetInstance();
  profiler.Start("total");

  info() << "N events: " << eventReader->GetNevents() << endl;

  auto singlePhotonCutFlow = make_shared<map<string, int>>(map<string, int>{
      {"00_initial", 0},
      {"01_conversionCuts", 0},
      {"02_etCuts", 0},
      {"03_swissCross", 0},
      {"04_etaCuts", 0},
      {"05_crackCuts", 0},
      {"06_hotSpotCuts", 0},
      {"07_HEMCuts", 0},
      {"08_showerShape", 0},
      {"09_hoverE", 0},
      {"10_seedTime", 0},
  });

  for (int iEvent = 0; iEvent < eventReader->GetNevents(); iEvent++) {
    auto event = eventReader->GetEvent(iEvent);
    lblObjectsManager->InsertGoodPhotonsCollection(event, singlePhotonCutFlow);
    lblObjectsManager->InsertGoodElectronsCollection(event);
    lblObjectsManager->InsertGoodTracksCollection(event);
    lblObjectsManager->InsertGoodMuonsCollection(event);

    cutFlowManager->UpdateCutFlow("initial");
    if (!lblSelections->PassesSinglePhotonSelection(event, cutFlowManager)) continue;
    if (!lblSelections->PassesChargedExclusivity(event, cutFlowManager)) continue;
    if (!lblSelections->PassesNeutralExclusivity(event, cutFlowManager)) continue;
    if (!lblSelections->PassesZDC(event, cutFlowManager)) continue;

    for (string eventsTreeName : eventsTreeNames) {
      eventWriter->AddCurrentEvent(eventsTreeName);
    }
  }

  profiler.Stop("total");

  cutFlowManager->SaveCutFlow();
  cutFlowManager->Print();
  eventWriter->Save();

  // print single photon cut flow
  int previousCount = 0;

  info() << "Single photon cut flow:" << endl;
  for (const auto& [cutName, count] : *singlePhotonCutFlow) {
    info() << "  " << cutName << ": " << count << "\t" << count / (float)previousCount << endl;
    // info() << count << "\t" << count/(float)previousCount << endl;
    previousCount = count;
  }

  auto& logger = Logger::GetInstance();
  logger.Print();

  profiler.Print();

  return 0;
}