#include "ArgsManager.hpp"
#include "ConfigManager.hpp"
#include "CutFlowManager.hpp"
#include "EventReader.hpp"
#include "ExtensionsHelpers.hpp"
#include "HistogramsFiller.hpp"
#include "HistogramsHandler.hpp"
#include "MonoHistogramsFiller.hpp"
#include "MonoObjectsManager.hpp"
#include "Logger.hpp"

using namespace std;

int main(int argc, char** argv) {
  vector<string> requiredArgs = {"config"};
  vector<string> optionalArgs = {"input_path", "output_hists_path"};
  auto args = make_unique<ArgsManager>(argc, argv, requiredArgs, optionalArgs);
  ConfigManager::Initialize(args);

  info() << "Creating objects" << endl;
  auto eventReader = make_shared<EventReader>();
  auto histogramsHandler = make_shared<HistogramsHandler>();
  auto cutFlowManager = make_shared<CutFlowManager>(eventReader);
  auto histogramsFiller = make_unique<HistogramsFiller>(histogramsHandler);
  auto monoHistogramsFiller = make_unique<MonoHistogramsFiller>(histogramsHandler);
  auto monoObjectsManager = make_unique<MonoObjectsManager>();

  cutFlowManager->RegisterCut("initial");

  info() << "Starting event loop" << endl;
  for (int iEvent = 0; iEvent < eventReader->GetNevents(); iEvent++) {
    auto event = eventReader->GetEvent(iEvent);

    monoObjectsManager->InsertGoodPhotonsCollection(event);
    monoObjectsManager->InsertGoodElectronsCollection(event);

    try {
      monoObjectsManager->InsertGenPhotonsCollection(event);
      monoObjectsManager->InsertGenElectronsCollection(event);
    } catch (const Exception& e) {
      warn() << "No gen-level information found in event. Skipping gen-level histograms filling." << endl;
    }

    cutFlowManager->UpdateCutFlow("initial");
    histogramsFiller->FillDefaultVariables(event);
    monoHistogramsFiller->Fill(event);
  }

  info() << "Finishing up" << endl;
  cutFlowManager->Print();
  histogramsFiller->FillCutFlow(cutFlowManager);

  info() << "Saving histograms" << endl;
  histogramsHandler->SaveHistograms();

  info() << "Logger:" << endl;
  auto& logger = Logger::GetInstance();
  logger.Print();
  return 0;
}