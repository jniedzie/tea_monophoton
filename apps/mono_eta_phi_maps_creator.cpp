#include "ConfigManager.hpp"
#include "CutFlowManager.hpp"
#include "EventReader.hpp"
#include "EventWriter.hpp"
#include "ExtensionsHelpers.hpp"
#include "UserExtensionsHelpers.hpp"
#include "HistogramsHandler.hpp"
#include "Profiler.hpp"
#include "HistogramsFiller.hpp"
#include "ArgsManager.hpp"
#include "LbLObjectsManager.hpp"

using namespace std;

int main(int argc, char **argv) {
  vector<string> requiredArgs = {"config"};
  vector<string> optionalArgs = {"input_path", "output_hists_path"};

  auto args = make_unique<ArgsManager>(argc, argv, requiredArgs, optionalArgs);
  ConfigManager::Initialize(args);

  auto eventReader = make_shared<EventReader>();
  auto histogramsHandler = make_shared<HistogramsHandler>();
  auto objectsManager = make_shared<LbLObjectsManager>();
  
  
  for (int iEvent = 0; iEvent < eventReader->GetNevents(); iEvent++) {
    auto event = eventReader->GetEvent(iEvent);
    objectsManager->InsertGoodPhotonsCollection(event);
    objectsManager->InsertGoodCaloTowerCollection(event);
    
    auto photonObjects = event->GetCollection("photon");    
    for (auto physicsObject : *photonObjects) {
      auto photon = asPhoton(physicsObject);
      histogramsHandler->Fill("photon_eta_phi", photon->GetEta(), photon->GetPhi());
    }

    auto goodPhotonObjects = event->GetCollection("goodPhoton");
    for (auto physicsObject : *goodPhotonObjects) {
      auto photon = asPhoton(physicsObject);
      histogramsHandler->Fill("goodPhoton_eta_phi", photon->GetEta(), photon->GetPhi());
    }

    auto caloTowerObjects = event->GetCollection("CaloTower");
    for (auto physicsObject : *caloTowerObjects) {
      auto caloTower = asCaloTower(physicsObject);
      histogramsHandler->Fill("caloTower_eta_phi", caloTower->GetEta(), caloTower->GetPhi());
    }

    auto goodCaloTowerObjects = event->GetCollection("goodCaloTower");
    for (auto physicsObject : *goodCaloTowerObjects) {
      auto caloTower = asCaloTower(physicsObject);
      histogramsHandler->Fill("goodCaloTower_eta_phi", caloTower->GetEta(), caloTower->GetPhi());
    }
  }
  
  histogramsHandler->SaveHistograms();

  auto &logger = Logger::GetInstance();
  logger.Print();

  return 0;
}