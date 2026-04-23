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

#include "TClass.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TTreeFormula.h"

using namespace std;

namespace {
tuple<string, string> SplitTreePath(const string& treePath) {
  auto pos = treePath.find_last_of('/');
  if (pos == string::npos) return {"", treePath};
  return {treePath.substr(0, pos), treePath.substr(pos + 1)};
}

TDirectory* GetOrCreateDirectory(TFile* file, const string& directoryPath) {
  if (directoryPath.empty()) return file;

  TDirectory* currentDir = file;
  for (const auto& part : split(directoryPath, '/')) {
    if (part.empty()) continue;

    auto* nextDir = currentDir->GetDirectory(part.c_str());
    if (!nextDir) nextDir = currentDir->mkdir(part.c_str());
    currentDir = nextDir;
  }
  return currentDir;
}

void CopyTreeWithSelection(TTree* inputTree, TDirectory* outputDir, const vector<long long>& selectedEntries) {
  if (!inputTree || !outputDir) return;

  outputDir->cd();
  auto* outputTree = inputTree->CloneTree(0);
  outputTree->SetName(inputTree->GetName());
  outputTree->SetTitle(inputTree->GetTitle());

  for (const auto entry : selectedEntries) {
    if (inputTree->GetEntry(entry) < 0) {
      fatal() << "Couldn't read entry " << entry << " from tree " << inputTree->GetName() << endl;
      continue;
    }
    outputTree->Fill();
  }

  outputTree->Write();
}

void CopyDirectoryContents(TDirectory* inputDir,
                           TDirectory* outputDir,
                           const vector<long long>& selectedEntries,
                           const string& directoryPath = "") {
  TIter nextKey(inputDir->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = dynamic_cast<TKey*>(nextKey()))) {
    TClass* objectClass = gROOT->GetClass(key->GetClassName());
    if (!objectClass) continue;

    if (objectClass->InheritsFrom(TDirectory::Class())) {
      auto* inputSubdir = dynamic_cast<TDirectory*>(inputDir->Get(key->GetName()));
      if (!inputSubdir) continue;
      auto* outputSubdir = outputDir->mkdir(inputSubdir->GetName(), inputSubdir->GetTitle());
      const string subdirPath = directoryPath.empty() ? inputSubdir->GetName() : directoryPath + "/" + inputSubdir->GetName();
      CopyDirectoryContents(inputSubdir, outputSubdir, selectedEntries, subdirPath);
      continue;
    }

    if (objectClass->InheritsFrom(TTree::Class())) {
      auto* inputTree = dynamic_cast<TTree*>(inputDir->Get(key->GetName()));
      if (!inputTree) continue;
      outputDir->cd();
      CopyTreeWithSelection(inputTree, outputDir, selectedEntries);
      continue;
    }

    auto* object = key->ReadObj();
    if (!object) continue;
    outputDir->cd();
    object->Write(object->GetName());
  }
}

int RunDirectOnHIForest(const vector<string>& triggerNames) {
  auto& config = ConfigManager::GetInstance();

  string inputFilePath, outputFilePath;
  int maxEvents = -1;
  config.GetValue("inputFilePath", inputFilePath);
  config.GetValue("treeOutputFilePath", outputFilePath);
  config.GetValue("nEvents", maxEvents);

  const string treePath = "hltanalysis/HltTree";
  info() << "Detected HIForest input. Running trigger selector directly on " << treePath << endl;
  info() << "Input file path: " << inputFilePath << endl;

  gSystem->RedirectOutput("/dev/null", "a");
  auto* inputFile = TFile::Open(inputFilePath.c_str());
  gSystem->RedirectOutput(0);
  if (!inputFile || inputFile->IsZombie()) {
    fatal() << "Couldn't open input file: " << inputFilePath << endl;
    return 1;
  }

  auto* inputTree = dynamic_cast<TTree*>(inputFile->Get(treePath.c_str()));
  if (!inputTree) {
    fatal() << "Couldn't load HIForest tree: " << treePath << endl;
    return 1;
  }

  vector<unique_ptr<TTreeFormula>> triggerFormulas;
  for (const auto& triggerName : triggerNames) {
    if (!inputTree->GetBranch(triggerName.c_str())) {
      warn() << "Trigger not present in HIForest tree: " << triggerName << endl;
      continue;
    }

    auto formula = make_unique<TTreeFormula>(triggerName.c_str(), triggerName.c_str(), inputTree);
    if (formula->GetNdim() <= 0) {
      warn() << "Couldn't evaluate trigger formula for branch: " << triggerName << endl;
      continue;
    }
    triggerFormulas.push_back(std::move(formula));
  }

  if (triggerFormulas.empty()) {
    fatal() << "None of the requested triggers are available in " << treePath << endl;
    return 1;
  }

  makeParentDirectories(outputFilePath);
  auto* outputFile = TFile::Open(outputFilePath.c_str(), "recreate");
  if (!outputFile || outputFile->IsZombie()) {
    fatal() << "Couldn't create output file: " << outputFilePath << endl;
    return 1;
  }

  const auto nEntries = inputTree->GetEntries();
  const long long nEvents = (maxEvents >= 0 && maxEvents < nEntries) ? maxEvents : nEntries;
  info() << "N events: " << nEvents << endl;

  vector<long long> selectedEntries;
  selectedEntries.reserve(nEvents);
  long long nPassed = 0;
  for (long long iEvent = 0; iEvent < nEvents; ++iEvent) {
    inputTree->GetEntry(iEvent);

    bool passesTrigger = false;
    for (auto& formula : triggerFormulas) {
      formula->GetNdata();
      if (formula->EvalInstance() != 0) {
        passesTrigger = true;
        break;
      }
    }
    if (!passesTrigger) continue;

    selectedEntries.push_back(iEvent);
    ++nPassed;
  }

  CopyDirectoryContents(inputFile, outputFile, selectedEntries);

  outputFile->Close();
  inputFile->Close();

  info() << "Saved " << nPassed << " / " << nEvents << " events to " << outputFilePath << endl;
  return 0;
}
}  // namespace

int main(int argc, char** argv) {
  vector<string> requiredArgs = {"config"}; 
  vector<string> optionalArgs = {"input_path", "output_trees_path"};
  auto args = make_unique<ArgsManager>(argc, argv, requiredArgs, optionalArgs);
  ConfigManager::Initialize(args);

  auto& config = ConfigManager::GetInstance();
  vector<string> eventsTreeNames;
  vector<string> triggerNames;
  config.GetVector("eventsTreeNames", eventsTreeNames);
  config.GetVector("triggerSelection", triggerNames);

  if (eventsTreeNames.size() == 1 && eventsTreeNames[0] == "hltanalysis/HltTree") {
    auto status = RunDirectOnHIForest(triggerNames);
    auto& logger = Logger::GetInstance();
    logger.Print();
    return status;
  }

  auto eventReader = make_shared<EventReader>();
  auto eventWriter = make_shared<EventWriter>(eventReader);
  auto eventProcessor = make_unique<EventProcessor>();

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
