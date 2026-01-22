#include "ArgsManager.hpp"
#include "ConfigManager.hpp"
#include "Helpers.hpp"

using namespace std;

map<string, vector<float>*> vectors;
map<string, vector<int>*> vectorsInt;
map<string, int> sizes;
map<string, uint> uints;
map<string, ULong64_t> longs;
map<string, float[9999]> floats;
map<string, int[9999]> ints;

void SetupBranchesForMerging(TTree* tree, TTree* outputTree, map<string, vector<string>>& renameMap) {
  TObjArray* branches = tree->GetListOfBranches();
  if (!branches) {
    fatal() << "No branches found in tree " << tree->GetName() << endl;
    exit(1);
  }

  for (int i = 0; i < branches->GetEntries(); ++i) {
    auto branch = (TBranch*)branches->At(i);
    string branchName = branch->GetName();

    if (renameMap.find(branchName) == renameMap.end()) continue;

    string newBranchName = renameMap[branchName][0];

    if (renameMap[branchName][1] == "vector<float>") {
      vectors[branchName] = nullptr;
      tree->SetBranchAddress(branchName.c_str(), &vectors[branchName]);
      outputTree->Branch(renameMap[branchName][0].c_str(), &vectors[branchName]);
    } else if (renameMap[branchName][1] == "vector<int>") {
      vectorsInt[branchName] = nullptr;
      tree->SetBranchAddress(branchName.c_str(), &vectorsInt[branchName]);
      outputTree->Branch(renameMap[branchName][0].c_str(), &vectorsInt[branchName]);
    } else if (branchName == "run" || branchName == "lumis") {
      uints[branchName] = 0;
      tree->SetBranchAddress(branchName.c_str(), &uints[branchName]);
      outputTree->Branch(newBranchName.c_str(), &uints[branchName], renameMap[branchName][1].c_str());
    } else if (branchName == "event") {
      longs[branchName] = 0;
      tree->SetBranchAddress(branchName.c_str(), &longs[branchName]);
      outputTree->Branch(newBranchName.c_str(), &longs[branchName], renameMap[branchName][1].c_str());
    } else if (renameMap[branchName][1].find("]/I") != string::npos) {
      // ints[branchName] = {0}
      tree->SetBranchAddress(branchName.c_str(), &ints[branchName]);
      outputTree->Branch(newBranchName.c_str(), &ints[branchName], renameMap[branchName][1].c_str());
    } else if (renameMap[branchName][1].find("/I") != string::npos) {
      sizes[branchName] = 0;
      tree->SetBranchAddress(branchName.c_str(), &sizes[branchName]);
      outputTree->Branch(newBranchName.c_str(), &sizes[branchName], renameMap[branchName][1].c_str());
    } else if (renameMap[branchName][1].find("/F") != string::npos) {
      // floats[branchName] = {0};
      tree->SetBranchAddress(branchName.c_str(), &floats[branchName]);
      outputTree->Branch(newBranchName.c_str(), &floats[branchName], renameMap[branchName][1].c_str());
    } else {
      error() << "Not clear how to handle branch: " << branchName << " with type: " << renameMap[branchName][1] << endl;
    }
  }
}

void MergeTrees(vector<TTree*> inputTrees, TTree* outputTree, map<string, vector<string>>& renameMap, int maxEntries = -1) {
  Long64_t nEntries = maxEntries > 0 ? maxEntries : inputTrees[0]->GetEntries();

  for (auto tree : inputTrees) SetupBranchesForMerging(tree, outputTree, renameMap);

  for (Long64_t i = 0; i < nEntries; ++i) {
    if (i % 1000 == 0) info() << "Processing event " << i << " out of " << nEntries << endl;
    for (auto tree : inputTrees) tree->GetEntry(i);
    outputTree->Fill();
  }
}

int main(int argc, char** argv) {
  vector<string> requiredArgs = {"config", "input_path", "output_trees_path"};
  vector<string> optionalArgs = {};
  auto args = make_unique<ArgsManager>(argc, argv, requiredArgs, optionalArgs);
  ConfigManager::Initialize(args);

  string inputPath = args->GetString("input_path").value();
  auto inputFile = TFile::Open(inputPath.c_str());

  vector<TTree*> trees;
  vector<string> treeNames = {"ggHiNtuplizer/EventTree", "rechitanalyzerpp/zdcrechit", "l1object/L1UpgradeFlatTree", "hltanalysis/HltTree"};

  for (auto treeName : treeNames) {
    auto tree = (TTree*)inputFile->Get(treeName.c_str());
    if (!tree || tree->IsZombie()) {
      error() << "Tree " << treeName << " not found in file " << inputPath << ". It will be skipped" << endl;
      continue;
    }
    trees.push_back(tree);
  }

  // Create the output file and tree
  string outputPath = args->GetString("output_trees_path").value();
  auto outputFile = new TFile(outputPath.c_str(), "recreate");
  auto outputTree = new TTree("Events", "Events");

  // Define the branch renaming map
  auto& config = ConfigManager::GetInstance();
  map<string, vector<string>> renameMap;
  config.GetMap("branchesNames", renameMap);

  int maxEntries;
  config.GetValue("nEvents", maxEntries);

  // Process each tree
  MergeTrees(trees, outputTree, renameMap, maxEntries);

  // Write the output tree to the file
  outputTree->Write();
  outputFile->Close();

  return 0;
}
