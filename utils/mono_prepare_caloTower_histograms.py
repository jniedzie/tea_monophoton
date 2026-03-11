import ROOT

inputPath = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/collisionData/initial_singleEG5/ntuple_0.root"
outputPath = "../plots/caloTower_etaphi.root"

maxEvents = 50000

nBinsEta = 500
nBinsPhi = 500


def main():
  inputFile = ROOT.TFile(inputPath)
  tree = inputFile.Get("Events")

  hist_entries = ROOT.TH2F("caloTower_etaphi_entries", "Calo Tower Energy vs. Eta and Phi;Eta;Phi",
                           nBinsEta, -5, 5, nBinsPhi, -3.14, 3.14)
  hist_energy = ROOT.TH2F("caloTower_etaphi_energy", "Calo Tower Energy vs. Eta and Phi;Eta;Phi",
                          nBinsEta, -5, 5, nBinsPhi, -3.14, 3.14)

  ROOT.gInterpreter.Declare("""
  void fillHistograms(TTree* tree, TH2F* hist_entries, TH2F* hist_energy, int maxEvents) {
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t iEvent = 0; iEvent < nEntries && iEvent < maxEvents; ++iEvent) {
      if (iEvent % 1000 == 0) {
        printf("Processing event %lld...\\n", iEvent);
      }

      static std::vector<float>* eta = nullptr;
      static std::vector<float>* phi = nullptr;
      static std::vector<float>* energy = nullptr;

      if (iEvent == 0) {
        tree->SetBranchAddress("CaloTower_eta", &eta);
        tree->SetBranchAddress("CaloTower_phi", &phi);
        tree->SetBranchAddress("CaloTower_energy", &energy);
      }

      tree->GetEntry(iEvent);

      if (eta && phi && energy) {
        auto nCaloTower = eta->size();  // Use the size of the vector directly

        for (size_t i = 0; i < nCaloTower; ++i) {
          hist_entries->Fill((*eta)[i], (*phi)[i]);
          hist_energy->Fill((*eta)[i], (*phi)[i], (*energy)[i]);
        }
      }
    }
  }
  """)

  fillHistograms = ROOT.fillHistograms
  fillHistograms(tree, hist_entries, hist_energy, maxEvents)

  outputFile = ROOT.TFile(outputPath, "RECREATE")
  hist_entries.Write()
  hist_energy.Write()
  outputFile.Close()


if __name__ == "__main__":
  main()
