import ROOT

input_path = "../maps.root"

rebins = {
  "caloTower": 10,
  "goodCaloTower": 100,
  "photon": 100,
  "goodPhoton": 1000,
}

def main():
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  
  input_file = ROOT.TFile.Open(input_path)
  
  hists = {}
  
  hists["caloTower"] = input_file.Get("caloTower_eta_phi")
  hists["goodCaloTower"] = input_file.Get("goodCaloTower_eta_phi")
  hists["photon"] = input_file.Get("photon_eta_phi")
  hists["goodPhoton"] = input_file.Get("goodPhoton_eta_phi")
  
  canvas = ROOT.TCanvas("canvas", "canvas", 1200, 1200)
  canvas.Divide(2, 2)
  
  for i, hist in enumerate(hists.values()):
    canvas.cd(i+1)
    
    ROOT.gPad.SetLogz()
    
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("Photon #eta")
    hist.GetYaxis().SetTitle("Photon #phi")
    
  canvas.SaveAs("../plots/eta_phi_maps.pdf")
    
if __name__ == "__main__":
  main()

    