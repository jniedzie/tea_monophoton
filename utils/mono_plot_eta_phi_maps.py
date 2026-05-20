import ROOT

input_path = "../maps.root"

rebins = {
  "caloTower_barrel": 1,
  "caloTower_endcap": 1,
  "goodCaloTower_barrel": 2,
  "goodCaloTower_endcap": 2,
  "photon_barrel": 5,
  "photon_endcap": 5,
  "goodPhoton_barrel": 5,
  "goodPhoton_endcap": 5,
}

def split_into_regions(hist):
  hist_barrel = hist.Clone(hist.GetName() + "_barrel")
  hist_barrel.Reset()
  
  hist_endcap_minus = hist.Clone(hist.GetName() + "_endcap_minus")
  hist_endcap_minus.Reset()
  
  hist_endcap_plus = hist.Clone(hist.GetName() + "_endcap_plus")
  hist_endcap_plus.Reset()
  
  for i in range(1, hist.GetNbinsX() + 1):
    eta = hist.GetXaxis().GetBinCenter(i)
    
    for j in range(1, hist.GetNbinsY() + 1):
      
      if abs(eta) < 1.479:
        hist_barrel.SetBinContent(i, j, hist.GetBinContent(i, j))
      elif 1.479 < eta < 3.0:
        hist_endcap_plus.SetBinContent(i, j, hist.GetBinContent(i, j))
      elif -1.479 > eta > -3.0:
        hist_endcap_minus.SetBinContent(i, j, hist.GetBinContent(i, j))
        
  return hist_barrel, hist_endcap_plus, hist_endcap_minus
  
  

def main():
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  
  input_file = ROOT.TFile.Open(input_path)
  
  input_hists = {}
  
  input_hists["caloTower"] = input_file.Get("caloTower_eta_phi")
  input_hists["goodCaloTower"] = input_file.Get("goodCaloTower_eta_phi")
  input_hists["photon"] = input_file.Get("photon_eta_phi")
  input_hists["goodPhoton"] = input_file.Get("goodPhoton_eta_phi")
  
  hists = {}
  
  for name, hist in input_hists.items():
    hist_barrel, hist_endcap_plus, hist_endcap_minus = split_into_regions(hist)
    
    hists[name + "_barrel"] = hist_barrel
    hists[name + "_endcap_plus"] = hist_endcap_plus
    hists[name + "_endcap_minus"] = hist_endcap_minus
  
  canvas = ROOT.TCanvas("canvas", "canvas", 1200, 800)
  canvas.Divide(6, 2)
  
  for i, (name, hist) in enumerate(hists.items()):
    
    canvas.cd(i+1)
    
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz()
    
    hist.Rebin2D(rebins[name.replace("_plus", "").replace("_minus", "")], rebins[name.replace("_plus", "").replace("_minus", "")])
    
    hist.SetTitle(name.replace("_", " "))
    
    hist.Draw("colz")
    hist.GetXaxis().SetTitle("Photon #eta")
    hist.GetYaxis().SetTitle("Photon #phi")
    
    hist.GetZaxis().SetRangeUser(1e-1, 1.5 * hist.GetMaximum())
    
    if "barrel" in name:
      hist.GetXaxis().SetRangeUser(-1.5, 1.5)
    elif "endcap_plus" in name:
      hist.GetXaxis().SetRangeUser(1.479, 3.0)
    elif "endcap_minus" in name:
      hist.GetXaxis().SetRangeUser(-3.0, -1.479)
    
  canvas.SaveAs("../plots/eta_phi_maps.pdf")
    
if __name__ == "__main__":
  main()

    