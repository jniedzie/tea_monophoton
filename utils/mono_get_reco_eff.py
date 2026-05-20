import ROOT

input_path = "../merged_skimmed_singleEG5_baseSelections_histograms.root"

rebin = 2

def main():
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  
  input_file = ROOT.TFile.Open(input_path)
  
  gen_et = input_file.Get("genPhoton_et")
  reco_et = input_file.Get("goodPhoton_et")
  
  gen_et.Sumw2()
  reco_et.Sumw2()
  
  gen_et.Rebin(rebin)
  reco_et.Rebin(rebin)
  
  gen_et.SetLineColor(ROOT.kBlue)
  reco_et.SetLineColor(ROOT.kRed)
  
  canvas = ROOT.TCanvas("canvas", "canvas", 800, 1200)
  canvas.Divide(1, 2)

  canvas.cd(1)
  ROOT.gPad.SetLogy()
  
  gen_et.Draw("")
  reco_et.Draw("same")
  
  gen_et.GetXaxis().SetTitle("Photon E_{T} (GeV)")
  gen_et.GetXaxis().SetRangeUser(0, 30)
  
  canvas.cd(2)
  
  # ratio:
  ratio = reco_et.Clone("ratio")
  ratio.Divide(reco_et, gen_et, 1, 1, "B")
  
  ratio.SetLineColor(ROOT.kBlack)
  
  ratio.Draw("pe")
  
  ratio.GetXaxis().SetTitle("Photon E_{T} (GeV)")
  ratio.GetYaxis().SetTitle("Reco efficiency")
  ratio.GetXaxis().SetRangeUser(0, 30)
  
  
  
  canvas.SaveAs("../plots/reco_efficiency_comparison.pdf")

if __name__ == "__main__":
  main()