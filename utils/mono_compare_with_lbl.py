import ROOT

input_path_mono = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/collisionData/merged_skimmed_singleEG5_baseSelections_haloFlags_histograms.root"
input_path_lbl = "/eos/cms/store/cmst3/group/lightbylight/tea_samples/collisionData/merged_skimmed_lblSelections_final_photonMatchingdeltaPhi0p15_test_noTimingCut_histograms.root"

hist_name_mono = "goodPhoton_withoutStandaloneMuon_seedTime"
hist_name_lbl = "goodPhotonSR_seedTime"

mono_rebin = 10

def prepare_hist(hist, color):
  hist.SetLineColor(color)
  hist.SetMarkerColor(color)
  hist.SetMarkerStyle(20)
  hist.SetMarkerSize(1.0)
  hist.Sumw2()

def main():
  ROOT.gStyle.SetOptStat(0)
  
  input_file_mono = ROOT.TFile(input_path_mono)
  input_file_lbl = ROOT.TFile(input_path_lbl)
  
  hist_mono = input_file_mono.Get(hist_name_mono)
  hist_lbl = input_file_lbl.Get(hist_name_lbl)
  
  prepare_hist(hist_mono, ROOT.kRed)
  prepare_hist(hist_lbl, ROOT.kBlue)
  
  hist_mono.Rebin(mono_rebin)
  
  canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
  legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
  
  ROOT.gPad.SetLogy()
  
  # hist_lbl.Scale(1.0 / hist_lbl.Integral())
  # hist_mono.Scale(1.0 / hist_mono.Integral())
  
  hist_mono.Scale(hist_lbl.GetBinContent(hist_lbl.GetXaxis().FindBin(0)) / hist_mono.GetBinContent(hist_mono.GetXaxis().FindBin(0)))
  
  hist_lbl.Draw("PE")
  hist_mono.Draw("PE same")
  
  hist_lbl.SetTitle("mono-#gamma scaled to match LbL at seed time = 0 ns")
  
  hist_lbl.GetXaxis().SetRangeUser(-30, 30)
  hist_lbl.GetXaxis().SetTitle("Photon seed time (ns)")
  # hist_lbl.GetYaxis().SetRangeUser(1e-3, 1)
  
  fit_fun_lbl = ROOT.TF1("fit_fun_lbl", "gaus(0)", -30, 30)
  fit_fun_lbl.SetParameters(hist_lbl.GetMaximum(), hist_lbl.GetMean(), hist_lbl.GetRMS())
  fit_fun_lbl.SetLineColor(ROOT.kBlue)
  
  hist_lbl.Fit(fit_fun_lbl, "R0")
  fit_fun_lbl.SetNpx(1000)
  fit_fun_lbl.Draw("same")
  
  fit_fun_mono = ROOT.TF1("fit_fun_mono", "gaus(0)+gaus(3)+gaus(6)", -30, 30)
  
  # set the first 3 parameters to the values from the LbL fit
  fit_fun_mono.SetParameter(0, fit_fun_lbl.GetParameter(0))  # amplitude
  fit_fun_mono.SetParameter(1, fit_fun_lbl.GetParameter(1))  # mean
  fit_fun_mono.SetParameter(2, fit_fun_lbl.GetParameter(2))  # sigma
  
  fit_fun_mono.SetParameter(3, 1)
  fit_fun_mono.SetParameter(4, 0)
  fit_fun_mono.SetParameter(5, 5)
  
  fit_fun_mono.SetParameter(6, 1)
  fit_fun_mono.SetParameter(7, 0)
  fit_fun_mono.SetParameter(8, 3)
  
  fit_fun_mono.SetLineColor(ROOT.kRed)
  hist_mono.Fit(fit_fun_mono, "R0")
  fit_fun_mono.SetNpx(1000)
  fit_fun_mono.Draw("same")
  
  mono_central_gaus = ROOT.TF1("mono_central_gaus", "gaus(0)", -30, 30)
  mono_central_gaus.SetParameters(fit_fun_mono.GetParameter(0), fit_fun_mono.GetParameter(1), fit_fun_mono.GetParameter(2))
  
  # scale central gausian to match the maximum of the complete mono fit
  scale = fit_fun_mono.GetMaximum() / mono_central_gaus.GetMaximum()
  mono_central_gaus.SetParameter(0, mono_central_gaus.GetParameter(0) * scale)
  
  mono_central_gaus.SetLineColor(ROOT.kRed)
  mono_central_gaus.SetLineStyle(ROOT.kDashed)
  mono_central_gaus.SetNpx(1000)
  mono_central_gaus.Draw("same")

  legend.AddEntry(hist_mono, "mono-#gamma", "l")
  legend.AddEntry(hist_lbl, "LbL", "l")
  legend.Draw()
  
  canvas.SaveAs("../plots/mono_vs_lbl_comparison.pdf")

if __name__ == "__main__":
  main()
