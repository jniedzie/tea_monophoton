import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

rebin = 5

# input_path = "../test_hists.root"
# input_path = "../test_hists_onlyTrigger.root"

base_path = "/pnfs/iihe/cms/store/user/jniedzie/upc/qed_superchic"
input_path = f"{base_path}/merged_initial_singleEG5_histograms.root"

input_file = ROOT.TFile(input_path, "READ")
hist_all = input_file.Get("allGenElectron_eta_vs_phi")
hist_unmatched = input_file.Get("unmatchedGenElectron_eta_vs_phi")

hist_all.Rebin2D(rebin, rebin)
hist_unmatched.Rebin2D(rebin, rebin)

canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

ratio = hist_unmatched.Clone("ratio")

ratio.GetXaxis().SetTitle("gen electron #eta")
ratio.GetYaxis().SetTitle("gen electron #phi")

ratio.Divide(hist_all)

ratio.Draw("colz")
# hist_all.Draw("colz")
# hist_unmatched.Draw("colz")

canvas.SaveAs("../plots/eta_phi_ratio.pdf")
