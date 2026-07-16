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

# draw lines at the boundaries of ECAL wedges

supermodule_width = 20  # degrees
supermodule_start = -10  # degrees, the first supermodule starts at -10 degrees in phi
n_supermodules = 18  # number of supermodules in the ECAL

lines = []

for i in range(n_supermodules):
    phi = supermodule_start + i * supermodule_width
    phi = phi * ROOT.TMath.Pi() / 180.0
    phi = ROOT.TVector2.Phi_mpi_pi(phi)
    
    lines.append(ROOT.TLine(-3, phi, 3, phi))
    lines[-1].SetLineColor(ROOT.kRed)
    lines[-1].SetLineStyle(ROOT.kDashed)
    lines[-1].Draw("same")


canvas.SaveAs("../plots/eta_phi_ratio.pdf")
