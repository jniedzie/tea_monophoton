import ROOT

input_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/collisionData/merged_skimmed_singleEG5_baseSelections_noTimeCut_zdc0n0n_histograms.root"


def main():
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)

  input_file = ROOT.TFile(input_path)
  hist2D = input_file.Get("goodPhoton_Barrel_et_vs_seedTime")

  canvas = ROOT.TCanvas("canvas", "canvas", 1000, 500)
  canvas.Divide(2, 1)

  # hist2D.Rebin2D(5, 1)

  canvas.cd(1)
  hist2D.Draw("COLZ")
  hist2D.GetXaxis().SetRangeUser(0, 100)
  hist2D.GetXaxis().SetTitle("E_{T} (GeV)")
  hist2D.GetYaxis().SetTitle("t_{seed} (ns)")
  

  canvas.cd(2)

  # make profiles along y axis and plot them on the same canvas with different colors

  # profile_ranges = [(0, 6), (6, 8), (8, 10), (10, 100)]
  profile_ranges = [(4, 6), (6, 10), (10, 100)]
  # profile_ranges = [(6, 8)]
  colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan]

  legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

  profiles = {}
  Xaxis = hist2D.GetXaxis()

  for i, range in enumerate(profile_ranges):
    first_bin = Xaxis.FindBin(range[0])
    last_bin = Xaxis.FindBin(range[1] - 1e-6)  # avoid jumping to next bin

    profiles[range] = hist2D.ProjectionY(f"profile_{i}", first_bin, last_bin)

    # profiles[range] = hist2D.ProfileY(f"profile_{i}", range[0], range[1])
    rebin = 10
    
    profiles[range].Rebin(rebin)
    profiles[range].Scale(1./rebin)
    profiles[range].SetLineColor(colors[i])
    profiles[range].SetMarkerColor(colors[i])
    profiles[range].SetMarkerStyle(20)
    profiles[range].SetMarkerSize(0.7)
    profiles[range].Draw("samePE" if i > 0 else "PE")
    profiles[range].GetYaxis().SetRangeUser(1e-2, 1e7)
    profiles[range].GetXaxis().SetTitle("t_{seed} (ns)")
    legend.AddEntry(profiles[range], f"ET {range[0]}-{range[1]} GeV", "l")

  legend.Draw()

  ROOT.gPad.SetLogy()

  canvas.Update()
  canvas.SaveAs("../plots/seedTime_vs_et.pdf")


if __name__ == "__main__":
  main()
