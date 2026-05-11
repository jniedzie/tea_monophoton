import ROOT

from Logger import info
from lbl_paths import base_path, skim

file_path = f"{base_path}/collisionData/merged_{skim}_histograms.root"

rebin = 20


def main():
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)
  
  info(f"Opening file: {file_path}")
  input_file = ROOT.TFile.Open(file_path)

  canvas = ROOT.TCanvas("canvas", "canvas", 500, 1500)
  canvas.Divide(1, 3)

  legend = ROOT.TLegend(0.2, 0.7, 0.9, 0.9)

  for i, prefix in enumerate(["", "Barrel_", "EndCap_"]):

    hist_after_collision = input_file.Get(f"goodPhoton_afterCollisionBX_{prefix}seedTime")
    hist_without_collision = input_file.Get(f"goodPhoton_withoutCollisionBX_{prefix}seedTime")

    hist_after_collision.SetTitle(f"{prefix[:-1]} photons")
    hist_after_collision.GetXaxis().SetTitle("Seed time (ns)")

    hist_after_collision.SetLineColor(ROOT.kRed)
    hist_after_collision.SetMarkerColor(ROOT.kRed)
    hist_after_collision.SetMarkerStyle(20)
    hist_after_collision.SetMarkerSize(1.0)
    
    hist_without_collision.SetLineColor(ROOT.kBlue)
    hist_without_collision.SetMarkerColor(ROOT.kBlue)
    hist_without_collision.SetMarkerStyle(20)
    hist_without_collision.SetMarkerSize(1.0)
    
    hist_after_collision.Rebin(rebin)
    hist_without_collision.Rebin(rebin)
    
    try:
      hist_after_collision.Scale(1.0 / hist_after_collision.Integral())
    except ZeroDivisionError:
      info(f"Warning: Integral of hist_after_collision is zero, skipping scaling.")
    
    try:
      hist_without_collision.Scale(1.0 / hist_without_collision.Integral())
    except ZeroDivisionError:
      info(f"Warning: Integral of hist_without_collision is zero, skipping scaling.")

    if i == 0:
      legend.AddEntry(hist_after_collision, "Had collision in previous N BXs'", "l")
      legend.AddEntry(hist_without_collision, "Didn't have collision in previous N BXs'", "l")

    canvas.cd(i + 1)
    ROOT.gPad.SetLogy()

    hist_after_collision.Draw("pe")
    hist_without_collision.Draw("pe same")
    
    hist_after_collision.GetXaxis().SetRangeUser(-30, 30)
    hist_after_collision.GetYaxis().SetRangeUser(1e-4, 30)
    
    legend.Draw()

  canvas.SaveAs("../plots/mono_compare_bx_timing.pdf")


if __name__ == "__main__":
  main()
