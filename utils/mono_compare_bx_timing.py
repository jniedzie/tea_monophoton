import ROOT

from lbl_paths import base_path, skim

file_path = f"{base_path}/collisionData/merged_{skim}_histograms.root"


def main():
  input_file = ROOT.TFile.Open(file_path)

  canvas = ROOT.TCanvas("canvas", "canvas", 1500, 500)
  canvas.Divide(3, 1)

  legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

  for i, prefix in enumerate(["", "Barrel_", "EndCap_"]):

    hist_after_collision = input_file.Get(f"goodPhoton/afterCollisionBX_{prefix}seedTime")
    hist_without_collision = input_file.Get(f"goodPhoton/withoutCollisionBX_{prefix}seedTime")

    hist_after_collision.SetLineColor(ROOT.kRed)
    hist_without_collision.SetLineColor(ROOT.kBlue)

    if i == 0:
      legend.AddEntry(hist_after_collision, "afterCollisionBX", "l")
      legend.AddEntry(hist_without_collision, "withoutCollisionBX", "l")

    canvas.cd(i + 1)
    hist_after_collision.Draw("hist")
    hist_without_collision.Draw("hist same")
    legend.Draw()

  canvas.SaveAs("../plots/mono_compare_bx_timing.pdf")


if __name__ == "__main__":
  main()
