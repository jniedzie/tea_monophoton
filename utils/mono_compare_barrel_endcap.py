from lbl_paths import base_path, trigger, merged_histograms_path
from lbl_plotter_config import histograms

import ROOT

# skim = f"skimmed_{trigger}_endCapTests"
skim = f"skimmed_{trigger}_lblSelections_et-gt4_widthCuts_shapeCuts_deltaT_tightSigmaIEtaIEta"

input_path = merged_histograms_path.format("collisionData", skim)

histograms_to_plot = [
  ("et", True, 2, 0, 100),
  ("eta", True, 1, -3, 3),
  ("phi", True, 1, -4, 4),
  ("seedTime", True, 1, -3, 3),
  ("SCEtaWidth", True, 2, 0, 0.04),
  ("SCPhiWidth", True, 10, 0, 0.2),
  ("sigmaIEtaIEta2012", True, 1, 0, 0.08),
]


def main():
  input_file = ROOT.TFile.Open(input_path, "READ")

  canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
  canvas.Divide(3, 3)

  hists = {}

  for i, (hist_name, logy, rebin, min_x, max_x) in enumerate(histograms_to_plot):

    hists[f"{hist_name}_barrel"] = input_file.Get(
      f"goodPhoton_Barrel_{hist_name}"
    )
    hists[f"{hist_name}_endcap"] = input_file.Get(
      f"goodPhoton_EndCap_{hist_name}"
    )

    hists[f"{hist_name}_barrel"].SetLineColor(ROOT.kBlue)
    hists[f"{hist_name}_endcap"].SetLineColor(ROOT.kRed)

    hists[f"{hist_name}_barrel"].Rebin(rebin)
    hists[f"{hist_name}_endcap"].Rebin(rebin)

    canvas.cd(i + 1)
    hists[f"{hist_name}_barrel"].GetXaxis().SetRangeUser(min_x, max_x)
    
    hists[f"{hist_name}_barrel"].DrawNormalized("HIST")
    
    
    
    hists[f"{hist_name}_endcap"].DrawNormalized("HIST SAME")

    ROOT.gPad.SetLogy(logy)

  canvas.SaveAs(
    f"../plots/{skim}/compare_barrel_endcap.pdf".replace("skimmed_", "")
  )


if __name__ == "__main__":
  main()
