from lbl_paths import base_path, trigger, merged_histograms_path
from lbl_plotter_config import histograms
from Logger import info, warn

import ROOT

skim = "skimmed_singleEG5_baseSelections"
input_path = merged_histograms_path.format("collisionData", skim)

histograms_to_plot = [
  # ("et", True, 2, 0, 30),
  ("eta", False, 2, -3, 3),
  # ("phi", False, 5, -4, 4),
  # ("seedTime", True, 5, -3, 3),
  # ("SCEtaWidth", True, 5, 0, 0.03),
  # ("SCPhiWidth", True, 20, 0, 0.15),
  # ("sigmaIEtaIEta2012", True, 2, 0, 0.06),
]


options = [
  ("", ROOT.kBlack, ROOT.kSolid),
  # ("_noMuonSegments", ROOT.kGreen+2, ROOT.kSolid),
  # ("_noDTsegments", ROOT.kBlue, ROOT.kSolid),
  # ("_noDTcosmicSegments", ROOT.kRed, ROOT.kDashed),
  ("_noCSCsegments", ROOT.kRed, ROOT.kDashed),
]

def main():
  ROOT.gStyle.SetOptStat(0)
  ROOT.gROOT.SetBatch(True)
  
  info(f"Opening file: {input_path}")
  
  input_file = ROOT.TFile.Open(input_path, "READ")

  canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
  # canvas.Divide(3, 3)

  hists = {}


  legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)

  for i, (hist_name, logy, rebin, min_x, max_x) in enumerate(histograms_to_plot):

    canvas.cd(i + 1)

    for (option, color, line_style) in options:
      key = f"{hist_name}_{option}"
      name = f"goodPhoton{option}_{hist_name}"
      hists[key] = input_file.Get(name)
      
      if hists[key] is None or type(hists[key]) == ROOT.TObject:
        warn(f"Warning: histogram {name} not found in file.")
        continue

      hists[key].Sumw2()
      
      hists[key].SetLineColor(color)
      hists[key].SetMarkerColor(color)
      hists[key].SetMarkerStyle(20)
      hists[key].SetLineStyle(line_style)
      hists[key].SetLineWidth(2)
      hists[key].Rebin(rebin)
    
      hists[key].GetXaxis().SetRangeUser(min_x, max_x)
      hists[key].Draw("" if option == "" else "same")
      
      if i == 0:
        legend.AddEntry(hists[key], option if option != "" else "with all segments", "l")

    ROOT.gPad.SetLogy(logy)
    legend.Draw()
    

  canvas.SaveAs(
    f"../plots/{skim}/compare_muon_segments.pdf".replace("skimmed_", "")
  )


if __name__ == "__main__":
  main()
