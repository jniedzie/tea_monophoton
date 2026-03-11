# load the root file with a histogram and turn it into a nice PDF plot
import ROOT
from array import array

inputPath = "../plots/caloTower_etaphi.root"
# histName = "caloTower_etaphi_energy"
histName = "caloTower_etaphi_entries"


def splitHists(hist, etaThreshold=3.0):
  # create two 2D histograms, one for |eta| < etaThreshold and one for |eta| >= etaThreshold
  hist_below = hist.Clone(hist.GetName() + "_below")
  hist_above = hist.Clone(hist.GetName() + "_above")
  
  for i in range(1, hist.GetNbinsX() + 1):
    for j in range(1, hist.GetNbinsY() + 1):
      eta = hist.GetXaxis().GetBinCenter(i)
      content = hist.GetBinContent(i, j)
      
      if abs(eta) < etaThreshold:
        hist_below.SetBinContent(i, j, content)
        hist_above.SetBinContent(i, j, 0)
      else:
        hist_below.SetBinContent(i, j, 0)
        hist_above.SetBinContent(i, j, content)
        
  return hist_below, hist_above


def main():
  inputFile = ROOT.TFile(inputPath)
  hist = inputFile.Get(histName)

  histBelow, histAbove = splitHists(hist)
  histAbove.Rebin2D(5, 10)

  ROOT.gStyle.SetLineScalePS(1)

  canvas = ROOT.TCanvas("canvas", "Calo Tower Energy vs. Eta and Phi", 800, 600)
  canvas.SetRightMargin(0.15)

  histBelow.SetStats(0)  # Hide statistics box
  histBelow.Draw("COLZ")  # Draw with color palette
  histAbove.Draw("SAME COLZ")  # Draw on top of the previous histogram

  # automatically adjust the z-axis range to fit the data in both histograms
  maxBelow = histBelow.GetMaximum()
  maxAbove = histAbove.GetMaximum()
  maxZ = max(maxBelow, maxAbove)
  histBelow.SetMaximum(maxZ)
  


  # Save the canvas as a PDF
  canvas.SaveAs("../plots/caloTower_etaphi.pdf")


if __name__ == "__main__":
  main()
