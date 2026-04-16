import ROOT
from lbl_paths import merged_histograms_path, skim

main_rebin = 10
ratio_rebin = 1


def get_ratio(histogram, split_point):

  histogram_above = histogram.Clone("histogram_above")
  histogram_below = histogram.Clone("histogram_below")

  for i in range(1, histogram.GetNbinsX() + 1):
    x = histogram.GetXaxis().GetBinCenter(i)
    if x < split_point:
      histogram_above.SetBinContent(i, 0)
      histogram_above.SetBinError(i, 0)
    else:
      histogram_below.SetBinContent(i, 0)
      histogram_below.SetBinError(i, 0)

  histogram_below_flipped = histogram_below.Clone("histogram_below_flipped")
  for i in range(1, histogram_below.GetNbinsX() + 1):
    x = histogram_below.GetXaxis().GetBinCenter(i)
    flipped_x = 2 * split_point - x
    flipped_bin = histogram_below_flipped.GetXaxis().FindBin(flipped_x)
    if flipped_bin >= 1 and flipped_bin <= histogram_below_flipped.GetNbinsX():
      histogram_below_flipped.SetBinContent(i, histogram_below.GetBinContent(flipped_bin))
      histogram_below_flipped.SetBinError(i, histogram_below.GetBinError(flipped_bin))
    else:
      histogram_below_flipped.SetBinContent(i, 0)
      histogram_below_flipped.SetBinError(i, 0)

  # Shift histograms so that both start from zero
  histogram_above.GetXaxis().SetLimits(histogram.GetXaxis().GetXmin() - split_point, histogram.GetXaxis().GetXmax() - split_point)
  histogram_below_flipped.GetXaxis().SetLimits(histogram.GetXaxis().GetXmin() - split_point, histogram.GetXaxis().GetXmax() - split_point)

  ratio_histogram = histogram_above.Clone("ratio_histogram")
  for i in range(1, ratio_histogram.GetNbinsX() + 1):
    above_content = histogram_above.GetBinContent(i)
    below_content = histogram_below_flipped.GetBinContent(i)
    above_error = histogram_above.GetBinError(i)
    below_error = histogram_below_flipped.GetBinError(i)

    if below_content > 0:
      ratio = above_content / below_content
      if above_content > 0:
        ratio_error = ratio * ((above_error / above_content)**2 + (below_error / below_content)**2)**0.5
      else:
        ratio_error = 0
      ratio_histogram.SetBinContent(i, ratio)
      ratio_histogram.SetBinError(i, ratio_error)
    else:
      ratio_histogram.SetBinContent(i, 0)
      ratio_histogram.SetBinError(i, 0)

  return ratio_histogram, histogram_above, histogram_below_flipped


def main():
  ROOT.gROOT.SetBatch(True)
  ROOT.gStyle.SetOptStat(0)

  file = ROOT.TFile(merged_histograms_path.format("collisionData", skim))
  histogram = file.Get("goodPhoton_Barrel_seedTime")

  canvas = ROOT.TCanvas("canvas", "", 1500, 1000)
  canvas.Divide(3, 2)

  canvas.cd(1)
  ROOT.gPad.SetLogy()

  histogram.SetMarkerStyle(20)
  histogram.SetMarkerSize(1.0)
  histogram.SetMarkerColor(ROOT.kBlack)
  histogram.SetLineColor(ROOT.kBlack)

  histogram.Rebin(main_rebin)

  histogram.Draw("PE")

  histogram.SetTitle("Original distribution")
  histogram.GetXaxis().SetTitle("Seed Time (ns)")
  histogram.GetYaxis().SetTitle("Entries")
  histogram.GetXaxis().SetRangeUser(-30, 30)
  histogram.GetYaxis().SetRangeUser(1e-1, 1e3)

  canvas.cd(2)
  ROOT.gPad.SetLogy()
  plus_minus_ratio, plus_hist, minus_hist = get_ratio(histogram, 0)

  minus_hist.SetMarkerColor(ROOT.kBlue)
  plus_hist.SetMarkerColor(ROOT.kRed)

  minus_hist.Draw("PE")
  plus_hist.Draw("PE same")

  minus_hist.SetTitle("Split at 0 ns")
  minus_hist.GetXaxis().SetRangeUser(0, 30)
  minus_hist.GetYaxis().SetRangeUser(1e-1, 1e3)

  canvas.cd(5)

  plus_minus_ratio.Rebin(ratio_rebin)
  plus_minus_ratio.Scale(1 / ratio_rebin)

  plus_minus_ratio.Draw("PE")
  plus_minus_ratio.SetTitle("Ratio of entries above 0 to below 0")
  plus_minus_ratio.GetXaxis().SetTitle("Time from 0 (ns)")
  plus_minus_ratio.GetYaxis().SetTitle("Ratio of entries above 0 to below 0")
  plus_minus_ratio.GetXaxis().SetRangeUser(0, 20)
  plus_minus_ratio.GetYaxis().SetRangeUser(0, 3)

  line = ROOT.TLine(0, 1, 20, 1)
  line.SetLineColor(ROOT.kRed)
  line.SetLineStyle(ROOT.kDashed)
  line.Draw("same")

  canvas.cd(3)
  ROOT.gPad.SetLogy()
  # find the x value where the original histogram has the maximum:
  max_bin = histogram.GetMaximumBin()
  max_x = histogram.GetXaxis().GetBinUpEdge(max_bin)
  print(f"Max bin: {max_bin}, Max x: {max_x}")
  peak_ratio, above_hist, below_hist = get_ratio(histogram, max_x)

  below_hist.SetMarkerColor(ROOT.kBlue)
  above_hist.SetMarkerColor(ROOT.kRed)

  below_hist.Draw("PE")
  above_hist.Draw("PE same")

  below_hist.SetTitle(f"Split at {max_x:.1f} ns")
  below_hist.GetXaxis().SetRangeUser(0, 30)
  below_hist.GetYaxis().SetRangeUser(1e-1, 1e3)

  canvas.cd(6)
  peak_ratio.Rebin(ratio_rebin)
  peak_ratio.Scale(1 / ratio_rebin)

  peak_ratio.Draw("PE")
  peak_ratio.SetTitle(f"Ratio of entries above {max_x:.1f} ns to below {max_x:.1f} ns")
  peak_ratio.GetXaxis().SetTitle("Time from peak (ns)")
  peak_ratio.GetYaxis().SetTitle("Ratio of entries above peak to below peak")
  peak_ratio.GetXaxis().SetRangeUser(0, 20)
  peak_ratio.GetYaxis().SetRangeUser(0, 20)

  line.Draw("same")

  canvas.cd(1)
  # draw vertical lines at 0 and max_x
  line_zero = ROOT.TLine(0, 1e-1, 0, 1e3)
  line_zero.SetLineColor(ROOT.kRed)
  line_zero.SetLineStyle(ROOT.kDashed)
  line_zero.Draw("same")
  line_peak = ROOT.TLine(max_x, 1e-1, max_x, 1e3)
  line_peak.SetLineColor(ROOT.kRed)
  line_peak.SetLineStyle(ROOT.kDashed)
  line_peak.Draw("same")

  canvas.SaveAs("../plots/mono_seedTime.pdf")


if __name__ == "__main__":
  main()
