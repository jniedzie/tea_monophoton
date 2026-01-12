import math
import ROOT

from Logger import info, warn, error, fatal

from lbl_params import luminosity, crossSections, nGenEvents, get_scale_factor, uncertainty_on_zero
from lbl_paths import processes, merged_histograms_path, qed_names

input_files = {}
input_aco_histograms = {}
input_mass_histograms = {}
input_cut_flow_histograms = {}


def silence_root():
  ROOT.gErrorIgnoreLevel = ROOT.kFatal


def unsilence_root():
  ROOT.gErrorIgnoreLevel = ROOT.kInfo


def load_histograms(skim, scale_to_integral=False):
  if len(input_files) > 0:
    return

  silence_root()

  for process in processes:
    file_path = merged_histograms_path.format(process, skim)

    try:
      input_files[process] = ROOT.TFile.Open(file_path)

    except OSError:
      warn(f"File not found: {file_path}")
      continue

    cut_flow_name = "cutFlow"

    input_cut_flow_histograms[process] = input_files[process].Get(cut_flow_name)

    if any([
      input_cut_flow_histograms[process] is None,
      type(input_cut_flow_histograms[process]) is ROOT.TObject
    ]):
      fatal(
        f"Some histograms for CEP normalization not found in file: {file_path}"
      )
      exit()

  if scale_to_integral:
    main_aco_hist = input_files["collisionData"].Get("diphoton_acoplanarity60")
    scale = main_aco_hist.Integral(
    ) / input_aco_histograms["collisionData"].Integral()

    for process in processes:
      input_aco_histograms[process].Scale(scale)

  unsilence_root()


def scale_non_cep_histograms():
  if not hasattr(scale_non_cep_histograms, "called"):
    scale_non_cep_histograms.called = True
  else:
    warn("scale_lbl_and_qed_histograms() called more than once!")
    return

  print("Scaling qed and lbl histograms...")

  photonScaleFactor = get_scale_factor(photon=True)[0]

  for process in processes:
    if process == "collisionData":
      continue

    if process not in input_aco_histograms:
      warn(
        f"Skipping scaling of process {process} because it was not found in the input files."
      )
      continue

    scale = luminosity * crossSections[process] * photonScaleFactor
    scale /= nGenEvents[process]

    info(f"Scaling {process} by {scale:.6f}")
    input_aco_histograms[process].Scale(scale)
    input_mass_histograms[process].Scale(scale)
    input_cut_flow_histograms[process].Scale(scale)


cep_scale = -1
cep_scale_err = -1


def get_cep_scale():
  return 294, 71


def get_alp_coupling(mass, cross_section):
  slope = 0.498809
  interstect = {
    5.0: 1.3352627780449595e-05,
    6.0: 1.459208443729578e-05,
    9.0: 1.8432339658251058e-05,
    11.0: 2.051773548082086e-05,
    14.0: 2.4387325416686648e-05,
    16.0: 2.643930742885861e-05,
    22.0: 3.323868088097294e-05,
    30.0: 4.4114845375054456e-05,
    50.0: 7.817673899208805e-05,
    90.0: 0.00019612635077675486,
    100.0: 0.000235,
  }

  return 10**(slope * math.log10(cross_section) + math.log10(interstect[mass]))


def limit_histogram(hist, max_x):
  # create a copy of the histogram, but cut at max_x. Keep bin width and min_x the same
  min_x = hist.GetXaxis().GetXmin()
  bin_width = hist.GetXaxis().GetBinWidth(1)
  n_bins = int(max_x / bin_width)
  new_hist = ROOT.TH1F(hist.GetName(), hist.GetTitle(), n_bins, min_x, max_x)
  for i in range(1, n_bins + 1):
    new_hist.SetBinContent(i, hist.GetBinContent(i))
    new_hist.SetBinError(i, hist.GetBinError(i))
  return new_hist


def add_uncertainties_on_zero(histogram):
  for i in range(1, histogram.GetNbinsX()):
    if histogram.GetBinContent(i) != 0:
      continue
    histogram.SetBinError(i, uncertainty_on_zero)

  return histogram


def sample_from_fit(hist, divide_bin_width, options=""):
  qed_sampling_n_events = 400
  qed_sampling_transition_point = 0.026
  qed_sampling_fit_max_aco = 0.2
  
  ROOT.gRandom.SetSeed(0)

  formula = f"[p0]*(x<={qed_sampling_transition_point})+(exp([p1]+[p2]*x))*(x>{qed_sampling_transition_point})"
  # formula = f"[p0]*(x<={qed_sampling_transition_point})+([p0]*exp([p2]*(x-{qed_sampling_transition_point})))*(x>{qed_sampling_transition_point})"

  if options == "superchic":
    fit = ROOT.TF1(
      "fit", "[p0]*(x<=0.032)+(exp([p1]+[p2]*x))*(x>0.032&&x<0.2)", 0, 0.1
    )
    fit.SetParameters(0.373012, 0.890578, -57.2145)
  elif options == "starlight":
    fit = ROOT.TF1(
      "fit", "[p0]*(x<=0.034)+(exp([p1]+[p2]*x))*(x>0.034&&x<0.2)", 0, 0.1
    )
    fit.SetParameters(0.265564, 0.0277834, -41.6268)
  else:
    fit = ROOT.TF1("fit", formula, 0, qed_sampling_fit_max_aco)
    fit.SetParameters(7, 3, -35)

  hist.Fit(fit, "WW")

  new_hist = hist.Clone()
  new_hist.Reset()
  new_hist.FillRandom("fit", qed_sampling_n_events)

  if divide_bin_width:
    new_hist.Scale(1.0, "width")
  new_hist.Scale(hist.Integral() / new_hist.Integral())

  canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
  # canvas.SetLogx()

  hist.SetMarkerColor(ROOT.kRed)
  hist.SetMarkerStyle(20)
  hist.SetLineColor(ROOT.kRed)
  hist.Draw("PE")

  hist.GetXaxis().SetRangeUser(1e-3, 0.2)
  hist.GetYaxis().SetRangeUser(0, 15)

  new_hist.SetMarkerColor(ROOT.kBlack)
  new_hist.SetMarkerStyle(20)
  new_hist.SetLineColor(ROOT.kBlack)
  new_hist.Draw("PEsame")

  canvas.SaveAs(f"../plots/qed_fit_{options}.pdf")

  return new_hist


def sample_lbl_from_fit(hist, start=0.06):
  ROOT.gRandom.SetSeed(0)

  fit = ROOT.TF1("fit", "[0]*exp([p1]+[p2]*x)", 0, 0.2)
  fit.SetParameters(7, 3, -70)
  hist.Fit(fit, "WW", "", 0.02, start)

  new_hist = hist.Clone()
  new_hist.Reset()
  new_hist.FillRandom("fit", 50000)

  new_hist.SetBinErrorOption(ROOT.TH1.kPoisson)

  new_hist.Scale(1.0, "width")
  new_hist.Scale(
    hist.Integral(
      hist.GetXaxis().FindFixBin(0.02),
      hist.GetXaxis().FindFixBin(start)
    ) / new_hist.Integral(
      hist.GetXaxis().FindFixBin(0.02),
      hist.GetXaxis().FindFixBin(start)
    )
  )

  final_hist = hist.Clone()

  final_hist.SetBinErrorOption(ROOT.TH1.kPoisson)

  for i in range(1, hist.GetNbinsX() + 1):
    bin_low_edge = hist.GetBinLowEdge(i)
    if bin_low_edge < start:
      continue

    final_hist.SetBinContent(i, new_hist.GetBinContent(i))
    final_hist.SetBinError(i, new_hist.GetBinError(i))

  canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
  canvas.SetLogy()

  hist.SetMarkerColor(ROOT.kRed)
  hist.SetMarkerStyle(20)
  hist.SetLineColor(ROOT.kRed)
  hist.Draw("PE")

  hist.GetXaxis().SetRangeUser(0, 0.2)
  hist.GetYaxis().SetRangeUser(1e-4, 1e3)

  final_hist.SetMarkerColor(ROOT.kBlack)
  final_hist.SetMarkerStyle(20)
  final_hist.SetLineColor(ROOT.kBlack)
  final_hist.Draw("PEsame")

  canvas.SaveAs("../plots/lbl_fit.pdf")

  return final_hist
