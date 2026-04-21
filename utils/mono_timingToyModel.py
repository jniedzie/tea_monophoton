import ROOT
from math import exp, tan, atan, sin, log

ecal_inner_radius = 1.29  # ECAL inner radius in meters
c = 3e8  # speed of light in m/s

ecal_half_length = 3.0  # half-length of the ECAL in meters

n_events = 10000


def sgn(x):
  return (x > 0) - (x < 0)


def get_random_float(min_value, max_value):
  return ROOT.gRandom.Uniform(min_value, max_value)


def get_timing_correction_function():
  correctionFunction = ROOT.TF1("correctionFunction", "sqrt([0]*[0]+x*x)/[1] * 1e9", -15, 15)
  correctionFunction.SetParameter(0, ecal_inner_radius)  # ECAL inner radius (l)
  correctionFunction.SetParameter(1, c)  # speed of light (c)
  correctionFunction.SetLineColor(ROOT.kBlue)
  correctionFunction.SetTitle("Timing Correction Function;z_{ECAL} (m);Timing Correction (ns)")

  return correctionFunction


def main():
  ROOT.gROOT.SetBatch(True)

  canvas = ROOT.TCanvas("canvas", "canvas", 1500, 1000)
  canvas.Divide(3, 2)

  canvas.cd(1)
  ROOT.gPad.SetLogy()

  correctionFunction = get_timing_correction_function()
  correctionFunction.Draw()

  canvas.cd(2)

  # define a gaussian to describe roughly the eta distribution of photons in the monophoton analysis
  eta_gaussian = ROOT.TF1("eta_gaussian", "gaus", -2.5, 2.5)
  eta_gaussian.SetParameters(1, 0, 1.0)  # amplitude, mean, sigma
  eta_gaussian.SetTitle("Assumed #eta_{#gamma} distribution;#eta_{#gamma};Entries")
  eta_gaussian.Draw()

  canvas.cd(3)

  time_smearing = ROOT.TF1("time_smearing", "gaus", -30, 30)
  time_smearing.SetParameters(1, 0, 1.0)  # amplitude, mean, sigma
  time_smearing.SetTitle("Assumed time resolution;Timing (ns);Entries")
  time_smearing.Draw()

  canvas.cd(4)

  # z_distribution = ROOT.TF1("z_distribution", "1", -30, 30)
  # z_distribution = ROOT.TF1("z_distribution", "abs(x) < 5 ? 1 : 0", -30, 30)
  # z_distribution = ROOT.TF1("z_distribution", "gaus", -30, 30)
  # z_distribution.SetParameters(1, 0, 5)  # amplitude, mean, sigma

  beam_intersection_width_nominal = 0.07  # width of the beam intersection region in meters
  beam_intersection_width_extra = 0.1
  beams_spacing = 2.5  # ns
  beams_distance = beams_spacing * c * 1e-9  # convert to meters
  beam_intersection_point = beams_distance / 2  # position of the beam intersection point in meters

  n_gausses = 5

  function = ""

  for i in range(n_gausses):
    if i > 0:
      function += " + "
    function += f"gaus({i*3})"

  z_distribution = ROOT.TF1("z_distribution", function, -3, 3)

  for i in range(n_gausses):
    # the middle gaussian should be centered at zero, others are spaced equally around it
    if i == n_gausses // 2:
      mean = 0  # center the middle gaussian at zero
      width = beam_intersection_width_nominal
    else:
      mean = (i - n_gausses // 2) * beams_distance  # space gaussians equally around zero
      width = beam_intersection_width_extra

    z_distribution.SetParameter(i * 3, 1)  # amplitude
    z_distribution.SetParameter(i * 3 + 1, mean)  # mean
    z_distribution.SetParameter(i * 3 + 2, width)  # sigma
  # z_distribution.SetParameters(
  #   1, beam_intersection_point, beam_intersection_width_extra, 1, 0, beam_intersection_width_nominal, 1, -beam_intersection_point,
  #   beam_intersection_width_extra
  # )  # amplitude1, mean1, sigma1, amplitude2, mean2, sigma2

  z_distribution.SetNpx(1000)

  z_distribution.SetTitle("Assumed z distribution of photons;Photon z (m);Entries")
  z_distribution.Draw()
  z_distribution.GetXaxis().SetRangeUser(-30, 30)

  canvas.cd(5)

  time_hist_uniform_z = ROOT.TH1F("time_hist", ";Timing (ns);Entries", 50, -30, 30)
  time_hist_zero_z = ROOT.TH1F("time_hist_zero_z", ";Timing (ns);Entries", 50, -30, 30)

  def get_corrected_time(z_gamma, eta_gamma):
    theta = 2 * atan(exp(-eta_gamma))
    z_gamma_prime = z_gamma - ecal_inner_radius / tan(theta)
    d = ecal_inner_radius / sin(theta)
    t_gamma = d / c * 1e9  # convert to ns
    t_smearing = time_smearing.GetRandom()
    t_gamma += t_smearing

    time_correction = correctionFunction.Eval(z_gamma_prime)
    time = t_gamma - time_correction
    return time

  def get_ecal_eta(z_gamma, eta_gamma):
    theta = 2 * atan(exp(-eta_gamma))
    z_gamma_prime = z_gamma - ecal_inner_radius / tan(theta)
    theta_zero = atan(ecal_inner_radius / z_gamma_prime)
    eta_ecal = -sgn(theta_zero) * (log(tan(abs(theta_zero) / 2)))
    return eta_ecal

  eta_zero = ROOT.TH1D("eta_zero", ";#eta_{#gamma};Entries", 200, -10, 10)
  eta_dist = ROOT.TH1D("eta_dist", ";#eta_{#gamma};Entries", 200, -10, 10)

  # simulate photons from z=0
  for _ in range(n_events):
    z_gamma = 0
    eta_gamma = eta_gaussian.GetRandom()
    time = get_corrected_time(z_gamma, eta_gamma)
    time_hist_zero_z.Fill(time)
    eta_zero.Fill(get_ecal_eta(z_gamma, eta_gamma))

  # simulate photons from anywhere along Z
  for _ in range(n_events):
    z_gamma = z_distribution.GetRandom()
    eta_gamma = eta_gaussian.GetRandom()

    z_gamma_prime = z_gamma - ecal_inner_radius / tan(2 * atan(exp(-eta_gamma)))
    z_gamma_prime_back = z_gamma - ecal_inner_radius / tan(2 * atan(exp(eta_gamma)))

    # if both photons miss ECAL:
    if abs(z_gamma_prime) > ecal_half_length and abs(z_gamma_prime_back) > ecal_half_length:
      continue

    # if both photons withing ECAL:
    if abs(z_gamma_prime) < ecal_half_length and abs(z_gamma_prime_back) < ecal_half_length:
      continue

    time = get_corrected_time(z_gamma, eta_gamma)
    time_hist_uniform_z.Fill(time)
    eta_dist.Fill(get_ecal_eta(z_gamma, eta_gamma))

  time_hist_uniform_z.SetLineColor(ROOT.kViolet + 2)
  time_hist_zero_z.SetLineColor(ROOT.kGreen + 2)
  time_hist_zero_z.SetFillColor(ROOT.kGreen + 2)

  ROOT.gPad.SetLogy()
  ROOT.gStyle.SetOptStat(0)

  time_hist_zero_z.Scale(0.1)
  time_hist_zero_z.Draw("hist")
  time_hist_zero_z.GetYaxis().SetRangeUser(1e-2, 1e6)

  time_hist_uniform_z.Draw("samePE")

  legend = ROOT.TLegend(0.3, 0.7, 0.9, 0.9)
  legend.AddEntry(time_hist_zero_z, "z=0", "f")
  legend.AddEntry(time_hist_uniform_z, "z from the assumed distribution", "lep")
  legend.Draw()

  canvas.cd(6)
  ROOT.gStyle.SetOptStat(111111)

  eta_zero.Rebin(5)
  eta_dist.Rebin(5)

  eta_zero.SetLineColor(ROOT.kGreen + 2)
  eta_zero.SetFillColor(ROOT.kGreen + 2)
  eta_zero.DrawNormalized("hist")

  eta_dist.SetLineColor(ROOT.kViolet + 2)
  eta_dist.DrawNormalized("samePE")

  canvas.SaveAs("../plots/timing_correction_function.pdf")


if __name__ == "__main__":
  main()
