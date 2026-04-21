import ROOT
from math import atan, cos, exp, log, pi, sin, sqrt, tan


ECAL_INNER_RADIUS = 1.29  # ECAL inner radius in meters
ECAL_HALF_LENGTH = 3.0  # half-length of the ECAL in meters
SPEED_OF_LIGHT = 3e8  # speed of light in m/s

N_EVENTS = 10000
HIST_BINS = 50
HIST_RANGE = (-1000, 1000)
Z_HIST_RANGE = (-100, 100)
ETA_RANGE = (-2.5, 2.5)
ETA_HIST_RANGE = (-10, 10)
ETA_HIST_BINS = 200
ETA_REBIN_FACTOR = 5

CANVAS_NAME = "canvas"
CANVAS_SIZE = (1500, 1000)
OUTPUT_PATH = "../plots/timing_pion_toy_model.pdf"

PI0_MASS_GEV = 0.1349768
PION_GAMMA_RANGE = (1.05, 50.0)
PION_ETA_RANGE = (-30, 30)
PHOTON_ENDCAP_MAX_ETA = 2.2


def sgn(x):
  return 1 if (x > 0) else -1


def get_theta(eta):
  return 2 * atan(exp(-eta))


def get_ecal_intersection_z(z_gamma, eta_gamma):
  theta = get_theta(eta_gamma)
  return z_gamma + ECAL_INNER_RADIUS / tan(theta)


def get_endcap_outer_radius():
  theta_limit = get_theta(PHOTON_ENDCAP_MAX_ETA)
  return ECAL_HALF_LENGTH * tan(theta_limit)


def get_endcap_intersection_radius(z_gamma, eta_gamma):
  theta = get_theta(eta_gamma)
  z_endcap = ECAL_HALF_LENGTH if eta_gamma >= 0 else -ECAL_HALF_LENGTH
  return abs(z_endcap - z_gamma) * tan(theta)


def photon_reaches_ecal(z_gamma, eta_gamma):
  return abs(get_ecal_intersection_z(z_gamma, eta_gamma)) < ECAL_HALF_LENGTH


def photon_reaches_allowed_endcap(z_gamma, eta_gamma):
  if eta_gamma == 0:
    return False

  z_endcap = ECAL_HALF_LENGTH if eta_gamma >= 0 else -ECAL_HALF_LENGTH
  delta_z = z_endcap - z_gamma

  if delta_z == 0 or delta_z * eta_gamma <= 0:
    return False

  return get_endcap_intersection_radius(z_gamma, eta_gamma) <= get_endcap_outer_radius()


def is_single_ecal_photon(z_gamma, eta_gamma):
  reaches_forward_ecal = photon_reaches_ecal(z_gamma, eta_gamma)
  reaches_backward_ecal = photon_reaches_ecal(z_gamma, -eta_gamma)
  return reaches_forward_ecal != reaches_backward_ecal


def get_random_uniform(min_value, max_value):
  return ROOT.gRandom.Uniform(min_value, max_value)


def get_corrected_time(z_gamma, eta_gamma, time_smearing, correction_function):
  theta = get_theta(eta_gamma)
  distance_to_ecal = ECAL_INNER_RADIUS / sin(theta)
  smeared_arrival_time = distance_to_ecal / SPEED_OF_LIGHT * 1e9 + time_smearing.GetRandom()
  time_correction = correction_function.Eval(get_ecal_intersection_z(z_gamma, eta_gamma))
  return smeared_arrival_time - time_correction


def get_ecal_eta(z_gamma, eta_gamma):
  z_ecal = get_ecal_intersection_z(z_gamma, eta_gamma)
  theta_ecal = atan(ECAL_INNER_RADIUS / z_ecal)
  return -sgn(theta_ecal) * log(tan(abs(theta_ecal) / 2))


def get_momentum_from_gamma(mass, gamma_factor):
  return mass * sqrt(gamma_factor * gamma_factor - 1.0)


def create_pion_eta_distribution():
  pion_eta_distribution = ROOT.TF1("pion_eta_distribution", "gaus(0) + gaus(3)", *PION_ETA_RANGE)
  pion_eta_distribution.SetParameters(1, -4.3, 0.1, 1, 4.3, 0.1)
  pion_eta_distribution.SetTitle("Assumed #eta_{#pi^{0}} distribution;#eta_{#pi^{0}};Entries")
  return pion_eta_distribution


def create_pion_gamma_distribution():
  pion_gamma_distribution = ROOT.TF1("pion_gamma_distribution", "gaus", *PION_GAMMA_RANGE)
  pion_gamma_distribution.SetParameters(1, 35, 1.0)
  pion_gamma_distribution.SetTitle("Placeholder #pi^{0} boost distribution;#gamma_{#pi^{0}};Entries")
  return pion_gamma_distribution


def generate_pion_four_vector(pion_eta_distribution, pion_gamma_distribution):
  pion_eta = pion_eta_distribution.GetRandom()
  pion_gamma = pion_gamma_distribution.GetRandom()
  pion_phi = get_random_uniform(-pi, pi)

  pion_theta = get_theta(pion_eta)
  pion_momentum = get_momentum_from_gamma(PI0_MASS_GEV, pion_gamma)
  pion_energy = pion_gamma * PI0_MASS_GEV

  pion_px = pion_momentum * sin(pion_theta) * cos(pion_phi)
  pion_py = pion_momentum * sin(pion_theta) * sin(pion_phi)
  pion_pz = pion_momentum * cos(pion_theta)

  pion = ROOT.TLorentzVector()
  pion.SetPxPyPzE(pion_px, pion_py, pion_pz, pion_energy)
  return pion


def decay_pi0_to_two_photons(pion):
  # Use an isotropic two-body decay in the pion rest frame as a simple baseline model.
  photon_energy_rest = PI0_MASS_GEV / 2.0
  photon_momentum_rest = photon_energy_rest

  decay_cos_theta = get_random_uniform(-1.0, 1.0)
  decay_sin_theta = sqrt(1.0 - decay_cos_theta * decay_cos_theta)
  decay_phi = get_random_uniform(-pi, pi)

  photon_1 = ROOT.TLorentzVector()
  photon_1.SetPxPyPzE(
    photon_momentum_rest * decay_sin_theta * cos(decay_phi),
    photon_momentum_rest * decay_sin_theta * sin(decay_phi),
    photon_momentum_rest * decay_cos_theta,
    photon_energy_rest,
  )

  photon_2 = ROOT.TLorentzVector()
  photon_2.SetPxPyPzE(-photon_1.Px(), -photon_1.Py(), -photon_1.Pz(), photon_energy_rest)

  boost_vector = pion.BoostVector()
  photon_1.Boost(boost_vector)
  photon_2.Boost(boost_vector)

  return photon_1, photon_2


def get_selected_photon_eta(z_gamma, photons):
  barrel_photons = [photon for photon in photons if photon_reaches_ecal(z_gamma, photon.Eta())]
  if len(barrel_photons) != 1:
    return None

  other_photon = photons[0] if barrel_photons[0] is photons[1] else photons[1]
  other_eta = other_photon.Eta()

  if photon_reaches_ecal(z_gamma, other_eta):
    return None

  if photon_reaches_allowed_endcap(z_gamma, other_eta):
    return None

  return barrel_photons[0].Eta()


def create_canvas():
  canvas = ROOT.TCanvas(CANVAS_NAME, CANVAS_NAME, *CANVAS_SIZE)
  canvas.Divide(3, 2)
  return canvas


def create_timing_correction_function():
  correction_function = ROOT.TF1(
    "correctionFunction",
    "sqrt([0]*[0]+x*x)/[1] * 1e9",
    -15,
    15,
  )
  correction_function.SetParameter(0, ECAL_INNER_RADIUS)
  correction_function.SetParameter(1, SPEED_OF_LIGHT)
  correction_function.SetLineColor(ROOT.kBlue)
  correction_function.SetTitle("Timing Correction Function;z_{ECAL} (m);Timing Correction (ns)")
  return correction_function


def create_eta_distribution():
  eta_distribution = ROOT.TF1("eta_gaussian", "gaus", *ETA_RANGE)
  eta_distribution.SetParameters(1, 0, 1.0)
  eta_distribution.SetTitle("Assumed #eta_{#gamma} distribution;#eta_{#gamma};Entries")
  return eta_distribution


def create_time_smearing():
  time_smearing = ROOT.TF1("time_smearing", "gaus", *HIST_RANGE)
  time_smearing.SetParameters(1, 0, 1.0)
  time_smearing.SetTitle("Assumed time resolution;Timing (ns);Entries")
  return time_smearing


def create_z_distribution():
  z_distribution = ROOT.TF1("z_distribution", "1", *Z_HIST_RANGE)

  # z_distribution = ROOT.TF1("z_distribution", "abs(x) < 5 ? 1 : 0", -30, 30)
  # z_distribution = ROOT.TF1("z_distribution", "gaus", *Z_HIST_RANGE)
  # z_distribution.SetParameters(1, 0, 2)  # amplitude, mean, sigma

  # beam_intersection_width_nominal = 0.1  # width of the beam intersection region in meters
  # beam_intersection_width_extra = 1.0
  # beams_spacing = 100  # ns
  # beams_distance = beams_spacing * SPEED_OF_LIGHT * 1e-9  # convert to meters
  # beam_intersection_point = beams_distance / 2  # position of the beam intersection point in meters

  # z_distribution = ROOT.TF1("z_distribution", "gaus(0) + gaus(3) + gaus(6)", -3, 3)
  # z_distribution.SetParameters(
  #   1, beam_intersection_point, beam_intersection_width_extra, 1, 0, beam_intersection_width_nominal, 1, -beam_intersection_point,
  #   beam_intersection_width_extra
  # )

  z_distribution.SetNpx(1000)
  z_distribution.SetTitle("Assumed z distribution of #pi^{0} production;z_{#pi^{0}} (m);Entries")
  return z_distribution


def create_time_histogram(name):
  return ROOT.TH1F(name, ";Timing (ns);Entries", HIST_BINS, *HIST_RANGE)


def create_eta_histogram(name):
  return ROOT.TH1D(name, ";#eta_{#gamma};Entries", ETA_HIST_BINS, *ETA_HIST_RANGE)


def simulate_zero_z_photons(time_histogram, eta_histogram, eta_distribution, time_smearing, correction_function):
  for _ in range(N_EVENTS):
    z_gamma = 0
    eta_gamma = eta_distribution.GetRandom()
    time_histogram.Fill(get_corrected_time(z_gamma, eta_gamma, time_smearing, correction_function))
    eta_histogram.Fill(get_ecal_eta(z_gamma, eta_gamma))


def simulate_displaced_photons(
  time_histogram,
  eta_histogram,
  z_distribution,
  pion_eta_distribution,
  pion_gamma_distribution,
  time_smearing,
  correction_function,
):
  i_generated = 0
  
  while True:
    if i_generated >= N_EVENTS:
      break
    z_gamma = z_distribution.GetRandom()
    pion = generate_pion_four_vector(pion_eta_distribution, pion_gamma_distribution)
    photons = decay_pi0_to_two_photons(pion)
    eta_gamma = get_selected_photon_eta(z_gamma, photons)

    if eta_gamma is None:
      continue

    time_histogram.Fill(get_corrected_time(z_gamma, eta_gamma, time_smearing, correction_function))
    eta_histogram.Fill(get_ecal_eta(z_gamma, eta_gamma))
    
    i_generated += 1


def draw_input_distributions(
  canvas,
  correction_function,
  eta_distribution,
  time_smearing,
  z_distribution,
):
  canvas.cd(1)
  ROOT.gPad.SetLogy()
  correction_function.Draw()

  canvas.cd(2)
  eta_distribution.Draw()

  canvas.cd(3)
  time_smearing.Draw()

  canvas.cd(4)
  z_distribution.Draw()
  z_distribution.GetXaxis().SetRangeUser(*HIST_RANGE)


def style_time_histograms(time_hist_uniform_z, time_hist_zero_z):
  time_hist_uniform_z.SetLineColor(ROOT.kViolet + 2)
  time_hist_zero_z.SetLineColor(ROOT.kGreen + 2)
  time_hist_zero_z.SetFillColor(ROOT.kGreen + 2)


def draw_time_comparison(canvas, time_hist_uniform_z, time_hist_zero_z):
  canvas.cd(5)
  ROOT.gPad.SetLogy()
  ROOT.gStyle.SetOptStat(0)

  time_hist_zero_z.Scale(0.01)
  time_hist_zero_z.Draw("hist")
  time_hist_zero_z.GetYaxis().SetRangeUser(1e-2, 1e6)

  time_hist_uniform_z.Draw("samePE")

  legend = ROOT.TLegend(0.3, 0.7, 0.9, 0.9)
  legend.AddEntry(time_hist_zero_z, "z=0", "f")
  legend.AddEntry(time_hist_uniform_z, "#pi^{0} decays from the assumed z distribution", "lep")
  legend.Draw()


def draw_eta_comparison(canvas, eta_zero, eta_dist):
  canvas.cd(6)
  ROOT.gStyle.SetOptStat(111111)

  eta_zero.Rebin(ETA_REBIN_FACTOR)
  eta_dist.Rebin(ETA_REBIN_FACTOR)

  eta_zero.SetLineColor(ROOT.kGreen + 2)
  eta_zero.SetFillColor(ROOT.kGreen + 2)
  eta_zero.Draw("hist")
  eta_zero.GetYaxis().SetRangeUser(0, 2000)

  eta_dist.SetLineColor(ROOT.kViolet + 2)
  eta_dist.Draw("samePE")


def main():
  ROOT.gROOT.SetBatch(True)

  canvas = create_canvas()
  correction_function = create_timing_correction_function()
  eta_distribution = create_eta_distribution()
  time_smearing = create_time_smearing()
  z_distribution = create_z_distribution()
  pion_eta_distribution = create_pion_eta_distribution()
  pion_gamma_distribution = create_pion_gamma_distribution()

  time_hist_uniform_z = create_time_histogram("time_hist")
  time_hist_zero_z = create_time_histogram("time_hist_zero_z")
  eta_zero = create_eta_histogram("eta_zero")
  eta_dist = create_eta_histogram("eta_dist")

  draw_input_distributions(
    canvas,
    correction_function,
    eta_distribution,
    time_smearing,
    z_distribution,
  )

  simulate_zero_z_photons(time_hist_zero_z, eta_zero, eta_distribution, time_smearing, correction_function)
  simulate_displaced_photons(
    time_hist_uniform_z,
    eta_dist,
    z_distribution,
    pion_eta_distribution,
    pion_gamma_distribution,
    time_smearing,
    correction_function,
  )

  style_time_histograms(time_hist_uniform_z, time_hist_zero_z)
  draw_time_comparison(canvas, time_hist_uniform_z, time_hist_zero_z)
  draw_eta_comparison(canvas, eta_zero, eta_dist)

  canvas.SaveAs(OUTPUT_PATH)


if __name__ == "__main__":
  main()
