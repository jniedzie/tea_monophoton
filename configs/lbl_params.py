# event cuts
eventCuts = {
  # "max_ZDCenergyPerSide": 10000.0,  # <4n
  "max_ZDCenergyPerSide": 7000.0,  # <3n
  # "max_ZDCenergyPerSide": 4000.0,  # <2n
  # "max_ZDCenergyPerSide": 1600.0,  # <1n
  "max_Nelectrons": 0,
  "max_Ntracks": 0,
  "max_Nmuons": 0,
  "max_Ntowers": 0,
}

# good object definitions
photonCuts = {
  # Common cuts:
  "max_swissCross": 0.95,
  "max_hOverE_barrel": 0.04596,
  "max_hOverE_endcap": 0.0590,
  "max_SCEtaWidth_barrel": 0.0106,
  "max_SCEtaWidth_endcap": 0.0272,
  "max_SCPhiWidth": 999999,  # try 0.01
  "min_sigmaIEtaIEta_barrel": 0,  # try 0.009
  "max_sigmaIEtaIEta_barrel": 0.02,
  "min_sigmaIEtaIEta_endcap": 0,  # try 0.009
  "max_sigmaIEtaIEta_endcap": 0.06,

  # LbL cuts:
  # "min_et": 2.0,
  # "max_absEta": 2.2,
  # "min_SCEtaWidth": 0.0,
  # "min_SCPhiWidth": 0.0,
  # "min_verticalOverCentral": 0.0,
  # "min_horizontalOverCentral": 0.0,
  # "max_seedTime": 3.0,

  # Tightened cuts:
  "min_et": 4.0,
  "min_SCEtaWidth": 0.001,
  "min_SCPhiWidth": 0.001,
  "min_verticalOverCentral": 0.03,
  "min_horizontalOverCentral": 0.03,
  "max_absEta": 1.2,
  "max_seedTime": 1.0,
}

photonHotSpots = {
  # eta_min, eta_max, phi_min, phi_max
  # "hotspot_1": (-1.87, -1.83, 2.20, 2.27),
  # "hotspot_2": (-1.62, -1.58, -2.78, -2.70),
  # "hotspot_3": (-1.60, -1.58, 2.17, 2.21),
  # "hotspot_4": (2.12, 2.14, 0.09, 0.13),
}

dataBlinding = {
  # "max_et": 10.0,  # blind data with photon ET > X GeV
  "max_et": 999999,  # blind data with photon ET > X GeV
}

electronCuts = {
  "min_pt": 2.0,
  "max_absEtaSC": 2.2,
  "max_nMissingHits": 1,
  "max_hOverE": 0.005,
  "max_deltaEtaAtVertex": 0.1,

  # we don't apply electron isolation:
  "max_PFChIso_barrel": 999999,
  "max_PFPhoIso_barrel": 999999,
  "max_PFNeuIso_barrel": 999999,
  "max_PFChIso_endcap": 999999,
  "max_PFPhoIso_endcap": 999999,
  "max_PFNeuIso_endcap": 999999,
}

trackCuts = {
  "min_pt": 0.3,
  "max_absEta": 2.4,
  "min_nValidHits": 4,

  # we don't apply these track selections:
  "max_normalizedChi2": 999999,
  "max_dxy": 999999,
  "max_dz": 999999,
  "max_dxyOverError": 999999,
  "max_dzOverError": 999999,
}

muonCuts = {
  "min_pt": 2.5,
  "max_absEta": 2.4,
}

# calorimeter cuts
caloNoiseThresholds = {
  "HFp": 6.0,
  "HFm": 6.0,

  # "HFp": 9.1, # + 0.8 - 1.3
  # "HFm": 8.8, # + 0.7 - 1.0
  "HB": 2.8,
  "HE": 1.0,
  "EB": 0.7,
  "EE": 3.0,
}

caloNoiseVariables = {
  "HFm": "energy",
  "HFp": "energy",
  "HB": "hadE",
  "HE": "hadE",
  "EB": "emE",
  "EE": "emE",
}

# detector parameters
deadEtas = {
  "HFp": (29, 30),  # 2.853 -- 3.139
  "HFm": (-29, -30),

  # "HFp": (29,), # 2.853 -- 3.000
  # "HFm": (-29,),
  "HE": (-16, 16),  # 1.305 -- 1.392
}

caloEtaEdges = {
  "maxEB": 1.479,
  "minEE": 1.479,
  "maxEE": 3.0,
  "maxHB": 1.305,
  "minHE": 1.305,
  "maxHE": 3.0,
  "minHF": 2.9,
  "maxHF": 5.2,
}

detectorParams = {
  "crack_start": 1.4442,
  "crack_end": 1.566,
  # "crack_end": 1.65,  # to kill weird monophotons
  "crackHadron_start": 1.305,
  "crackHadron_end": 1.41,
  "hem_etaStart": -3.0,
  "hem_etaEnd": -1.39,
  "hem_phiStart": -1.6,
  "hem_phiEnd": -0.9,
  "caloTower_etaMax": 2.4,
}

# matching between calorimeters and photons/electrons
caloMatching = {
  "maxDeltaEta_EB": 0.15,
  "maxDeltaPhi_EB": 0.15,
  "maxDeltaEta_EE": 0.15,
  "maxDeltaPhi_EE": 0.15,
}

# matching between photons and electrons
electronPhotonMatching = {
  "maxDeltaEta": 0.5,  # ???
  "maxDeltaPhi": 0.5,  # ???
}

# matching between tracks and electrons
electronTrackMatching = {
  "maxDeltaEta": 0.15,
  "maxDeltaPhi": 0.7,
}

#  scaling parameters

luminosity = 1647.180726  # μb^-1, with ZDC
# luminosity = 1647.2  # μb^-1, without ZDC
luminosity_err = luminosity * 0.015  # 1.5% uncertainty

reference_alp_cross_section = 10e-3  # μb

# QED needs to be scaled up, depending on the ZDC cut
# qed_scaling = 1.0  # 4n OR
# qed_scaling = 2.2  # 1n AND
# qed_scaling = 2.3  # 2n AND
# qed_scaling = 2.7  # 3n AND
# qed_scaling = 2.5  # 4n AND
# qed_scaling = 1.0  # Xn0n, no scaling
qed_scaling = 0.5  # adding SC + SL
# qed_scaling = 2.047  # no CEP

# LbL may be scaled up due to NLO corrections
lbl_scaling = 1.05  # inclusive
# lbl_scaling = 1.05 * 0.74  # 0n0n
# lbl_scaling = 1.05 * (0.74 + 0.046 + 0.006)  # 0n0n + 0n1n + 1n0n + 1n1n

mc_scale = 1.0
# mc_scale = 10

crossSections = {
  "lbl": mc_scale * 2.59 * lbl_scaling,  # μb
  "qed_superchic": mc_scale * 8827.220 * qed_scaling,  # μb
  "qed_starlight": mc_scale * 7920.0 * qed_scaling,  # μb
  "cep": mc_scale * 5.8e-3,  # we scale it to data
  "alps_5": reference_alp_cross_section,
  "alps_30": reference_alp_cross_section,
  "alps_90": reference_alp_cross_section,

  # "alps_5": mc_scale * 2e2 * 1e-3,  # nb -> μb, limit cross section
  # "alps_30": mc_scale * 5 * 1e-3,  # nb -> μb, limit cross section
  # "alps_90": mc_scale * 5 * 1e-3,  # nb -> μb, limit cross section

  # nb -> μb, g = 0.2 TeV-1
  # "alps_14": 70.21369385210362 * 1e-3,
  # "alps_30": 21.396925059153958 * 1e-3,

  # nb -> μb, g = 0.25 TeV-1
  # "alps_14": 109.82586410834259 * 1e-3,
  # "alps_30": 33.46834007669329 * 1e-3,

  # nb -> μb, g = 0.3 TeV-1
  # "alps_14": 158.28699751455125 * 1e-3,
  # "alps_30": 48.23638862799777 * 1e-3,
}

# photon ET > 2.0 GeV, diphoton pt < 1 GeV
scale_factors = {
  "photonReco": 0.9758,
  "photonID": 0.946,
  "electronRecoID": 0.943,
  "l1eg": 1.0089,
  "l1hf": 0.8716,
  "che": 0.9252,
  "nee": 0.8487,  # old SF
  # "nee": 0.953,  # new SF (higher HF thresholds)
}

scale_factor_errors = {
  "photonReco": 0.0314,
  "photonID": 0.049,
  "electronRecoID": 0.0085,
  "l1eg": 0.002,
  "l1hf": 0.054,
  "che": 0.0087,
  "nee": 0.0085,
}


def get_scale_factor(photon=True, single_photon=False):
  value = 1
  error = 0

  to_skip = "electron" if photon else "photon"
  squared = "photon" if photon else "electron"

  if single_photon:
    squared = "noSquaring"

  for variable in scale_factors:
    if to_skip in variable:
      continue

    error += (scale_factor_errors[variable] / scale_factors[variable])**2
    value *= scale_factors[variable]

    if squared in variable:
      error += (scale_factor_errors[variable] / scale_factors[variable])**2
      value *= scale_factors[variable]

  sf_error = value * error**(1 / 2)

  return value, sf_error


nGenEvents = {
  "lbl": 466000,
  "cep": 668000,  # we scale it to data
  "qed_superchic": 59260000,
  "qed_starlight": 66750000,
  "alps_5": 754000,
  "alps_30": 719000,
  "alps_90": 449000,
}

uncertainty_on_zero = 1.84  # 95% CL
# uncertainty_on_zero = 1.14  # 68% CL

total_uncertainty_qed = 1.068
total_uncertainty_lbl_run2 = 1.23
non_stat_uncertainty_lbl_run2 = 1.18
stat_uncertainty_lbl_run2 = 1.15

total_uncertainty_lbl_run1 = 1.24
alp_mc_uncertainty = 1.03

total_diphoton_efficiency = 0.1352
total_diphoton_efficiency_err = 0.0030
