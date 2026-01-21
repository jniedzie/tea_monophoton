import ROOT
from Sample import Sample, SampleType
from Legend import Legend
from Histogram import Histogram, Histogram2D
from HistogramNormalizer import NormalizationType
from lbl_helpers import get_cep_scale
from lbl_params import luminosity, crossSections, nGenEvents, get_scale_factor, total_uncertainty_qed, total_uncertainty_lbl_run2
from lbl_paths import base_path, processes, skim


output_path = f"../plots/{skim.replace('skimmed_', '')}/"

do_photons = True
do_alps = True

single_photon = True

lbl_error = total_uncertainty_lbl_run2 - 1
qed_error = total_uncertainty_qed - 1

scale = 1.0

print(f"{get_scale_factor(do_photons, single_photon)=}")

samples = [
    Sample(
        name="qed_superchic",
        file_path=f"{base_path}/qed_superchic/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=scale*crossSections["qed_superchic"]*get_scale_factor(do_photons, single_photon)[0],
        initial_weight_sum=nGenEvents["qed_superchic"],
        fill_color=ROOT.kYellow,
        fill_alpha=1.0,        
        marker_size=0.0,
        legend_description="#gamma#gamma#rightarrowe^{+}e^{-}",
        custom_legend=Legend(0.62, 0.70, 0.82, 0.75, "FL"),
    ),
    Sample(
        name="qed_starlight",
        file_path=f"{base_path}/qed_starlight/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=scale*crossSections["qed_starlight"]*get_scale_factor(do_photons, single_photon)[0],
        initial_weight_sum=nGenEvents["qed_starlight"],
        fill_color=ROOT.kYellow,
        line_color=ROOT.kYellow,
        fill_alpha=1.0,
        marker_size=0.0,
        legend_description=""
    ),
    Sample(
        name="data",
        file_path=f"{base_path}/collisionData/merged_{skim}_histograms.root",
        type=SampleType.data,
        line_color=ROOT.kBlack,
        line_style=ROOT.kSolid,
        marker_style=20,
        marker_size=1.0,
        marker_color=ROOT.kBlack,
        fill_alpha=0.0,
        legend_description="Data",
        custom_legend=Legend(0.62, 0.80, 0.82, 0.90, "pl", ""),
    ),
    Sample(
        name="ds_from_lbl",
        file_path=f"{base_path}/ds_from_lbl/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=crossSections["lbl"]*get_scale_factor(do_photons, single_photon)[0]*82,
        initial_weight_sum=nGenEvents["lbl"],
        fill_color=ROOT.kGreen+2,
        fill_alpha=1.0,
        legend_description="DS (from LbL)",
        custom_legend=Legend(0.62, 0.75, 0.82, 0.80, "FL"),
    ),
    Sample(
        name="lbl",
        file_path=f"{base_path}/lbl/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=crossSections["lbl"]*get_scale_factor(do_photons, single_photon)[0],
        initial_weight_sum=nGenEvents["lbl"],
        fill_color=ROOT.kOrange+1,
        fill_alpha=1.0,
        legend_description="#gamma#gamma#rightarrow#gamma#gamma",
        custom_legend=Legend(0.62, 0.65, 0.82, 0.70, "FL"),
    ),
    Sample(
        name="cep",
        file_path=f"{base_path}/cep/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=crossSections["cep"]*get_scale_factor(do_photons, single_photon)[0]*get_cep_scale()[0],
        initial_weight_sum=nGenEvents["cep"],
        fill_color=ROOT.kAzure-4,
        fill_alpha=1.0,
        legend_description="gg#rightarrow#gamma#gamma",
        custom_legend=Legend(0.62, 0.60, 0.82, 0.65, "FL"),
    ),
    Sample(
        name="qed_mg1gamma",
        file_path=f"{base_path}/qed_mg1gamma/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=scale*crossSections["qed_mg1gamma"]*get_scale_factor(do_photons, single_photon)[0],
        initial_weight_sum=nGenEvents["qed_mg1gamma"],
        fill_color=ROOT.kRed-1,
        fill_alpha=1.0,
        marker_size=0.0,
        legend_description="#gamma#gamma#rightarrowe^{+}e^{-}#gamma",
        custom_legend=Legend(0.62, 0.55, 0.82, 0.60, "FL"),
    ),
]

custom_stacks_order = ["qed_mg1gamma", "cep", "lbl", "qed_starlight", "qed_superchic", "ds_from_lbl", "data"]

alp_colors = (
    ROOT.kGray+2,
    ROOT.kCyan+1,
    # ROOT.kCyan,
    ROOT.kMagenta,
    ROOT.kViolet,
    ROOT.kBlue,
    ROOT.kGreen,
    ROOT.kYellow+1,
    ROOT.kOrange,
    ROOT.kRed,
)

legend_width = 0.05
legend_min_x = 0.72
legend_max_x = 0.82

legend_height = 0.05
legend_max_y = 0.85


if do_alps:
    alp_scale = 1.0
    alp_index = 0
    for process in processes:
        if "alps" not in process:
            continue

        legend_x_0 = 0.40 if alp_index < 5 else 0.55
        legend_x_1 = 0.47 if alp_index < 5 else 0.62
        
        legend_y_0 = 0.85 - alp_index*legend_height if alp_index < 5 else 0.85 - (alp_index-5)*legend_height
        legend_y_1 = 0.90 - alp_index*legend_height if alp_index < 5 else 0.90 - (alp_index-5)*legend_height

        samples.append(
            Sample(
                name=process,
                file_path=f"{base_path}/{process}/merged_{skim}_histograms.root",
                type=SampleType.signal,
                cross_section=crossSections[process]*get_scale_factor(do_photons, single_photon)[0]*alp_scale,
                initial_weight_sum=nGenEvents[process],
                line_color=alp_colors[alp_index],
                line_style=ROOT.kSolid,
                line_width=3,
                # fill_style=0,
                fill_color=alp_colors[alp_index],
                fill_alpha=0.2,
                marker_size=0.0,
                legend_description=process.replace("alps_", "m_{a} = ")+" GeV",
                custom_legend=Legend(legend_x_0, legend_y_0, legend_x_1, legend_y_1, "l")
            )
        )
        custom_stacks_order.append(process)
        alp_index += 1

y_label = "Events"

default_lumi = NormalizationType.to_lumi
# default_lumi = NormalizationType.to_data

histograms = (
  #           name              title  logx logy    norm_type               rebin xmin   xmax  ymin    ymax,    xlabel                ylabel            suffix
  Histogram("cutFlow"             , "", False, True, NormalizationType.to_data, 1, 0, 20, 1e-5, 1e13, "Selection", "#sum genWeight"),
  Histogram("event_ZDCenergyPlus" , "", False, True, NormalizationType.to_lumi, 100, 0, 10000, 1e-1, 1e5, "#sum E_{ZDC}^{+} (GeV)", y_label),
  Histogram("event_ZDCenergyMinus" , "", False, True, NormalizationType.to_lumi, 100, 0, 10000, 1e-1, 1e5,"#sum E_{ZDC}^{-} (GeV)", y_label),
)

histograms2D = (
  #           name                      title logs              norm          rebins  x_range  y_range  z_range   labels
  Histogram2D("egamma_et_vs_goodPhoton_et", "", False, False, False, default_lumi, 1,  1,  0,  12,  0,  12,  0, 50,  "e/#gamma E_{T} (GeV)", "reco-#gamma E_{T} (GeV)", "Counts"),
)

# for prefix in ["", "Barrel_", "EndCap_"]:
for prefix in [""]:
  histograms += (
    #           name                  title logx logy    norm_type                    rebin xmin   xmax  ymin    ymax,    xlabel                ylabel            suffix
    # Histogram(f"goodPhoton_{prefix}et", "", False, True, default_lumi, 1,   0, 20, 1e-2, 5e5, "E_{T}^{#gamma} (GeV)", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}et", "", False, True, default_lumi, 5,   0, 100, 1e-2, 5e5, "E_{T}^{#gamma} (GeV)", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}et", "", False, True, default_lumi, 40,   0, 20, 1e-2, 5e5, "E_{T}^{#gamma} (GeV)", y_label, "_normCheck", lbl_error),
    # Histogram(f"goodPhoton_{prefix}logEt", "", False, True, default_lumi, 5,   0.3, 2.6, 1e-2, 5e3, "log_{10}[E_{T}^{#gamma} (GeV)]", y_label, "", lbl_error),

    Histogram(f"goodPhoton_{prefix}eta", "", False, True, default_lumi, 1,   -3, 3, 1e-2, 5e5, "#eta^{#gamma}", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}phi", "", False, True, default_lumi, 2,   -4, 4, 1e-2, 5e5, "#phi^{#gamma}", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}seedTime", "", False, True, default_lumi, 2,   -4, 4, 1e-2, 5e5, "Photon seed time", y_label, "", lbl_error),

    Histogram(f"goodPhoton_{prefix}SCEtaWidth" , "", False, True, default_lumi, 1,   0, 0.01, 1e-2, 5e5, "#eta^{SC} width", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}SCPhiWidth" , "", False, True, default_lumi, 1,   0, 0.01, 1e-2, 5e5, "#phi^{SC} width", y_label, "", lbl_error),
    # Histogram(f"goodPhoton_{prefix}SCPhiWidth" , "", False, True, default_lumi, 20,   0, 0.1, 1e-2, 5e5, "#phi^{SC} width", y_label, "", lbl_error),

    Histogram(f"goodPhoton_{prefix}verticalOverCentral", "", False, True, default_lumi, 15,   0, 0.5, 1e-2, 3e3, "E_{right+left}/(2*E_{max})", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}horizontalOverCentral", "", False, True, default_lumi, 15,   0, 0.5, 1e-2, 3e3, "E_{top+bottom}/(2*E_{max})", y_label, "", lbl_error),
    # Histogram(f"goodPhoton_{prefix}verticalOverCentral", "", False, True, default_lumi, 5,   0, 0.5, 1e-2, 3e3, "E_{right+left}/(2*E_{max})", y_label, "", lbl_error),
    # Histogram(f"goodPhoton_{prefix}horizontalOverCentral", "", False, True, default_lumi, 5,   0, 0.5, 1e-2, 3e3, "E_{top+bottom}/(2*E_{max})", y_label, "", lbl_error),

    Histogram(f"goodPhoton_{prefix}topOverCentral"   , "", False, True, default_lumi, 1,   0, 0.5, 1e-2, 3e3, "E_{top}/E_{max}", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}bottomOverCentral", "", False, True, default_lumi, 1,   0, 0.5, 1e-2, 3e3, "E_{bottom}/E_{max}", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}leftOverCentral"  , "", False, True, default_lumi, 1,   0, 0.5, 1e-2, 3e3, "E_{left}/E_{max}", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}rightOverCentral" , "", False, True, default_lumi, 1,   0, 0.5, 1e-2, 3e3, "E_{right}/E_{max}", y_label, "", lbl_error),

    Histogram(f"goodPhoton_{prefix}verticalImbalance"  , "", False, True, default_lumi, 5,   -2, 2, 1e-3, 3e5, "E_{top-bottom}/E_{top+bottom}", y_label, "", lbl_error),
    Histogram(f"goodPhoton_{prefix}horizontalImbalance", "", False, True, default_lumi, 5,   -2, 2, 1e-3, 3e5, "E_{left-right}/E_{left+right}", y_label, "", lbl_error),

    Histogram(f"goodPhoton_{prefix}sigmaIEtaIEta2012"  , "", False, True, default_lumi, 1,   0, 0.06, 1e-1, 3e5, "#sigma_{i#eta i#eta, 2012}", y_label, "", lbl_error),    
  )

  histograms2D += (
    #           name                      title logs              norm          rebins  x_range  y_range  z_range   labels
    Histogram2D(f"goodPhoton_{prefix}eta_vs_phi", "", False, False, False, default_lumi, 5,  5,  -3,  3,  -4,  4,  0, 1e3,  "#eta", "#phi", "Counts"),  
  )

histogramsRatio = []

n_signal = len([s for s in samples if s.type ==
               SampleType.signal and s.custom_legend is None])
n_data = len([s for s in samples if s.type ==
             SampleType.data and s.custom_legend is None])
n_background = len([s for s in samples if s.type ==
                   SampleType.background and s.custom_legend is None])

# here default legends per sample type are defined. If you want to override them, specify custom_legend in the sample
legends = {
    SampleType.signal: Legend(legend_min_x, legend_max_y - n_signal*legend_height, legend_min_x+legend_width, legend_max_y, "l"),
    SampleType.data: Legend(legend_max_x-legend_width, legend_max_y - 2*legend_height, legend_max_x, legend_max_y-legend_height, "pl", title="#gamma#gamma"),
    SampleType.background: Legend(legend_min_x, legend_max_y - (n_background+1)*legend_height, legend_max_x, legend_max_y-legend_height, "f"),
}

plotting_options = {
    SampleType.background: "hist",
    # SampleType.background: "hist nostack e",
    SampleType.signal: "hist nostack",
    SampleType.data: "nostack pe0",
}

canvas_size = (800, 600)
canvas_size_2Dhists = (800, 600)
show_ratio_plots = True
ratio_limits = (0.0, 5.0)

show_cms_labels = True
extraText = "Preliminary"

beam_label = " PbPb @ 5.02 TeV"
lumi_unit = "nb"
# lumi_label_offset = -0.2
lumi_label_offset = 0.1
lumi_label_value = luminosity


# output_formats = ["pdf", "C"]
output_formats = ["pdf"]
