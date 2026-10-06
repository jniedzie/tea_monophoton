from mono_paths import skim, qed_names, base_path, processes
from mono_params import luminosity, crossSections, nGenEvents, get_scale_factor, total_uncertainty_mono_run2, alp_mc_uncertainty
from mono_helpers import get_cep_scale

from Sample import Sample, SampleType
from Histogram import Histogram
from HistogramNormalizer import NormalizationType

qed_superchic, qed_starlight = qed_names

do_alps = True
alp_mass = 5

output_path = f"../datacards/{skim.replace('skimmed_', '')}/"

scale_factor = get_scale_factor(True)[0]
mono_error = total_uncertainty_mono_run2 - 1
add_uncertainties_on_zero = False

samples = [
    Sample(
        name="data",
        file_path=f"{base_path}/collisionData/merged_{skim}_histograms.root",
        type=SampleType.data,
    ),
    Sample(
        name="lbl" if do_alps else "signal_lbl",
        file_path=f"{base_path}/lbl/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=crossSections["lbl"]*scale_factor,
        initial_weight_sum=nGenEvents["lbl"],
    ),
    Sample(
        name=qed_superchic,
        file_path=f"{base_path}/{qed_superchic}/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=crossSections[qed_superchic]*scale_factor,
        initial_weight_sum=nGenEvents[qed_superchic],
    ),
    Sample(
        name=qed_starlight,
        file_path=f"{base_path}/{qed_starlight}/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=crossSections[qed_starlight]*scale_factor,
        initial_weight_sum=nGenEvents[qed_starlight],
    ),
    Sample(
        name="cep",
        file_path=f"{base_path}/cep/merged_{skim}_histograms.root",
        type=SampleType.background,
        cross_section=get_cep_scale()[0],
        initial_weight_sum=luminosity,
    )
]

if do_alps:
    samples.append(
        Sample(
            name=f"signal_alps_{alp_mass}",
            file_path=f"{base_path}/alps_{alp_mass}/merged_{skim}_histograms.root",
            type=SampleType.background,
            cross_section=crossSections[f"alps_{alp_mass}"]*scale_factor,
            initial_weight_sum=nGenEvents[f"alps_{alp_mass}"],
        )
    )

histograms = []

if do_alps:
    histograms.append(Histogram("diphotonSR_mass200"      , "", False, False, NormalizationType.to_lumi, 1, 0, 0, 0, 0, "", "", "", mono_error))
else:
    histograms.append(Histogram("diphoton_acoplanarity40" , "", False, False, NormalizationType.to_lumi, 1, 0, 0, 0, 0, "", "", "", mono_error))

nuisances = {
    "bck_syst": {
        "lbl": total_uncertainty_mono_run2,
        "qed": total_uncertainty_mono_run2,
        "qed_starlight": total_uncertainty_mono_run2,
        "cep": total_uncertainty_mono_run2,
    },
    "alp_mc": {
        "signal_alps_5": alp_mc_uncertainty,
        "signal_alps_6": alp_mc_uncertainty,
        "signal_alps_9": alp_mc_uncertainty,
        "signal_alps_11": alp_mc_uncertainty,
        "signal_alps_14": alp_mc_uncertainty,
        "signal_alps_16": alp_mc_uncertainty,
        "signal_alps_22": alp_mc_uncertainty,
        "signal_alps_30": alp_mc_uncertainty,
        "signal_alps_50": alp_mc_uncertainty,
        "signal_alps_90": alp_mc_uncertainty,
    }
}
