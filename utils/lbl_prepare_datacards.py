from lbl_helpers import load_histograms, get_cep_scale, add_uncertainties_on_zero
from lbl_helpers import input_histograms, scale_non_cep_histograms
from lbl_params import total_uncertainty_lbl_run2, alp_mc_uncertainty
from lbl_paths import processes, skim, qed_names
from Logger import info, warn, error, fatal
import ROOT
import os

output_path = f"../combine/histograms_{skim}.root"
signal_name = "ds_from_lbl"


def save_output_histograms():
  output_dir = os.path.dirname(output_path)
  os.system(f"mkdir -p {output_dir}")
  output_file = ROOT.TFile(output_path, "recreate")

  scale_non_cep_histograms()

  for process in processes:
    if process not in input_histograms:
      continue

    if "alps" in process:
      continue

    if process == "cep":
      scale, _ = get_cep_scale()
      print(f"CEP {scale=}")
      input_histograms[process].Scale(scale)

    input_histograms[process] = add_uncertainties_on_zero(
      input_histograms[process]
    )
    name = process if process != "collisionData" else "data_obs"
    input_histograms[process].SetName(name)
    output_file.cd()
    input_histograms[process].Write()

  output_file.Close()


def add_datacard_header(file, rates):
  histograms_path = output_path.replace("../combine/", "")

  file += "imax 1  number of channels\n"
  file += "jmax * number of backgrounds\n"
  file += "kmax * number of nuisance parameters\n"
  file += f"shapes * * {histograms_path} "
  file += " $PROCESS $PROCESS_$SYSTEMATIC\n"
  file += "bin bin1\n"
  file += f"observation {rates['collisionData']}\n"

  file += "bin            "
  for name, rate in rates.items():
    if name == "collisionData":
      continue

    file += "bin1 "

  file += "\n"

  file += "process        "
  for name, rate in rates.items():
    if name == "collisionData":
      continue

    file += f"{name} "

  file += "\n"

  file += "process        "

  counter = 0
  for name, rate in rates.items():
    if name == "collisionData":
      continue

    file += f" {counter}"
    counter += 1

  file += "\n"

  return file


def add_datacard_rates(file, rates):
  file += "rate             "

  for rate in rates:
    if rate == "collisionData":
      continue

    file += f"{rates[rate]} "

  file += "\n"

  return file


def add_datacard_nuisances(file, rates):

  file += "bck_syst     lnN    "

  for name, rate in rates.items():
    if name == "collisionData":
      continue

    if name == signal_name:
      file += "- "
    else:
      file += f"{total_uncertainty_lbl_run2} "

  file += "\n"

  file += "bin1   autoMCStats  10\n"
  return file


def save_datacard():
  rates = {}
  for process in processes:
    if process not in input_histograms:
      warn(f"Process {process} not found in input histograms!")
      continue

    if "alps" in process:
      continue

    info(
      f"Integral for process {process}: {input_histograms[process].Integral()}"
    )

    rates[process] = input_histograms[process].Integral()

  # bring the entry at "signal_name" to the first position in the rates dict
  rates = {signal_name: rates.pop(signal_name), **rates}

  output_file = ""
  output_file = add_datacard_header(output_file, rates)
  output_file = add_datacard_rates(output_file, rates)
  output_file = add_datacard_nuisances(output_file, rates)
  outfile = open(output_path.replace(".root", ".txt"), "w")
  outfile.write(output_file)


def get_data_background_chi2():
  data = input_histograms["collisionData"]
  backgrounds = [input_histograms["lbl"], input_histograms["cep"]]

  for qed_name in qed_names:
    backgrounds.append(input_histograms[qed_name])

  chi2 = 0
  for i in range(1, data.GetNbinsX() + 1):
    data_bin = data.GetBinContent(i)
    background_bin = backgrounds[0].GetBinContent(i) + backgrounds[
      1].GetBinContent(i) + backgrounds[2].GetBinContent(i)

    if data_bin != 0:
      chi2 += (data_bin - background_bin)**2 / data_bin
  return chi2


def main():
  ROOT.gROOT.SetBatch(True)

  info(f"Storing datacard/root file in: {output_path}")

  load_histograms(skim)

  for process in processes:
    if process not in input_histograms:
      continue
    input_histograms[process].SetBinErrorOption(ROOT.TH1.kPoisson)

  for process in processes:
    if process not in input_histograms:
      continue
    input_histograms[process].SetBinErrorOption(ROOT.TH1.kPoisson)

  save_output_histograms()
  save_datacard()


if __name__ == "__main__":
  main()
