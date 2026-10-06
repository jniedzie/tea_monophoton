from teaHelpers import get_facility

import os
import sys
from Logger import info

facility = get_facility()

trigger = "singleEG5"
# trigger = "doubleEG2"
# trigger = "UnpairedBptx"
# trigger = "noTrigger"

processes = (
    # "collisionData",
    # "ds_from_lbl",
    # "qed_superchic",
    # "qed_starlight",
    # "lbl",
    # "cep",
    # "gamma_y",
    # "upsilon",
    # "tautau",
    "tautau_allDecays",

    # "qed_mg1gamma",
    # "alps_5",
    # "alps_30",
    # "alps_90",
    # "emptyBX",
    # "zeroBias",
)

qed_names = ["qed_superchic", "qed_starlight"]

# Check if this is to preselect on trigger:
try:
  executable = os.readlink("/proc/self/exe")
except OSError:
  executable = sys.executable
command = [executable]
command.extend(sys.argv)
command = " ".join(command)

if "trigger_selector" in command:
  info("\n\nRunning trigger selection\n\n")
  do_trigger_selection = True
  bad_names_input = True
else:
  info("\n\nRunning non-trigger apps\n\n")
  do_trigger_selection = False
  bad_names_input = None


# Figure out input/output paths
if do_trigger_selection:
  # input_skim = "initial_noTrigger"
  # skim = f"initial_{trigger}"

  # input_skim = "initial_noTrigger_unmerged"
  # skim = f"initial_{trigger}_unmerged"

  input_skim = "bad_names_noTrigger"
  skim = f"bad_names_{trigger}"

else:
  input_skim = f"initial_{trigger}"
  skim = f"skimmed_{trigger}_baseSelections"

  # skim = f"initial_{trigger}"
  # skim = f"skimmed_{trigger}_monoElectronWithIDselections"
  # skim = f"skimmed_{trigger}_gammaPlusElectronWithIDselections"


if facility == "naf":
  input_base_path = "/data/dust/user/jniedzie/monophoton/"
  output_base_path = "/data/dust/user/jniedzie/monophoton/"
elif facility == "lxplus":
  # base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/"

  # input_base_path = "/eos/home-j/jniedzie/mono_cmssw/output_ntuples/"
  # output_base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples"

  # input_base_path = "/eos/cms/store/cmst3/group/lightbylight/new_tea_samples/"
  # output_base_path = "/eos/cms/store/cmst3/group/lightbylight/new_tea_samples/"

  input_base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples"
  output_base_path = input_base_path

  redirector = "eoscms.cern.ch"
elif facility == "vub":
  # input_base_path = "/store/group/phys_diffraction/lbyl_2018/ntuples_2026_muonSegments/HIForward/ntuples_muonSegments/260526_120651/0000"
  # input_base_path = "/store/group/phys_diffraction/lbyl_2018/ntuples_2026_muonSegments/HIForward/ntuples_muonSegments/260526_120651/0001"
  # input_base_path = "/store/group/phys_diffraction/lbyl_2018/ntuples_2026_muonSegments/HIForward/ntuples_muonSegments/260526_120651/0002"
  # input_base_path = "/store/group/phys_diffraction/lbyl_2018/ntuples_2026_muonSegments/HIForward/ntuples_muonSegments/260526_120651/0003"
  # input_base_path = "/store/group/phys_diffraction/lbyl_2018/ntuples_2026_muonSegments/HIForward/ntuples_muonSegments/260526_120651/0004"
  # input_base_path = "/store/group/phys_diffraction/lbyl_2018/ntuples_2026_muonSegments/HIForward/ntuples_muonSegments/260526_120651/0005"
  input_base_path = "/pnfs/iihe/cms/store/user/jniedzie/upc"
  output_base_path = "/pnfs/iihe/cms/store/user/jniedzie/upc"
  base_path = "/pnfs/iihe/cms/store/user/jniedzie/upc"
  redirector = "eoscms.cern.ch"
else:
  # Outside CERN/VUB/NAF, read the shared samples through XRootD.
  # input_base_path = "root://eoscms.cern.ch//eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples"
  # output_base_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ntuples")

  input_base_path = "/Users/jeremi/Documents/Physics/tea_monophoton.nosync/samples"
  output_base_path = "/Users/jeremi/Documents/Physics/tea_monophoton.nosync/"

  redirector = "eoscms.cern.ch"

# Allow local copies of the samples without changing site-specific defaults.
input_base_path = os.environ.get("MONO_INPUT_BASE_PATH", input_base_path)
output_base_path = os.environ.get("MONO_OUTPUT_BASE_PATH", output_base_path)

# Utilities use this alias on every supported facility.
base_path = input_base_path
merged_histograms_path = output_base_path + "/{}/merged_{}_histograms.root"
