import socket
from Logger import fatal
import sys
from lbl_params import eventCuts

hostname = socket.gethostname()
if "lxplus" in hostname:
  facility = "lxplus"
elif "naf" in hostname:
  facility = "NAF"
else:
  fatal(f"Unknown facility for hostname: {hostname}")
  sys.exit(1)
  exit(1)


# trigger = "doubleEG2"
trigger = "singleEG5"

processes = (
  # "collisionData",
  "ds_from_lbl",
  "qed_superchic",
  "qed_starlight",
  "lbl",
  "cep",
  "qed_mg1gamma",
  # "alps_5",
  # "alps_30",
  # "alps_90",
  # "emptyBX",
)

qed_names = ["qed_superchic", "qed_starlight"]

# trigger selection
# input_skim = "initial_noTrigger"
skim = f"initial_{trigger}"

# skimming
input_skim = f"initial_{trigger}"

zdcCutNames = {
  0: "None",
  1: "LbLstyle",
  2: "0n0n",
  3: "Le1n1n",
  4: "TargetNucleusBreaking",
}

zdcCut = zdcCutNames[eventCuts["ZDC_cut"]]

# skim = f"skimmed_{trigger}_baseSelections_zdc{zdcCut}"

if facility == "NAF":
  base_path = "/data/dust/user/jniedzie/monophoton/"
elif facility == "lxplus":
  base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/"

merged_histograms_path = base_path + "/{}/merged_{}_histograms.root"
