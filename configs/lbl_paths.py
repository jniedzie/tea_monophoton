from teaHelpers import get_facility

from lbl_params import eventCuts, zdcCutNames, photonCuts

# trigger = "doubleEG2"
trigger = "singleEG5"

processes = (
  "collisionData",
  "ds_from_lbl",
  "qed_superchic",
  "qed_starlight",
  "lbl",
  "cep",
  "qed_mg1gamma",
  "alps_5",
  "alps_30",
  "alps_90",
  # "emptyBX",
)

qed_names = ["qed_superchic", "qed_starlight"]

# trigger selection
# input_skim = "initial_noTrigger"
# skim = f"initial_{trigger}"

# skimming
input_skim = f"initial_{trigger}"

zdcCut = zdcCutNames[eventCuts["ZDC_cut"]]

noTimeCut = photonCuts["min_seedTime"] < -900 and photonCuts["max_seedTime"] > 900
timeCut = "_noTimeCut" if noTimeCut else ""

skim = f"skimmed_{trigger}_baseSelections{timeCut}_zdc{zdcCut}"

if get_facility() == "NAF":
  base_path = "/data/dust/user/jniedzie/monophoton/"
elif get_facility() == "lxplus":
  base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/"

merged_histograms_path = base_path + "/{}/merged_{}_histograms.root"
