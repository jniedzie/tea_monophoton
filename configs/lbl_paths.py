from teaHelpers import get_facility

from lbl_params import eventCuts, zdcCutNames, photonCuts

# trigger = "doubleEG2"
trigger = "singleEG5"
# trigger = "UnpairedBptx"
# trigger = "noTrigger"

do_trigger_selection = True

processes = (
  "collisionData",
  # "ds_from_lbl",
  # "qed_superchic",
  # "qed_starlight",
  # "lbl",
  # "cep",
  # "qed_mg1gamma",
  # "alps_5",
  # "alps_30",
  # "alps_90",
  # "emptyBX",
  # "zeroBias",
)

qed_names = ["qed_superchic", "qed_starlight"]

if do_trigger_selection:
  # input_skim = "initial_noTrigger"
  # skim = f"initial_{trigger}"
  
  # input_skim = "initial_noTrigger_unmerged"
  # skim = f"initial_{trigger}_unmerged"
  
  input_skim = "bad_names_noTrigger"
  skim = f"bad_names_{trigger}"
  
else:
  input_skim = f"initial_{trigger}"

  # figure out cuts combination
  zdcCut = zdcCutNames[eventCuts["ZDC_cut"]]
  timeCut = "_noTimeCut" if (photonCuts["min_seedTime"] < -900 and photonCuts["max_seedTime"] > 900) else ""
  swissCrossCut = "_noSwissCrossCut" if photonCuts["max_swissCross"] >= 1.0 else ""
  scPhiWidthCut = "_tightPhiWidth" if photonCuts["max_SCPhiWidth_barrel"] < 0.15 else ""
  VHfractionsCut = "_noVHfractionsCut" if (photonCuts["min_verticalOverCentral"] == 0.0 and photonCuts["min_horizontalOverCentral"] == 0) else ""
  hOverEcut = "_tightHOverE" if (photonCuts["max_hOverE_barrel"] < 0.005 and photonCuts["max_hOverE_endcap"] < 0.005) else ""

  # suffix = ""
  # suffix = "_goodNEE"
  suffix = "_superCleanNEE_superCleanCHE"

  # suffix += "_noPhotonCuts"

  # build skim name
  skim = f"skimmed_{trigger}_baseSelections{timeCut}{swissCrossCut}{scPhiWidthCut}{VHfractionsCut}{hOverEcut}_zdc{zdcCut}{suffix}"

if get_facility() == "naf":
  base_path = "/data/dust/user/jniedzie/monophoton/"
elif get_facility() == "lxplus":
  base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/"

merged_histograms_path = base_path + "/{}/merged_{}_histograms.root"
