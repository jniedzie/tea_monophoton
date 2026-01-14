# facility = "lxplus"
facility = "NAF"

# trigger = "doubleEG2"
trigger = "singleEG5"

processes = (
  "collisionData",
  "lbl",
  "ds_from_lbl",
  "cep",
  "qed_superchic",
  "qed_starlight",
  "alps_5",
  "alps_30",
  "alps_90",
)

qed_names = ["qed_superchic", "qed_starlight"]

input_skim = f"initial_{trigger}"

# skim = f"skimmed_{trigger}_lblSelections"
# skim = f"skimmed_{trigger}_lblSelections_et-gt4"
# skim = f"skimmed_{trigger}_lblSelections_et-gt4_widthCuts"
# skim = f"skimmed_{trigger}_lblSelections_et-gt4_widthCuts_shapeCuts"
# skim = f"skimmed_{trigger}_lblSelections_et-gt4_widthCuts_shapeCuts_barrelOnly"
# skim = f"skimmed_{trigger}_lblSelections_et-gt4_widthCuts_shapeCuts_deltaT"
skim = f"skimmed_{trigger}_lblSelections_et-gt4_widthCuts_shapeCuts_barrelOnly_deltaT"


if facility == "NAF":
  base_initial_path = "/data/dust/user/jniedzie/light_by_light/ntuples/"
  base_path = "/data/dust/user/jniedzie/monophoton/"
elif facility == "lxplus":
  base_initial_path = "/eos/cms/store/cmst3/group/lightbylight/tea_samples"
  base_path = "/eos/cms/store/cmst3/group/lightbylight/tea_samples"

merged_histograms_path = base_path + "/{}/merged_{}_histograms.root"
