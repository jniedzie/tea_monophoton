from lbl_paths import trigger, facility

nEvents = -1

if facility == "naf":
  base_path = "/data/dust/user/jniedzie/monophoton/"
elif facility == "lxplus":
  base_path = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/"


input_sample = "lbl"
output_sample = f"ds_from_lbl_{trigger}"


file_name = "mc_HiForestAOD_1.root"


inputFilePath = f"{base_path}/{input_sample}/initial_{trigger}/{file_name}"
treeOutputFilePath = f"{base_path}/{output_sample}/initial_{trigger}/{file_name}"


# weightsBranchName = "genWeight"
eventsTreeNames = ["Events",]

# define simple event-level selections
eventSelections = {}

specialBranchSizes = {}
branchesToKeep = ["*"]
branchesToRemove = []
