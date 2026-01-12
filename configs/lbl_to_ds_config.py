from lbl_paths import trigger

nEvents = -1

base_path = "/data/dust/user/jniedzie/monophoton/"

input_sample = "lbl"
output_sample = "ds_from_lbl"

file_name = "ntuples_loose_selections_3.root"

inputFilePath = f"{base_path}/{input_sample}/initial_{trigger}/{file_name}"
treeOutputFilePath = f"{base_path}/{output_sample}/initial_{trigger}/{file_name}"


# weightsBranchName = "genWeight"
eventsTreeNames = ["Events",]

# define simple event-level selections
eventSelections = {}

specialBranchSizes = {}
branchesToKeep = ["*"]
branchesToRemove = []
