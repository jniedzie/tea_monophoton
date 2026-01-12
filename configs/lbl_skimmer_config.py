from lbl_params import *

nEvents = -1

base_path = "/data/dust/user/jniedzie/monophoton/"

# trigger = "doubleEG2"
trigger = "singleEG5"

# sample = "collisionData"
# sample = "lbl"
# sample = "cep"
# sample = "qed"
# sample = "qed_starlight"

# sample = "alps_5"
# sample = "alps_6"
# sample = "alps_9"
sample = "alps_11"
# sample = "alps_14"
# sample = "alps_16"
# sample = "alps_22"
# sample = "alps_30"
# sample = "alps_50"
# sample = "alps_90"

input_skim = "initial"
output_skim = f"skimmed_baselineSelections_{trigger}"

inputFilePath = f"{base_path}/{sample}/initial_{trigger}/ntuple_0.root"
treeOutputFilePath = inputFilePath.replace("initial", output_skim)

# weightsBranchName = "genWeight"
eventsTreeNames = ["Events",]

# define simple event-level selections
eventSelections = {}

specialBranchSizes = {}
branchesToKeep = ["*"]
branchesToRemove = []
