from lbl_params import *

nEvents = -1

triggerSelection = (
    # "SingleEG5",
    # "SingleEG5_Prescaled",
    # "UnpairedBptxMinus",
    # "UnpairedBptxPlus",
    "HLT_HIUPC_SingleEG5_NotMBHF2AND_v1",
)

eventCuts = {}

# eventsTreeNames = ["Events",]
eventsTreeNames = ["hltanalysis/HltTree",]

eventSelections = {}

specialBranchSizes = {}
branchesToKeep = ["*"]
branchesToRemove = []
