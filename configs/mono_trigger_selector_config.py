from lbl_params import *
from lbl_paths import redirector, trigger, bad_names_input
from Logger import fatal

nEvents = -1


if bad_names_input and trigger == "singleEG5":
  triggerSelection = ("HLT_HIUPC_SingleEG5_NotMBHF2AND_v1",)
elif bad_names_input and trigger == "doubleEG2":
  triggerSelection = ("HLT_HIUPC_DoubleEG2_NotMBHF2AND_v1",)
elif not bad_names_input and trigger == "singleEG5":
  triggerSelection = ("SingleEG5",)
elif not bad_names_input and trigger == "doubleEG2":
  triggerSelection = ("DoubleEG2",)
elif not bad_names_input and trigger == "UnpairedBptx":
  triggerSelection = ("UnpairedBptxMinus", "UnpairedBptxPlus",)
else:
  fatal(f"Wrong trigger selection - check mono_trigger_selector_config.py")

eventCuts = {}

if bad_names_input:
  eventsTreeNames = ["hltanalysis/HltTree",]
else:
  eventsTreeNames = ["Events",]

eventSelections = {}

specialBranchSizes = {}
branchesToKeep = ["*"]
branchesToRemove = []

customRedirector = redirector
