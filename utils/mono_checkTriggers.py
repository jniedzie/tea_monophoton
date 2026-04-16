import os
import ROOT

maxFilesToCheck = 10

# directory = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/emptyBX/initial_badNames"
# directory = "../"
directory = "/eos/cms/store/group/phys_diffraction/lbyl_2018/HIEmptyBX/ntuples_3_11_2026/HIEmptyBX/ntuples_emptyBx/260409_135249/0000"


# trigger_name = "HLT_HIUPC_SingleEG5_NotMBHF2AND_v1"
trigger_name = "HLT_HIUPC_SingleEG5_NotMBHF2AND_v1_Prescl"


def countPassingEvents(hlt_tree, total_events):
  if hlt_tree.GetBranch(trigger_name) is None:
    return total_events

  total_events += hlt_tree.GetEntries(f"{trigger_name} == 1")

  return total_events


def count_trigger_events():
  total_events = 0

  for i_files, filename in enumerate(os.listdir(directory)):
    if i_files >= maxFilesToCheck:
      break

    if not filename.endswith(".root"):
      continue

    file_path = os.path.join(directory, filename)

    print(f"Checking file: {file_path}")

    root_file = ROOT.TFile(file_path)
    if root_file.IsZombie():
      print(f"Error opening file: {file_path}")
      continue
    hlt_tree = root_file.Get("hltanalysis/HltTree")
    if hlt_tree is None:
      print(f"HLT tree not found in file: {file_path}")
      continue

    passing_events_this_file = countPassingEvents(hlt_tree, total_events)

    root_file.Close()
  return total_events


def main():
  ROOT.gROOT.SetBatch(True)

  events_with_trigger = count_trigger_events()
  print(f"Total events with trigger '{trigger_name}': {events_with_trigger}")


if __name__ == "__main__":
  main()
