import ROOT
import glob
import os

from lbl_paths import base_path

# directory = f"{base_path}/collisionData/initial"
# pattern = "ntuple_*.root"

# Empty BX samples:
# directory = "/eos/cms/store/cmst3/group/lightbylight/EmptyBx_HIForest/HIEmptyBX/pbpb_Emptybx_2018/190527_113416/0000/"  # 1'690'049
# directory = "/eos/cms/store/group/phys_diffraction/lbyl_2018/HIEmptyBX/pbpb_Emptybx_2018/190527_113416/0000"  # 1'690'049
# directory = "/eos/cms/store/group/phys_diffraction/lbyl_2018/HIEmptyBX/ntuples_3_11_2026/HIEmptyBX/ntuples_emptyBx/260409_135249/0000"  # 55'974'067
# directory = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/emptyBX/initial_noTrigger"  # data_*: 56'853'028, HiForestAOD_*: 1'685'816
# directory = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/emptyBX/initial_UnpairedBptx"  # 31'264'750

# Collision data samples:
# directory = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/collisionData/bad_names_noTrigger"  # 577'628'726
# directory = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/collisionData/bad_names_singleEG5_old"  # 59'780'048
# directory = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/collisionData/bad_names_singleEG5"  # 59`780`048
directory = "/eos/cms/store/cmst3/group/lightbylight/upc_monophoton/ntuples/collisionData/initial_singleEG5"  #

pattern = "*.root"


def count_tree_entries(directory, pattern):

  tree_name_1 = "ggHiNtuplizer/EventTree"
  tree_name_2 = "Events"

  total_entries = 0
  path_pattern = os.path.join(directory, pattern)

  file_paths = glob.glob(path_pattern)

  print(f"Found {len(file_paths)} files matching pattern {path_pattern}")

  for filename in file_paths:
    print(f"Processing file: {filename}")
    try:
      root_file = ROOT.TFile.Open(filename, "READ")
    except OSError:
      print(f"Error opening file {filename}")
      continue
    if root_file.IsOpen():
      tree = root_file.Get(tree_name_1)
      if tree:
        total_entries += tree.GetEntries()
      else:
        tree = root_file.Get(tree_name_2)
        if tree:
          total_entries += tree.GetEntries()
        else:
          print(f"Failed to open tree {tree_name_1} or {tree_name_2}")
      root_file.Close()
    else:
      print(f"Failed to open {filename}")

  return total_entries


def main():
  total = count_tree_entries(directory, pattern)
  print(f"Total entries in all trees: {total}")


if __name__ == "__main__":
  main()
