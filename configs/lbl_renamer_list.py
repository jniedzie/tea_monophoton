from lbl_paths import processes, base_path
samples = processes
sample_path = ""

input_dir_name = "bad_names_singleEG5"
output_dir_name = "initial_singleEG5"

# input_dir_name = "bad_names_noTrigger"
# output_dir_name = "initial_noTrigger_unmerged"

# input_dir_name = "noTrigger_merged"
# output_dir_name = "initial_noTrigger"

# input_directory = "/eos/cms/store/group/phys_diffraction/lbyl_2018/HIEmptyBX/ntuples_3_11_2026/HIEmptyBX/ntuples_emptyBx/260409_135249/0000"
# input_directory = "/eos/cms/store/group/phys_diffraction/lbyl_2018/HIZeroBias/ntuples_16_04_2026/CRAB_UserFiles/ntuples_zeroBias/260420_084846/0000/"

input_directory = f"{base_path}/{sample_path}/{input_dir_name}/"
output_trees_dir = f"{base_path}/{sample_path}/{output_dir_name}/"
