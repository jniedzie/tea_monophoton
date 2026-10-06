from lbl_paths import processes, skim, input_base_path, output_base_path

samples = processes
sample_path = ""

# condor:
input_directory = f"{input_base_path}/{sample_path}/{skim}/"
output_hists_dir = f"{output_base_path}/{sample_path}/{skim}/histograms/"


# local:
# input_output_file_list = (
#     (f"{base_path}/ntuples/{sample_path}/merged_{skim}.root",
#      f"{base_path}/ntuples/{sample_path}/merged_{skim}_histograms.root"
#      ),
# )
