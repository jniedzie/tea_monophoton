import glob
import os

# input_path = "/eos/cms/store/cmst3/group/lightbylight/Pranati/Final_afterTrigger/Data/Data_29thJune"
# input_pattern = "ntuples_loose_selections_*.root"

base_path = "/data/dust/user/jniedzie/monophoton/"

# sample = "qed_superchic"
# sample = "collisionData"
sample = "qed_mg1gamma"

dir = "noTrigger"
# dir = "initial_singleEG5"


input_path = f"{base_path}/{sample}/{dir}"
input_pattern = "*.root"

output_path = f"{base_path}/{sample}/{dir}_merged"
output_pattern = "ntuple_{}.root"

n_files_to_merge = 20


def get_file_paths():
  path_pattern = os.path.join(input_path, input_pattern)
  file_paths = glob.glob(path_pattern)
  print(f"Found {len(file_paths)} files matching pattern {path_pattern}")
  return file_paths


def get_file_name(index):
  output_filename = output_pattern.format(index)
  output_filename = os.path.join(output_path, output_filename)
  return output_filename


def merge_batch_of_files(files_to_merge, output_ntuple_counter):
  output_filename = get_file_name(output_ntuple_counter)
  command = f"hadd -f -j -k {output_filename} {' '.join(files_to_merge)}"
  print("Executing commnand:", command)
  os.system(command)


def merge_n_files():
  file_paths = get_file_paths()
  files_to_merge = []
  output_ntuple_counter = 0

  # merge files in batches of n_files_to_merge:
  for filename in file_paths:

    # if full, merge and move to the next batch:
    if len(files_to_merge) == n_files_to_merge:
      merge_batch_of_files(files_to_merge, output_ntuple_counter)
      files_to_merge = []
      output_ntuple_counter += 1

    files_to_merge.append(filename)

  # merge any remaining files:
  if len(files_to_merge) != 0:
    merge_batch_of_files(files_to_merge, output_ntuple_counter)


def main():
  # create output directory if it doesn't exist:
  if not os.path.exists(output_path):
    os.makedirs(output_path)

  merge_n_files()


if __name__ == "__main__":
  main()
