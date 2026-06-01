from lbl_paths import processes, skim, base_path
from Logger import info

import os

base_command = "./lbl_histogramer --config lbl_histogramer_config.py --input_path {}/{}/merged_{}.root --output_hists_path {}/{}/merged_{}_histograms.root"

for process in processes:
  command = base_command.format(base_path, process, skim, base_path, process, skim)
  info(f"Running histogramming for {process} with command:\n{command}\n")
  
  os.system(command)
  

# ./lbl_histogramer --config lbl_histogramer_config.py --input_path /pnfs/iihe/cms/store/user/jniedzie/upc/collisionData/merged_skimmed_singleEG5_baseSelections.root --output_hists_path /pnfs/iihe/cms/store/user/jniedzie/upc/collisionData/merged_skimmed_singleEG5_baseSelections_histograms.root
# ./lbl_histogramer --config lbl_histogramer_config.py --input_path /pnfs/iihe/cms/store/user/jniedzie/upc/ds_from_lbl/merged_skimmed_singleEG5_baseSelections.root --output_hists_path /pnfs/iihe/cms/store/user/jniedzie/upc/ds_from_lbl/merged_skimmed_singleEG5_baseSelections_histograms.root
# ./lbl_histogramer --config lbl_histogramer_config.py --input_path /pnfs/iihe/cms/store/user/jniedzie/upc/qed_superchic/merged_skimmed_singleEG5_baseSelections.root --output_hists_path /pnfs/iihe/cms/store/user/jniedzie/upc/qed_superchic/merged_skimmed_singleEG5_baseSelections_histograms.root
# ./lbl_histogramer --config lbl_histogramer_config.py --input_path /pnfs/iihe/cms/store/user/jniedzie/upc/qed_starlight/merged_skimmed_singleEG5_baseSelections.root --output_hists_path /pnfs/iihe/cms/store/user/jniedzie/upc/qed_starlight/merged_skimmed_singleEG5_baseSelections_histograms.root
# ./lbl_histogramer --config lbl_histogramer_config.py --input_path /pnfs/iihe/cms/store/user/jniedzie/upc/lbl/merged_skimmed_singleEG5_baseSelections.root --output_hists_path /pnfs/iihe/cms/store/user/jniedzie/upc/lbl/merged_skimmed_singleEG5_baseSelections_histograms.root
# ./lbl_histogramer --config lbl_histogramer_config.py --input_path /pnfs/iihe/cms/store/user/jniedzie/upc/cep/merged_skimmed_singleEG5_baseSelections.root --output_hists_path /pnfs/iihe/cms/store/user/jniedzie/upc/cep/merged_skimmed_singleEG5_baseSelections_histograms.root
