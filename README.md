# Welcome to CMS UPC-monophoton analysis (Run 2)

## Prerequisites

The analysis framework is based on `tea` - please have a look at the general setup instructions here: [cern.ch/tea](https://jniedzie.github.io/tea/docs/build/).

**You don't need to install `tea` though, it will come automatically when you install `tea_monophoton`!!**

If you're going to run on lxplus, just remember to use lxplus9, install the `correctionlib`, and use `python3`.

## Installation

To clone the repo, together with all submodules, run:

```bash
git clone --recurse-submodules git@github.com:jniedzie/tea_monophoton.git
```

Then, simply run `. tea/build.sh` from the `tea_monophoton` directory.

## Example samples

To test things on lxplus, we share LbL MC ntuples that can be found here:
`/eos/cms/store/cmst3/group/lightbylight/tea_samples/lbl/bad_names/`

## How to run

The analysis involves a few steps:
1. Preparing ntuples with names that `tea` can understand out of the HIForest ntuples.
2. Skimming the ntuples (applying selections).
3. Creating histograms.
4. Plotting.
5. Statistical analysis.

All apps and configs are in the `bin` directory - that's where you should run all commands.

### Renaming

The first step is to format ntuples in a way that `tea` can easily understand, also getting rid of some ambiguity in the original HIForest branch names.
You can see the input and output branch names, together with their types, in `lbl_renamer_config.py` - typically there's no need to modify this file.

The input/output paths are defined in `lbl_renamer_list.py` - you may need to adjust them, especially the output path.

Finally, have a look at `lbl_paths.py` - this is where we decide, among other things, which samples to run on. For instance, to use the example samples,
you can comment out everything except for "lbl" in the `processes` tuple.

When all paths look good, you can run the renaming. It's always a good idea to first try locally (and for this tiny sample it will be fast enough to just run it locally anyway):

```
python3 submitter.py --app lbl_renamer --config lbl_renamer_config.py --files_config lbl_renamer_list.py --local
```

You can now check that new files have been created in the output directory, which contain the `Events` tree with updated branch names.

As an excercise, you can also replace the `--local` flag with `--condor` to see how to run things on the grid. It will automatically schedule 1 job per file and produce the same output files as a result.

### Skimming

Now we want to apply selections to skim out ntuples. There are a few files to look at:

1. `lbl_skimmer_list.py`: input and output paths, you may need to adjust those.
2. `lbl_skimmer_config.py`: here you decide which groups of selections to use (e.g. diphoton, neutral exclusivity, charged exclusivity, etc.). At first, no need to modify anything.
4. `lbl_paths.py`: this is used in every step. In skimming what's especially important is that you'll define the skim name here (could be anything, just something that will explain which selections are applied).
5. `lbl_params.py`: this is where values of different cuts, thresholds, etc. are defined. At first, no need to modify anything.

Once you update your paths and give the skim some name, you can run the skimmer (you can also replace `--local` with `--condor` to run on the grid):

```
python3 submitter.py --app lbl_skimmer --config lbl_skimmer_config.py --files_config lbl_skimmer_list.py --local
```

If everything went well, you should now have skimmed ntuples (so with events passing all selections) in the output path, in a directory with your skim name.

### Histogramming

Now it's time to produce histograms. Have a look at these files:

1. `lbl_histogramer_list.py`: update input/output paths here.
2. `lbl_histogramer_config.py`: this is where histograms and their binning are defined - no need to change anything at first.


With the updated paths, you can run histogramming (again, use `--condor` instead of `--local` to run on the grid):

```
python3 submitter.py --app lbl_histogramer --config lbl_histogramer_config.py --files_config lbl_histogramer_list.py --local
```

After this step, in your output directory, inside of the skim directory, a new one will be created called `histograms` - it contains root files with histograms.
Usually it's good to merge these histograms at this stage:

```
python3 lbl_merge_tea_files.py
```

This will produce a histograms file in the main output directory like this: `merged_skimmed_lblSelections_histograms.root`


### Plotting

The last thing to do is to plot these histograms. Have a look at `lbl_plotter_config.py` to see where the style of the plots, legends, etc. is defined.
Since you only produced LbL MC histograms for now, you will need to comment out all other samples in this config (otherwise it will complain that it cannot
find histograms for collision data, etc.). Then, you can just run:

```
python3 plotter.py lbl_plotter_config.py
```

This will save plots in `plots/your_skim_name` directory, which is parallel to the `bin` directory. 
