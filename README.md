# Welcome to CMS UPC-monophoton analysis (Run 2)

## Installation

Create a directory for your project:
- on macOS: `mkdir tea_monophoton.nosync`
- on other platforms: `mkdir tea_monophoton`

For simplicity, the examples below will assume `tea_monophoton`, but adjust the commands if you're on macOS.

To clone the repo, together with all submodules, run:

```bash
cd tea_monophoton
git clone --recurse-submodules git@github.com:jniedzie/tea_monophoton.git .
```

Then, simply run `. tea/build.sh` from the `tea_monophoton` directory.

After updating an existing checkout, run `. tea/build.sh` again to rebuild the
`mono_*` executables and refresh the config links in `bin`. Analysis classes use
the `Mono` prefix, and Python configs use `mono_`. Update any personal scripts
that call the former executables or import the former configs to these names.
LbL sample keys, physics labels, and existing ntuple/histogram paths are retained
so the same input data can still be used.

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

For Condor on lxplus, activate the environment and run the submitter from `bin`:

```bash
source tea/setup.sh
cd bin
python submitter.py --app mono_trigger_selector --config mono_trigger_selector_config.py --files_config mono_trigger_selector_list.py --condor --max_materialize 500 --save_logs
```

The submitter keeps each submission's wrapper and input-file list under
`~/.local/state/tea/condor/` on AFS. `--save_logs` puts stdout, stderr, and the
Condor event log in that submission directory, whose path is printed by the
submitter. Set `TEA_CONDOR_DIR` to another AFS directory to change this location.
Keep the submission directory until its jobs finish. The workers use the EOS
analysis checkout and installed environment, so keep those available as well.

CERN submissions use `condor_submit` with a stable AFS working directory and
transfer the input-file list and any valid VOMS proxy into the worker's scratch
directory. ROOT outputs go to the directories configured in the files config.
This supports deferred job materialization without using `-spool`; an EOS
working directory with `-spool` left factory jobs held before startup. See
[CERN's EOS submission documentation](https://batchdocs.web.cern.ch/troubleshooting/eos.html).
Add `--dry` to prepare and inspect a submission without queuing jobs.

### Renaming

The first step is to format ntuples in a way that `tea` can easily understand, also getting rid of some ambiguity in the original HIForest branch names.
You can see the input and output branch names, together with their types, in `mono_renamer_config.py` - typically there's no need to modify this file.

The input/output paths are defined in `mono_renamer_list.py` - you may need to adjust them, especially the output path.

Finally, have a look at `mono_paths.py` - this is where we decide, among other things, which samples to run on. For instance, to use the example samples,
you can comment out everything except for "lbl" in the `processes` tuple.

When all paths look good, you can run the renaming. It's always a good idea to first try locally (and for this tiny sample it will be fast enough to just run it locally anyway):

```
python3 submitter.py --app mono_renamer --config mono_renamer_config.py --files_config mono_renamer_list.py --local
```

You can now check that new files have been created in the output directory, which contain the `Events` tree with updated branch names.

As an excercise, you can also replace the `--local` flag with `--condor` to see how to run things on the grid. It will automatically schedule 1 job per file and produce the same output files as a result.

### Skimming

Now we want to apply selections to skim out ntuples. There are a few files to look at:

1. `mono_skimmer_list.py`: input and output paths, you may need to adjust those.
2. `mono_skimmer_config.py`: here you decide which groups of selections to use (e.g. diphoton, neutral exclusivity, charged exclusivity, etc.). At first, no need to modify anything.
4. `mono_paths.py`: this is used in every step. In skimming what's especially important is that you'll define the skim name here (could be anything, just something that will explain which selections are applied).
5. `mono_params.py`: this is where values of different cuts, thresholds, etc. are defined. At first, no need to modify anything.

Once you update your paths and give the skim some name, you can run the skimmer (you can also replace `--local` with `--condor` to run on the grid):

```
python3 submitter.py --app mono_skimmer --config mono_skimmer_config.py --files_config mono_skimmer_list.py --local
```

If everything went well, you should now have skimmed ntuples (so with events passing all selections) in the output path, in a directory with your skim name.

### Histogramming

Now it's time to produce histograms. Have a look at these files:

1. `mono_histogramer_list.py`: update input/output paths here.
2. `mono_histogramer_config.py`: this is where histograms and their binning are defined - no need to change anything at first.


With the updated paths, you can run histogramming (again, use `--condor` instead of `--local` to run on the grid):

```
python3 submitter.py --app mono_histogramer --config mono_histogramer_config.py --files_config mono_histogramer_list.py --local
```

After this step, in your output directory, inside of the skim directory, a new one will be created called `histograms` - it contains root files with histograms.
Usually it's good to merge these histograms at this stage:

```
python3 mono_merge_tea_files.py
```

This will produce a histograms file in the main output directory like this: `merged_<skim>_histograms.root`


### Plotting

The last thing to do is to plot these histograms. Have a look at `mono_plotter_config.py` to see where the style of the plots, legends, etc. is defined.
Since you only produced LbL MC histograms for now, you will need to comment out all other samples in this config (otherwise it will complain that it cannot
find histograms for collision data, etc.). Then, you can just run:

```
python3 plotter.py mono_plotter_config.py
```

This will save plots in `plots/your_skim_name` directory, which is parallel to the `bin` directory. 
