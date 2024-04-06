# qml-pub-bias
Estimates t-hurdles, bias adjustments, and local FDRs on the Chen-Zimmermann Dataset for the paper "Do t-stat hurdles need to be raised?"

This version uses QML, replacing several previous versions that use SMM.  QML is much cleaner in this setting due to the huge amount of estimation noise that comes with truncated simulations (unless you do the truncated simulations very carefully).  This lets you run everything without in a few days without a supercomputer.

Having said that the results using QML and SMM are not very different.  The most recent SMM version is here https://github.com/chenandrewy/t-hurdles

# Instructions / description of files

To replicate, run files with prefixes 1-6. It should auto-download any data you need. 

## Code 
All code is in R. Required packages are listed in `O_Settings.r`.

- `0_Settings.r`: Declares libraries, functions, globals, except for the parallel bootstrap
    - `0b_Bootstrap_Parallel.R`: Function for bootstrapping qml estimates. This is a workaround for estimating in parallel on my Windows machine. This way I can source 0_settings.r inside the foreach loop.
- `1_CheckSet.r`: Downloads data from Chen-Zimmermann dataset from web. Also runs estimations on simulated data for checking that this stuff works.
- `2_Baseline.r`: Does the bootstrapped estimation. Saves many Rdata files to disk in `intermediate/`.
- `3_Robustness.r`: Runs bootstrapped estimates for robustness tests. Using 10 cores, takes about 3 hours total. Also outputs Rdata files to `intermediate/`
- `4_Exhibits.r`: Loads output of `2_Baseline.r` and then outputs to `output/`
    - `4b_Robustness_Exhibits`: Outputs robustness exhibits to `output/`
- `5_Demos.r`: Outputs the easy HLZ and intro exhibits to `output/`
- `6_Demos_Other_Evidence.r`: Simulates a model with supplementary evidence, outputs figures to `output/`

## Saved Output / Data
The intermediate data files and output of the last runs of the code are in `intermediate/` and `output/`. Mapping file names to the figures in the paper should be straightforward.