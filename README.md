# Distinguishing empty and cell-containing droplets

This repository contains analysis scripts and simulation code for the manuscript, _Distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data_.
To reproduce the results in the manuscript:

1. Enter `data` and run `download.sh`, which will download publicly available data fromn the 10X Genomics website.
2. Enter `simulations` and run `simrun.R`, which will perform simulations based on the real datasets to assess cell detection methods.
Run `plotsim.R` to recreate the plots in the manuscript for the simulation results.
3. Enter `simulations` and run `realrun.R`, which will apply cell detection methods to the real datasets.
4. Enter `simulations` and run `negcheck.R`, which will examine the p-value distribution reported by _EmptyDrops_ for low-count barcodes.
5. Enter `analysis/pbmc4k` and run `analysis.Rmd`, which will perform the analysis of the PBMC dataset.
Compare to the expected results in `analysis.md`.
6. Enter `analysis/neurons_900` and run `analysis.Rmd`, which will perform the analysis of the 900 brain cell dataset.
Compare to the expected results in `analysis.md`.

You can also examine the LaTeX source of the manuscript in `manuscript`, if you enjoy that sort of thing.
