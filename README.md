# Distinguishing empty and cell-containing droplets

## Overview

This repository contains analysis scripts and simulation code for the manuscript **Distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data**
by [Lun _et al._ (2018)](https://doi.org/10.1101/234872).
It also contains the manuscript files themselves, which can be compiled with `pdflatex` if you're into that sort of thing.

## Downloading files

Enter `data` and run:

- `download_tenx.sh`, which will download publicly available data from the 10X Genomics website.
- `download_placenta.sh`, which will download unprocessed 10X data from [the Vento-Tormo _et al._ study](https://doi.org/10.1038/s41586-018-0698-6).

## Simulations

Enter `simulations` and run:

- `simrun.R`, which will perform simulations based on the real datasets to assess cell detection methods.
- `plotsim.R`, to recreate the plots in the manuscript for the simulation results.

## Real data

Enter `real` and run:

- `realrun.R`, which will apply cell detection methods to the real datasets.
- `negcheck.R`, which will examine the p-value distribution reported for low-count barcodes.

Each subdirectory of `analysis` contains self-contained analysis files for each data set.
Run:

- `analysis.Rmd`, which will perform the analysis of each dataset.
- `plot_maker.R`, which will create the plots used in the manuscript.

