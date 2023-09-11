# mimosoid_phylodiversity
Data and analyses for the manuscript entitled, "Phylogenetic diversity and regionalization of root nodule symbiosis." The contents are as follows, summarized by directory structure: 

## Mimosoid_CSVs_ToShare_WGS84
This folder contains diversity calculations represented as CSVs, which were the data products used for model building. For instance, `CANAPE_20km.csv` represents the CANAPE test results. The format is longitude, latitude, value.

## R_figs
Miscellaneous raw figure outputs. These were used for exploratory analysis and to prepare main-text figures.

## environmental_data
Grid-cell-level environmental data extracted for the diversity calculations. All environmental preditors were extracted for all diversity statistics to ensure grid cell matching. For instance, `BIOCLIM_10_CANAPE_20km.csv` represents BIO10 data (Mean Temperature of Warmest Quarter) for the CANAPE results. More predictors were prepared than used, with variable selection performed according to the details noted in the manuscript. The format is value, longitude, latitude.

## sdmMaps_Mimosoid
PNGs representing inferred range maps per species. These were inferred from SDMs as detailed in the main manuscript, with `1` representing presence. One file is included for each species.

## update1
`mimosoid_update.html` is a browsable report on diversity statistics (refer to GitHub documentation for rendering directly from the repository) generated with `mimosoid_update.Rmd`. The subfolder `mimosoidBiodiverse` contains raw output from Biodiverse, including GEOTiffs of the diversity calculations used in the study as well as additional statistics.
