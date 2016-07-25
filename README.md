# Description

This directory contains the R/RMarkdown code and associated files necessary to reproduce the paper by Marchand *et al.* "Reserves and trade jointly determine exposure to food supply shocks" (currently a preprint submitted to Environmental Research Letters).


# List of files

- `cereals-network-article.Rmd`: RMarkdown script for the article
- `cereals-network-supplementary.Rmd`: RMarkdown script for the supplementary materials
- `cereals-network-funcs.R`: R script containing functions for model simulation and summary output generation, used by the scripts above
- `model-flow-chart.pdf`: non-script generated flow chart of the simulation model
- `fig_dCext_net.png`: manually edited version of a script-generated figure (adding country names key)
- `ref.bib`: BibTeX references
- `institute-of-physics-harvard.csl`: style file to format citations
- `ciso3.txt`: reference list of ISO3 country codes for the countries in the data files

## in `data/` directory
- `get_cereals_prod_trade.R`: code used to obtain cereals production and trade from public FAOSTAT data
- `cereals_prod_trade.RData`: output file from the script above used in the article, based on FAOSTAT data acquired in January 2016
- `get_cereals_stocks.R`: code used to obtain cereals stocks from public USDA-PSD data
- `cereals_stocks.RData`: output file from the script above used in the article, based USDA-PSD data acquired in October 2015
- `country_list.csv`: list of countries in data with FAOSTAT codes and population information
- `crop_list.tsv`: list of cereal commodities considered in this article, with FAOSTAT codes and supplementary information (e.g. kcal/ton conversion factors)


# R packages dependencies

- Data manipulation: `tidyr`, `dplyr`, `reshape2`, `stringr`
- Visualization: `ggplot2`, `maptools`, `cowplot`
- LaTeX/PDF generation: `knitr`, `stringr`
- Others: `igraph`, `FAOSTAT`


# Acknowledgments

This work was supported by the National Socio-Environmental Synthesis Center (SESYNC) under funding received from the National Science Foundation DBI-1052875.
