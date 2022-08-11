#!/bin/bash
##
## installs zUMIs on R 4.2x, bioC 3.15
## on sherborne
##
## Izaskun Mallona
## 11 Aug 2022


mkdir -p ~/soft/zumis_alt/
cd $_

git clone https://github.com/sdparekh/zUMIs.git
cd zUMIs



## R based deps for R binaries of the 3.6x series
export R_LIBS=~/R/x86_64-pc-linux-gnu-library/4.2/

## R based deps
/usr/local/R/R-4.2.0/bin/R -e '
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")
                                                
cran_pcks <- c("inflection","yaml","shiny","shinythemes","shinyBS","ggplot2","mclust","dplyr","cowplot","Matrix","BiocManager","devtools","stringdist","data.table","stringr","extraDistr")
                                                
install.packages(cran_pcks, repos = "https://cloud.r-project.org/")
bioc_pcks <- c("GenomicRanges","GenomicFeatures","GenomicAlignments","AnnotationDbi","GenomeInfoDb","plyranges")
BiocManager::install(bioc_pcks)
devtools::install_github("VPetukhov/ggrastr")
'
