#!/bin/bash
##
## zumis install, including
## R-based dependencies
##
## Izaskun Mallona
## 18th Oct 2021
## GPLv2
##
## ran in sherborne


mkdir -p ~/soft

git clone https://github.com/sdparekh/zUMIs.git
cd zUMIs

# mkdir -p ~/soft/icu
# cd $_
# wget https://github.com/unicode-org/icu/archive/refs/tags/release-58-3.tar.gz
# tar -xzvf release-58-3.tar.gz

# cd ~/soft/icu/icu-release-58-3/icu4c/source
# ./configure
# make


## R based deps
/usr/local/R/R-4.1.0/bin/R -e '
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1")
                                                
cran_pcks <- c("inflection","yaml","shiny","shinythemes","shinyBS","ggplot2","mclust","dplyr","cowplot","Matrix","BiocManager","devtools","stringdist","data.table","stringr","extraDistr")
                                                
install.packages(cran_pcks, repos = "https://cloud.r-project.org/")
bioc_pcks <- c("GenomicRanges","GenomicFeatures","GenomicAlignments","AnnotationDbi","GenomeInfoDb","plyranges")
BiocManager::install(bioc_pcks)
devtools::install_github("VPetukhov/ggrastr")
'
