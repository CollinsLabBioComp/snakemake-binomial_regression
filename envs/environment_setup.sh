#!/bin/bash

mamba env create -f environment.yml 

# activate environment and download GMMAT
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate snakemake_binomial_regression

R -e "library(devtools); devtools::install_github('scverse/anndataR')"
R -e "library(devtools); devtools::install_version('car', version = '3.1-3', repos = 'http://cran.us.r-project.org')"
R -e "install.packages('metRology', repos='http://cran.us.r-project.org')"