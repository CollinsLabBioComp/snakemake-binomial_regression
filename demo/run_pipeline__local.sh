#!/bin/sh

## Activate QTL conda environment
conda activate snakemake_binomial_regression

## Set any necessary env variables
export SNK_DIR="../"

## ADD ANY SPECIFIC ENV NEEDS
snakemake \
     --snakefile ${SNK_DIR}/Snakefile \
     --configfile config_analysis.yml \
     --printshellcmds \
     --cores 1 \
     $1
