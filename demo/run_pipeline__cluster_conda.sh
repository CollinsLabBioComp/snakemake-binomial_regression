#!/bin/sh

## Activate conda environment
source activate snakemake_binomial_regression

## Set any necessary env variables
export SNK_DIR=$(realpath ../.)

snakemake \
     --snakefile ${SNK_DIR}/Snakefile \
     --configfile config_analysis.yml \
     --use-conda \
     --printshellcmds \
     --latency-wait 600 \
     --jobs 999 \
     --cluster-config ${SNK_DIR}/lib/configs/config_cluster_slurm.json \
     --cluster ${SNK_DIR}/lib/wrappers/cluster/slurm.py \
     $1
