#!/bin/sh

# Check for docker
docker --version

# Get repo to pass to Docker
SNK_REPO=$(pwd)/../.

# Run command inside Docker image -- snakemake only uses singularity natively
mkdir -p results # make results dir to mount 
docker run --name snakemake_binomial_regression \
     --mount type=bind,source="$SNK_REPO",target=/home/container_user/snakemake_binomial_regression \
     --mount type=bind,source="${PWD}",target=/home/container_user/analysis \
     --rm -t \
     snakemake_binomial_regression \
     /bin/bash -c "cd /home/container_user/analysis &&
          source activate snakemake_binomial_regression &&
          snakemake \
               --snakefile /home/container_user/snakemake_binomial_regression/Snakefile \
               --configfile config_analysis.yml \
               --printshellcmds \
               --cores 1 \
               --use-conda \
               $1
          "
