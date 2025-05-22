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
     -it \
     snakemake_binomial_regression \
     bash
