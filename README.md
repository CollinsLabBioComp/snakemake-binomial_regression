Snakemake RNA modification regression pipeline
======================

Snakemake pipeline for performing regression in RNA modification data



Quickstart
----------

Quick demo shown below.

1. Set up conda environments shown in `envs`.

```bash
# Starting from this directory
SNK_REPO=$(pwd)

cd envs

# Option 1: direct conda
bash environment_setup.sh

# Option 2: docker (use if running on mac M*)
docker build -t rna_mod_gsis --platform linux/amd64 .
```

2. Run demo

```bash
# Starting from this directory
SNK_REPO=$(pwd)

# set up the directory with the required data dir
cd demo
ln -s data/demo_data/demo_data-random1k.h5ad data/data.h5ad

# Run analysis
bash run_gwas__local_docker.sh
```

Working directory structure
---------------------------

* **Required input files**: in the `data` directory.
* **Final output data**: in the `results` directory.
* **Logs of all jobs**: in the `logs` directory.


Input
---------------------------

AnnData .h5 file with the following attributes:
- Layers contains those specified under left hand side of `formula` configuration. Example: `counts_modified` and `counts_unmodified` in the demo dataset
- Observation slot (.obs) contains columns specified on right hand side of `formula` configuration. Example: `glucose_condition`, `cell_type`, and `library_size` in the demo dataset
- X slot contains total counts of the two test layers. Example: `X` contains `counts_modified` + `counts_unmodified` in the demo dataset