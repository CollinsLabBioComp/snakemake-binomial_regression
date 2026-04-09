#!/usr/bin/env python

import argparse
import os
import random

import anndata as ad
import numpy as np


seed_value = 0
os.environ["PYTHONHASHSEED"] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)


def main():
    parser = argparse.ArgumentParser(
        description="Sample features from an AnnData file for beta-binomial parameter estimation."
    )
    parser.add_argument(
        "--anndata_file",
        required=True,
        help="Input AnnData file."
    )
    parser.add_argument(
        "--n_mods_sample",
        type=int,
        default=1000,
        help="Number of features to sample."
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Output AnnData file."
    )
    parser.add_argument(
        "--anndata_compression_opts",
        type=int,
        default=4,
        help="Compression level for output AnnData."
    )
    args = parser.parse_args()

    adata = ad.read_h5ad(args.anndata_file)
    n_vars = adata.n_vars
    if n_vars == 0:
        raise ValueError("Input AnnData has no features to sample.")

    n_sample = min(args.n_mods_sample, n_vars)
    sampled_idx = np.sort(np.random.choice(n_vars, size=n_sample, replace=False))
    sampled = adata[:, sampled_idx].copy()
    sampled.write(
        args.output_file,
        compression="gzip",
        compression_opts=args.anndata_compression_opts
    )


if __name__ == "__main__":
    main()
