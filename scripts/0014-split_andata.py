#!/usr/bin/env python

import argparse
import os
import random
import numpy as np
import anndata as ad

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Subset an AnnData object based on chunks.
            """
    )

    parser.add_argument(
        '-af', '--anndata_file',
        action='store',
        dest='af',
        required=True,
        help='Anndata file to be split on var.'
    )

    parser.add_argument(
        '-nchunk', '--number_chunk',
        action='store',
        dest='nchunk',
        default=2,
        type=int,
        help='Total number chunks to split the data into.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ob', '--output_base',
        action='store',
        dest='output_base',
        default='chunk',
        help='Basename of output anndata file, assuming output in current \
            working directory. Consider adding .h5ad to the end.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--anndata_compression_opts',
        action='store',
        dest='anndata_compression_opts',
        default=4,
        type=int,
        help='Compression level in anndata. A larger value decreases disk \
            space requirements at the cost of compression time. \
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the file details, without reading the full file into memory
    adata_stream = ad.read_h5ad(options.af)

    # get different chunks
    chunks = np.array_split(np.arange(len(adata_stream.var)), options.nchunk)
    for i in range(options.nchunk):
        # Read in the chunk
        adata = adata_stream[:, chunks[i]]
        
        # Save the adata matrix
        # output_file = output_dir 
        adata.write(
            '{}__chunk{}_{}.h5'.format(options.output_base, str(i+1), options.nchunk),
            compression='gzip',
            compression_opts=options.anndata_compression_opts
        )

    # close file
    adata_stream.file.close()

if __name__ == '__main__':
    main()
