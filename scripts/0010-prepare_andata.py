#!/usr/bin/env python


__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
from packaging.version import Version
import os
import random
import warnings
import numpy as np
import pandas as pd
import anndata as ad

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def make_anndata(
    count_data,
    metadata_smpls,
    metadata_pos,
    output_file,
    anndata_compression_opts=4
):
    """Merge 10x data.

    Parameters
    ----------
    count_data : pandas.DataFrame
        Description of parameter `count_data`.
    metadata_smpls : pandas.DataFrame
        Description of parameter `metadata_smpls`.
    metadata_pos : pandas.DataFrame
        Description of parameter `metadata_pos`.
    output_file : string
        Description of parameter `output_file`.
    anndata_compression_opts : int
        Anndata gzip compression level.

    Returns
    -------
    output_file : string
        output_file
    """
    # Get compression opts for pandas
    compression_opts = 'gzip'
    if Version(pd.__version__) > Version('1.0.0'):
        compression_opts = dict(
            method='gzip',
            compresslevel=anndata_compression_opts
        )

    # Pivot the count_data to make it wide
    count_data_wide_mod = count_data.pivot_table(
        index='sample_id',     # Rows
        columns='position_id',     # Columns 
        values='modified_counts'
    )
    count_data_wide_unmod = count_data.pivot_table(
        index='sample_id',     # Rows
        columns='position_id',     # Columns
        values='unmodified_counts'
    )
    count_data_wide_sum = count_data_wide_mod + count_data_wide_unmod

    # check no duplicate sample ids
    vals, counts = np.unique(count_data_wide_sum.index, return_counts=True)
    if np.sum(counts > 1):
        raise Exception('Error {} duplicate sample_id:\t{}'.format(
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))
    # print(count_data_wide_sum.index)
    # metadata_smpls = metadata_smpls.transpose()
    # print(metadata_smpls.head())
    # print(metadata_smpls.loc[:, count_data_wide_sum.index])
    # Initialize the AnnData object
    adata = ad.AnnData(
        X=count_data_wide_sum,
        obs=metadata_smpls.loc[count_data_wide_sum.index, :],
        var=metadata_pos.loc[count_data_wide_sum.columns, :]
    )
    # Add in count data as layers
    adata.layers['counts_modified'] = count_data_wide_mod
    adata.layers['counts_unmodified'] = count_data_wide_unmod
    #adata.layers['counts_total'] = count_data_wide_sum
    
    # Record details of X which are total counts
    adata.uns['X_details'] = 'Total counts'

    # Tidy up index
    adata.obs.index.name = None
    adata.obs['sample_id'] = adata.obs.index
    adata.var.index.name = None
    adata.var["position_id"] = adata.var.index
    

    # Possible additional basic filtering on the full dataset.
    # sc.pp.filter_cells(adata, min_genes=200)
    # sc.pp.filter_genes(adata, min_cells=1)
    print('[adata] {} obs, {} vars'.format(
        adata.n_obs,
        adata.n_vars
    ))

    # Save the adata matrix
    # output_file = output_dir 
    adata.write(
        '{}.h5ad'.format(output_file),
        compression='gzip',
        compression_opts=anndata_compression_opts
    )
    # adata.write_csvs(output_file)
    # adata.write_loom(output_file+".loom")
    
    return(output_file)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Creates an AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-cf', '--counts_file',
        action='store',
        dest='counts_file',
        required=True,
        help='Tsv file with the following headers: sample_id,\
            position_id, modified_counts, unmodified_counts.'
    )

    parser.add_argument(
        '-smf', '--sample_metadata_file',
        action='store',
        dest='smf',
        required=True,
        help='File with metadata on samples matching sample_id column in\
            counts_file.'
    )

    parser.add_argument(
        '-pmf', '--position_metadata_file',
        action='store',
        dest='pmf',
        required=True,
        help='File with metadata on samples matching position_id column in\
            counts_file.'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=2,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata',
        help='Basename of output anndata file, assuming output in current \
            working directory. Will have .h5ad appended.\
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
    
    # Load a file of the samples to analyse
    count_data = pd.read_csv(options.counts_file, sep='\t')
    data_check = [
        'sample_id',
        'position_id',
        'modified_counts',
        'unmodified_counts'
    ]
    if not all(k in count_data.columns for k in data_check):
        raise Exception(
            'Error invalid count file. Missing columns.'
        )

    # Load the metadata
    metadata_smpls = pd.read_csv(
        options.smf, 
        sep='\t',
        index_col='sample_id'
    )
    metadata_smpls.columns = metadata_smpls.columns.str.strip(
        ).str.replace(' ', '_').str.lower()

    # Delete the metadata columns that we do not want.
    # if options.mcd != '':
    #     for i in options.mcd.split(','):
    #         if i in metadata.columns:
    #             metadata = metadata.drop(i, axis=1)

    # Load position metadata
    metadata_pos = pd.read_csv(
        options.pmf,
        sep='\t',
        index_col='position_id'
    )

    # Run the merge function.
    out_file = make_anndata(
        count_data,
        metadata_smpls,
        metadata_pos,
        output_file=options.of,
        anndata_compression_opts=options.anndata_compression_opts
    )
    print(out_file)


if __name__ == '__main__':
    main()
