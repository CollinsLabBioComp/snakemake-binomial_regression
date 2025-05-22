#!/usr/bin/env snakemake

"""
Binomial regression pipeline
----------------------------

For modifications
"""

__version__ = '0.1.0'


################################################################################
# Helper functions
################################################################################
def format_list_params(arg, pval):
    if type(pval) is list:
        if len(pval) == 0:
            return('')

        return('--{} "{}"'.format(arg, ','.join(pval)))

    # assume string
    if pval == "":
        return('')

    return('--{} "{}"'.format(arg, pval))

################################################################################
rule all:
    input:
        'results/binomial_regression-results.tsv.gz',
        'results/plots/binomial_regression-volcano.png',
        'results/plots/binomial_regression-ma-count_mean.png',
        expand(
            'results/plots/model={model}/plots.done',
            model = list(config['models'])
        )

rule binomial_regression__split:
    """
    Split out files for running regression
    """
    input:
        adata = 'data/data.h5ad'
    output:
        temp(
            expand(
                'results/scatter_data/adata__chunk{j_cur}_{j_total}.h5', 
                j_cur=range(1, config['n_chunks']+1), 
                j_total=config['n_chunks']
            )
        )
    params:
        n_chunks = config['n_chunks'],
        out_base = 'results/scatter_data/adata',
        script = srcdir('scripts/0014-split_andata.py')
    shell:
        '{params.script} '
        '--anndata_file {input.adata} '
        '--number_chunk {params.n_chunks} '
        '--anndata_compression_opts 5 '
        '--output_base {params.out_base}'

################################################################################
# Per-model runs
################################################################################
def get_regression_input(wildcards):
    in_dict = {}
    in_dict['adata'] = 'results/scatter_data/adata__chunk{}_{}.h5'.format(
        wildcards.j_cur,
        wildcards.j_total
    )

    dist = config['models'][wildcards.model]['distribution']
    return(in_dict)
    
rule binomial_regression__run:
    """
    Run binomial regression
    """
    input:
        unpack(get_regression_input)
    output:
        temp('results/model={model}/scatter_analysis/regression_results__chunk{j_cur}_{j_total}.tsv')
    threads: config['n_threads']
    params:
        formula = lambda wildcards: config['models'][wildcards.model]['formula'],
        targ_var = lambda wildcards: config['models'][wildcards.model]['target_variable'],
        cont_covs = format_list_params('continuous_covariates', config['covariates_continuous']),
        disc_covs = format_list_params('factor_covariates', config['covariates_discrete']),
        filt = lambda wildcards: format_list_params('filter', config['models'][wildcards.model]['filter']),
        dist = lambda wildcards: config['models'][wildcards.model]['distribution'],
        n_perms = lambda wildcards: config['models'][wildcards.model]['n_null_permutations'],
        na_fill = lambda wildcards: format_list_params('fill_y_na_value', config['models'][wildcards.model]['fill_y_na_value']),
        script = srcdir('scripts/0024-binomial_regression.R')
    shell:
        '{params.script} '
        '--anndata_file {input.adata} '
        '--formula "{params.formula}" '
        '--model_id "{wildcards.model}" '
        '--target_variable "{params.targ_var}" ' 
        '--distribution "{params.dist}" '
        '--n_permutations {params.n_perms} '
        '{params.cont_covs} '
        '{params.disc_covs} '
        '{params.filt} '
        '{params.na_fill} '
        '--threads {threads} '
        '--verbose '
        '--output_file {output}'

rule binomial_regression__concat:
    """
    Combine the output
    
    Notes:
    - 1 based so first iteration is 1 not 0.
    """
    input:
        expand(
            'results/model={{model}}/scatter_analysis/regression_results__chunk{j_cur}_{j_total}.tsv', 
            j_cur=range(1, config['n_chunks']+1), 
            j_total=config['n_chunks']
        )
    output:
        out = 'results/model={model}/binomial_regression-results.tsv.gz'
    shell:
        """
        set +o pipefail;

        # Get header then do all files
        cat {input[0]} | head -n 1 | gzip > {output.out};

        for file in {input}; do
            cat $file | tail -n +2 | gzip >> {output.out};
        done;
        """

rule binomial_regression__update_results:
    """
    Update results with shrinkage and corrected p-values
    """
    input:
        rez = 'results/model={model}/binomial_regression-results.tsv.gz'
    output:
        flag = 'results/model={model}/update.done'
    params:
        sdist = lambda wildcards: config['models'][wildcards.model]['shrinkage_distribution'],
        script = srcdir('scripts/0026-update_results.R')
    shell:
        """
        mv {input.rez} {input.rez}.tmp;

        {params.script} \
        --results {input.rez}.tmp \
        --shrink_distribution {params.sdist} \
        --output_file {input.rez} \
        --verbose;

        touch {output.flag};
        """

################################################################################
# Final files
################################################################################
rule binomial_regression__merge_final:
    """
    Combine the output
    """
    input:
        flags = expand(
            'results/model={model}/update.done', 
            model = list(config['models'])
        ),
        rez = expand(
            'results/model={model}/binomial_regression-results.tsv.gz', 
            model = list(config['models'])
        )
    output:
        out = 'results/binomial_regression-results.tsv.gz'
    shell:
        """
        set +o pipefail;

        # Get header then do all files
        zcat < {input.rez[0]} | head -n 1 | gzip > {output.out};
        for file in {input.rez}; do
            zcat < $file | tail -n +2 | gzip >> {output.out};
        done;
        """

rule binomial_regression__overview_plots:
    """
    Generate overview plots
    """
    input:
        rez = 'results/binomial_regression-results.tsv.gz'
    output:
        volcano = 'results/plots/binomial_regression-volcano.png',
        mass = 'results/plots/binomial_regression-ma-count_mean.png'
    params:
        out_base = 'results/plots/binomial_regression',
        script = srcdir('scripts/0033-plot_results.R')
    shell:
        '{params.script} '
        '--results {input.rez} '
        '--grouping_col model_id '
        '--output_base {params.out_base} '
        '--verbose;'

rule binomial_regression__association_plots:
    """
    Generate specific association plots
    """
    input:
        adata = 'data/data.h5ad',
        rez = 'results/binomial_regression-results.tsv.gz'
    output:
        flag = 'results/plots/model={model}/plots.done'
    params:
        out_base = lambda wildcards: 'results/plots/model={}'.format(
            wildcards.model
        ),
        sig_col = config['plots']['significance_column'],
        sig_thresh = config['plots']['significance_threshold'],
        sort_col = config['plots']['sort_column'],
        top_n = config['plots']['top_n'],
        bot_n = config['plots']['bottom_n'],
        targ_var = lambda wildcards: config['models'][wildcards.model]['target_variable'],
        cont_covs = format_list_params('continuous_covariates', config['covariates_continuous']),
        disc_covs = format_list_params('factor_covariates', config['covariates_discrete']),
        dist = lambda wildcards: config['models'][wildcards.model]['distribution'],
        filt = lambda wildcards: format_list_params('filter', config['models'][wildcards.model]['filter']),
        script = srcdir('scripts/0029-plot_associations.R')
    shell:
        '{params.script} '
        '--anndata_file {input.adata} '
        '--results {input.rez} '
        '--model_id {wildcards.model} '
        '--significance_column {params.sig_col} '
        '--significance_threshold "{params.sig_thresh}" '
        '--sort_column {params.sort_col} '
        '--plot_top_n_associations {params.top_n} '
        '--plot_bottom_n_associations {params.bot_n} '
        '--target_variable "{params.targ_var}" '
        '--distribution "{params.dist}" '
        '{params.cont_covs} '
        '{params.disc_covs} '
        '{params.filt} '
        '--out_base {params.out_base}; '
        'touch {output.flag}'
        
