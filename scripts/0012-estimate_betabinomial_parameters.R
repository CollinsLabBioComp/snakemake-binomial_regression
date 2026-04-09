#!/usr/bin/env Rscript

SCRIPT_NAME <- "estimate_betabinomial_parameters.R"

old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(anndataR))
suppressPackageStartupMessages(library(gamlss))

set.seed(0)

fit_betabinomial <- function(df, formula_string) {
    if (!requireNamespace("gamlss.dist", quietly = TRUE)) {
        stop("Package 'gamlss.dist' is required for beta-binomial estimation.")
    }
    model <- gamlss::gamlss(
        as.formula(formula_string),
        family = gamlss.dist::BB,
        data = df,
        control = gamlss::gamlss.control(
            n.cyc = 1e5,
            trace = FALSE
        )
    )
    if (!isTRUE(model$converged)) {
        stop("Model did not converge.")
    }
    return(model)
}

command_line_interface <- function() {
    option_list <- list(
        optparse::make_option(
            c("--anndata_file"),
            type = "character",
            help = "AnnData file with counts_modified and counts_unmodified layers."
        ),
        optparse::make_option(
            c("--formula"),
            type = "character",
            default = "cbind(counts_modified,counts_unmodified) ~ 1",
            help = "Model formula. Only the response is used; the RHS is replaced with ~ 1."
        ),
        optparse::make_option(
            c("--filter"),
            type = "character",
            default = "",
            help = "Optional filter to apply to obs before estimation."
        ),
        optparse::make_option(
            c("--fill_y_na_value"),
            type = "numeric",
            default = NA,
            help = "Value to fill NA counts with before estimation."
        ),
        optparse::make_option(
            c("--per_sample"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = "Estimate separate parameters for each sample."
        ),
        optparse::make_option(
            c("--output_mu_file"),
            type = "character",
            default = "estimated_mu.txt",
            help = "Output file for estimated mu."
        ),
        optparse::make_option(
            c("--output_sigma_file"),
            type = "character",
            default = "estimated_sigma.txt",
            help = "Output file for estimated sigma."
        ),
        optparse::make_option(
            c("--verbose"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = "Verbose mode."
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = option_list,
        description = "Estimate beta-binomial parameters from an AnnData file."
    )

    arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }

    response <- strsplit(param[["formula"]], split = "~", fixed = TRUE)[[1]][1]
    param[["formula"]] <- sprintf("%s ~ 1", trimws(response))

    adata <- anndataR::read_h5ad(
        param[["anndata_file"]],
        as = "InMemoryAnnData"
    )

    if (!("counts_modified" %in% names(adata$layers)) ||
        !("counts_unmodified" %in% names(adata$layers))) {
        stop("AnnData must contain counts_modified and counts_unmodified layers.")
    }

    if (!is.na(param[["fill_y_na_value"]])) {
        adata$X[is.na(adata$X)] <- param[["fill_y_na_value"]]
        for (layer_name in names(adata$layers)) {
            adata$layers[[layer_name]][is.na(adata$layers[[layer_name]])] <- param[["fill_y_na_value"]]
        }
    }

    if (!is.na(param[["filter"]]) && (param[["filter"]] != "")) {
        obs_tmp <- dplyr::filter(
            adata$obs,
            !! rlang::parse_expr(param[["filter"]])
        )
        obs_ret <- rownames(obs_tmp)
        obs_ix <- match(obs_ret, rownames(adata$obs))
        adata <- anndataR::AnnData(
            X = adata$X[obs_ix, ],
            obs = adata$obs[obs_ret, ],
            var = adata$var,
            layers = lapply(adata$layers, function(x) x[obs_ix, ]),
            obsm = lapply(adata$obsm, function(x) x[obs_ret, ]),
            varm = adata$varm,
            obsp = lapply(adata$obsp, function(x) x[obs_ret, ]),
            varp = adata$varp,
            uns = adata$uns
        )
    }

    sample_ids <- rownames(adata$obs)
    mod_mat <- as.matrix(adata$layers[["counts_modified"]])
    unmod_mat <- as.matrix(adata$layers[["counts_unmodified"]])

    if (param[["verbose"]]) {
        message("Estimating beta-binomial parameters [", SCRIPT_NAME, "]...")
        message("\tper_sample = ", param[["per_sample"]])
        message("\tn_obs = ", nrow(mod_mat))
        message("\tn_vars = ", ncol(mod_mat))
    }

    mu_df <- data.frame(sample_id = character(0), val = numeric(0))
    sigma_df <- data.frame(sample_id = character(0), val = numeric(0))

    build_counts_df <- function(mod_counts, unmod_counts) {
        df <- data.frame(
            counts_modified = as.numeric(mod_counts),
            counts_unmodified = as.numeric(unmod_counts)
        )
        total_counts <- df[["counts_modified"]] + df[["counts_unmodified"]]
        df <- df[stats::complete.cases(df) & !is.na(total_counts), , drop = FALSE]
        return(df)
    }

    if (param[["per_sample"]]) {
        for (i in seq_along(sample_ids)) {
            df <- build_counts_df(mod_mat[i, ], unmod_mat[i, ])
            if (nrow(df) == 0) {
                next
            }
            model <- fit_betabinomial(df, param[["formula"]])
            mu_df <- rbind(mu_df, data.frame(sample_id = sample_ids[i], val = unname(fitted(model, what = "mu"))[1]))
            sigma_df <- rbind(sigma_df, data.frame(sample_id = sample_ids[i], val = unname(fitted(model, what = "sigma"))[1]))
        }
    } else {
        df <- build_counts_df(mod_mat, unmod_mat)
        if (nrow(df) == 0) {
            stop("No complete rows available to estimate beta-binomial parameters.")
        }
        model <- fit_betabinomial(df, param[["formula"]])
        mu_df <- rbind(mu_df, data.frame(sample_id = "global", val = unname(fitted(model, what = "mu"))[1]))
        sigma_df <- rbind(sigma_df, data.frame(sample_id = "global", val = unname(fitted(model, what = "sigma"))[1]))
    }

    if ((nrow(mu_df) == 0) || (nrow(sigma_df) == 0)) {
        stop("Failed to estimate beta-binomial parameters.")
    }

    write.table(
        mu_df,
        file = param[["output_mu_file"]],
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
    write.table(
        sigma_df,
        file = param[["output_sigma_file"]],
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
}

main <- function() {
    run_time <- system.time(command_line_interface())
    message(
        "Analysis execution time [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]] / 3600,
        " hours."
    )
}

main()
