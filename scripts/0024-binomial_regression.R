#!/usr/bin/env Rscript

SCRIPT_NAME <- "binomial_regression.R"

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(anndataR))
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(fitdistrplus))
suppressPackageStartupMessages(library(metRology))
suppressPackageStartupMessages(library(gamlss))

# for parallel processing
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(iterators))


set.seed(0)

split_into_chunks <- function(x, n) {
  split(x, cut(seq_along(x), n, labels = FALSE))
}

#' Inverse normalize.
#'
#' @param x List.
#'     Vector data to inverse normalize
#' 
#' @return List.
#'     Rank inverse normalized vector
#'
#' @export
rank_inverse_normal <- function(x) {
    return(qnorm((rank(x, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(x))))
}

#' Combines dataframes in a list
#'
#' @importFrom data.table rbindlist
rbindlist_df <- function(...) {
    return(data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE))
}

rbind_fillna <- function(...) {
    return(data.table::rbindlist(
        list(...),
        use.names = TRUE,
        fill = TRUE
    ))
}

#' Original code from roland-hochmuth/gamlssdiag. Modified to fit the needs of this script.
#' 
#' Evaluates Cook's Distance.
#'
#' @note Although lm/glm stores the original data in the model at m$data, that is not the case in GAMLSS.
#' Therefore, unlike the function cooks.distance which only requires a model parameter, this function requires
#' the formula and data.
#'
#' @return A vector of distances
#'
#' @references https://www.ime.usp.br/~abe/lista/pdf1USQwcGBX1.pdf
#' @references Robust Diagnostic Regression Analysis, Anthony Atkinson and Marco Riani
#' @references http://www.gamlss.com/
#'
#' @export
cooksd <- function(model, data, param) {
  coefficients <- coef(model)
  pred <- exp(predict(model, newdata=data, data=data, what=c("mu")))
  num_coefficients <- length(coefficients) + 1
  num_observations <- nrow(data)
  resids <- residuals(model, type = "weighted")
  # Evaluate the residual mean square estimate of the variance using 
  # Formula 2.12 in Robust Diagnostic Regression Analysis.
  s2 <- sum(resids^2)/(num_observations - num_coefficients)
  d <- vector(mode = "double", length = num_observations)

  for (i in seq_along(d)) {
    # Create the leave-one-out (loo) dataframe.
    loo_data <- dplyr::filter(data, !dplyr::row_number() %in% i)
    loo_m <- regress(loo_data, param)
    loo_p <- exp(predict(loo_m, newdata=data, data=loo_data, what=c("mu")))
    # Evaluate Cook's Distance using Formula 2.41 in Robust Diagnostic Regression Analysis.
    d[i] <- sum((loo_p - pred)^2) / (num_coefficients*s2)
  }
  return(d)
}

empirical_distribution <- function(
        permuted_statistics,
        observed_statistic,
        distribution = "t.scaled",
        custom_params = list()
    ) {
    if (distribution == "t.scaled") {
        if (!("df" %in% names(custom_params))) {
            custom_params[["df"]] <- 1
        }
        custom_params[["mean"]] <- mean(permuted_statistics, na.rm = T)
        custom_params[["sd"]] <- sd(permuted_statistics, na.rm = T)
        fit_dist <- fitdistrplus::fitdist(
            permuted_statistics, 
            distribution,
            start = custom_params
        )
        p_lower <- metRology::pt.scaled(
            observed_statistic,
            df = fit_dist$estimate["df"],
            mean = fit_dist$estimate["mean"],
            sd = fit_dist$estimate["sd"],
            lower.tail = TRUE
        )
    } else if (distribution == "norm") {
        fit_dist <- fitdistrplus::fitdist(
            permuted_statistics, 
            distribution
        )
        p_lower <- pnorm(
            observed_statistic,
            mean = fit_dist$estimate["mean"],
            sd = fit_dist$estimate["sd"],
            lower.tail = TRUE
        )
    } 
    # NOTE: computing p-value this way, by taking min
    # handles for case where observed_statistic is 
    # less than the mean
    empirical_p_value <- 2 * min(p_lower, 1 - p_lower)
    return(empirical_p_value)
}

run_permutations <- function(
        df_reg,
        param,
        original_result,
        n_permutations = 1e5,
        early_stop_threshold = 0
    ) {
    # Store p-values from permutations
    # permuted_p_values <- rep(NA, n_permutations)
    # observed_p_value <- original_result[["p_value"]]
    permuted_statistics <- rep(NA, n_permutations)
    observed_statistic <- original_result[["test_statistic"]]
    count_more_extreme <- 0  # Counter for permutations with p-value below the observed
    df_reg[['target_var_original']] <- df_reg[[param[["target_variable"]]]]

    for (i in 1:n_permutations) {
        # Shuffle the target variable
        df_reg[[param[["target_variable"]]]] <- sample(
            df_reg[[param[["target_variable"]]]]
        )
        
        # Fit the model with permuted target variable
        permuted_model <- tryCatch({
            regress(df_reg, param)
        }, error = function(e) {
            #message("Error: when fitting permuted model: ", e$message, "\n")
            NULL  # If the model fails, return NULL
        })

        if (!is.null(permuted_model)) {
            # Extract p-value for the target variable
            coef_table <- summary(permuted_model)$coefficients
            target_row <- grep(
                pattern = param[["target_variable"]],
                x = rownames(coef_table),
                fixed = TRUE
            )
            if (length(target_row) == 1) {
                # permuted_p_values[i] <- coef_table[
                #     target_row,
                #     grepl(
                #         pattern = "Pr(",
                #         x = names(coef_table[target_row,]),
                #         fixed = TRUE
                #     )
                # ]
                permuted_statistics[i] <- coef_table[
                    target_row,
                    param[["test_statistic"]]
                ]
            }
        }
        
        # Early stopping criterion
        if (early_stop_threshold > 0) {
            empirical_p_value <- (
                sum(abs(permuted_statistics) >= abs(observed_statistic), na.rm = T) + 1
            ) / (sum(!is.na(permuted_statistics), na.rm = T) + 1)
            if (i >= early_stop_threshold && empirical_p_value >= 0.75) {
                print(sum(abs(permuted_statistics) >= abs(observed_statistic), na.rm = T))
                print(sum(!is.na(permuted_statistics), na.rm = T))
                message(
                    "Stopping permutations early after ", 
                    early_stop_threshold,
                    " permutations. Empirical p-value: ",
                    empirical_p_value
                )
                break
            }
        }
    }
    empirical_p_value <- NA
    empirical_p_value_dist <- NA
    if (!all(is.na(permuted_statistics))) {
        #permuted_p_values <- permuted_p_values[!is.na(permuted_p_values)]
        permuted_statistics <- permuted_statistics[!is.na(permuted_statistics)]
        if (length(permuted_statistics) >= 3) {
            # empirical_p_value <- (
            #     sum(permuted_p_values <= observed_p_value, na.rm = T) + 1
            # ) / (length(permuted_p_values) + 1)
            empirical_p_value <- (
                sum(abs(permuted_statistics) >= abs(observed_statistic), na.rm = T) + 1
            ) / (length(permuted_statistics) + 1)
            # If the empirical p-value is 1, then the method to estimate
            # the p-value from the distribution is almost always off, so handle this
            # case. 
            if (empirical_p_value == 1) {
                empirical_p_value_dist <- 1
            } else {
                distribution <- "norm"
                if (param[["test_statistic"]] == "t value") {
                    distribution <- "t.scaled"
                }
                empirical_p_value_dist <- tryCatch({
                    empirical_distribution(
                        permuted_statistics,
                        observed_statistic,
                        distribution,
                        custom_params = list(
                            df = max(1, original_result[["degrees_freedom"]])
                        )
                    )
                }, error = function(e) {
                    message("Error: when fitting null model for empirical p-values: ", e$message, "\n")
                    NA  # If the model fails, return NA
                })
            }
        }
    }
    
    # Reset the target variable back to the original
    df_reg[[param[["target_variable"]]]] <- df_reg[['target_var_original']]
    df_reg[['target_var_original']] <- NULL

    #cat(empirical_p_value, empirical_p_value_dist, "\n")
    return_list <- list(
        "p_value_empirical" = empirical_p_value,
        "p_value_empirical_distribution" = empirical_p_value_dist,
        "n_permutations" = length(permuted_statistics)
    )
    return(return_list)
}

#' Perform regression
#'
#' @param df data.frame.
#'     Data frame with regression information
#'     
#' @param p List.
#'     List of paramters
#' 
#' @return List.
#'     Rank inverse normalized vector
#'
#' @export
regress <- function(df, p) {
    dist_family <- binomial
    if (p[["distribution"]] == 'quasibinomial') {
        dist_family <- quasibinomial
    }
    
    if (!grepl("|", p[["formula"]], fixed = T)) {
        model <- glm(
            p[["formula"]],
            family = dist_family,
            data = df,
            control = list(
                #epsilon = 1e-10,
                maxit = 1e5
            )
        )
    } else {
        # LME4 does not support quasibinomial
        # https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#fitting-models-with-overdispersion
        model <- lme4::glmer(
            p[["formula"]],
            family = dist_family,
            data = df,
            verbose = FALSE,
            control = lme4::glmerControl(
                #tolPwrss = 1e-10,
                optCtrl = list(iter.max = 1e5, eval.max = 1e5)
            )
        )
    }
    return(model)
}

# #' Combines dataframes in a list
# #'
# #' @importFrom data.table rbindlist
# rbindlist_dflist <- function(...) {
#     df_coloc <- data.table::rbindlist(
#         lapply(list(...), FUN = function(x) {return(x[["df_coloc"]])}),
#         use.names = TRUE,
#         fill = TRUE
#     )
#     df_harmonized_ss <- data.table::rbindlist(
#         lapply(
#             list(...),
#             FUN = function(x) {return(x[["df_harmonized_ss"]])}
#         ),
#         use.names = TRUE,
#         fill = TRUE
#     )

#     return(list(
#         "df_coloc" = df_coloc,
#         "df_harmonized_ss" = df_harmonized_ss
#     ))
# }


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("--anndata_file"),
            type = "character",
            help = paste0(
                "Anndata file with the following layers: counts_modified, counts_unmodified",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--formula"),
            type = "character",
            default = "cbind(counts_unmodified,counts_modified) ~ glucose_condition + library_size",
            help = paste0(
                "Formula.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--target_variable"),
            type = "character",
            default = "glucose_condition",
            help = paste0(
                "Target variable.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--continuous_covariates"),
            type = "character",
            default = "library_size",
            help = paste0(
                "Comma separated list of continuous covariates.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--factor_covariates"),
            type = "character",
            default = "glucose_condition",
            help = paste0(
                "Comma separated list of discrete covariates.",
                " [default %default]"
            )
        ),
        
        optparse::make_option(c("--model_id"),
            type = "character",
            default = "model",
            help = paste0(
              "ID for model.",
              " [default %default]"
            )
        ),
        
        optparse::make_option(c("--filter"),
            type = "character",
            default = "",
            help = paste0(
              "Filter for anndata.",
              " [default %default]"
            )
        ),
        
        optparse::make_option(c("--distribution"),
            type = "character",
            default = 'binomial',
            help = paste0(
              "Distribution to use to model data. Should be one of:",
              "`binomial` or `quasibinomial`",
              " [default %default]"
            )
        ),

        optparse::make_option(c("--n_permutations"),
            type = "numeric",
            default = 0,
            help = paste0(
                "Number permutations. If 0, then no permutations run.",
                " [default %default]"
            )
        ),


        optparse::make_option(c("--output_file"),
            type = "character",
            default = "binomial_regression-results.tsv",
            help = paste0(
                "Output file.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--threads"),
            type = "numeric",
            default = 1,
            help = paste0(
                "Number of threads to use.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--fill_y_na_value"),
            type = "numeric",
            default = NA,
            help = paste0(
                "Value to fill NAs in y variable with. If parameter set to NA does nothing.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--verbose"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "Verbose mode (write extra info to std.err).",
                " [default %default]"
            )
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Runs binomial regression for anndata files."
        )
    )

    # a hack to fix a bug in optparse that won"t let you use positional args
    # if you also have non-boolean optional args:
    getOptionStrings <- function(parserObj) {
        optionStrings <- character()
        for (item in parserObj@options) {
            optionStrings <- append(optionStrings,
                                    c(item@short_flag, item@long_flag))
        }
        optionStrings
    }
    optStrings <- getOptionStrings(parser)
    arguments <- optparse::parse_args(parser, positional_arguments = TRUE)

    # read in the parameters
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }
    for (i in c("factor_covariates", "continuous_covariates")) {
        param[[i]] <- strsplit(param[[i]], split = ",")[[1]]
    }

    # Add in the test statistic column based on the regression model
    # NOTE: this is the output in coef_table
    if (param[["distribution"]] == 'binomial') {
        param[["test_statistic"]] <- "z value"
    } else if (param[["distribution"]] == 'quasibinomial') {
        param[["test_statistic"]] <- "t value"
    } else {
        stop("ERROR: Invalid distribution.")
    }

    # NOTE: in theory one could only read in chunk, but anndataR 
    # does not seem to support a view
    #chunk_indices <- split_into_chunks(adata$var$position_id, 3)
    df_results <- run_analysis(
      param,
      fill_y_na_value = param[["fill_y_na_value"]],
      inverse_normalize_covariates = FALSE
    )
    
    write.table(
      df_results,
      file = param[['output_file']],
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t",
      na = ""
    )

    return(df_results)
}


#' Run analysis function
#'
#' @importFrom data.table fread
#' @importFrom parallel detectCores
#' @importFrom parallel makeForkCluster
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom doMC registerDoMC
#' @importFrom iterators isplit
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
run_analysis <- function(
    param,
    fill_y_na_value = NA,
    inverse_normalize_covariates = FALSE,
    parallel_type = NA,
    size_factor_column = "library_size"
) {

    # Read in the data and fill in NAs if needed
    adata <- anndataR::read_h5ad(
        param[["anndata_file"]] ,
        to = "InMemoryAnnData"
    )
    if (!is.na(fill_y_na_value)) {
        adata$X[is.na(adata$X)] <- fill_y_na_value
        for (i in names(adata$layers)) {
            adata$layers[[i]][is.na(adata$layers[[i]])] <- fill_y_na_value
        }
    }

    # Apply filter that may exist
    if (!is.na(param[['filter']]) & (param[['filter']] != "")) {
        obs_ret <- adata$obs %>%
            dplyr::filter(!! rlang::parse_expr(param[['filter']])) %>%
            rownames(.)
        obx_ix <- match(obs_ret, rownames(adata$obs))

        # re-create anndata with subset
        adata <- anndataR::AnnData(
            X = adata$X[obx_ix, ],
            obs = adata$obs[obs_ret, ],
            var = adata$var,
            layers = lapply(adata$layers, function(x) x[obx_ix, ]),
            obsm = lapply(adata$obsm, function(x) x[obs_ret, ]),
            varm = adata$varm,
            obsp = lapply(adata$obsp, function(x) x[obs_ret, ]),
            varp = adata$varp,
            uns = adata$uns
        )
    }

    # set the parallel backend: fork, snow, domc
    if (is.na(parallel_type)) {
        param[["parallel_type"]] <- "domc"
    } else {
        param[["parallel_type"]] <- parallel_type
    }

    # print out parameters
    if (param[["verbose"]]) {
        message(paste0("Parameters [", SCRIPT_NAME, "]:"))
        for (i in names(param)) {
            pstr <- param[[i]]
            if (is.vector(param[[i]]) | is.list(param[[i]])) {
                pstr <- paste(param[[i]], collapse = ', ')
            }
            
            message("\t", i, " = ", pstr)
        }
    }

    # for parallel processing with foreach
    # more info here: https://ljdursi.github.io/beyond-single-core-R
    time_init <- Sys.time()
    if (param[["threads"]] > parallel::detectCores()) {
        warning(paste0(
            "User requested more threads than available. Setting number of",
            " threads to doParallel::detectCores() - 1, that is",
            " ", parallel::detectCores() - 1, "."
        ))
        param[["threads"]] <- parallel::detectCores() - 1
    }
    if (param[["parallel_type"]] == "fork") {
        sink("/dev/null") # use sink to grab annoying output when outfile = ""
            # (v1) Forking = copy the R session in its current state. Fast
            # because only copies objects if they are modified; however it takes
            # much more memory as files are copied
            cl <- parallel::makeForkCluster(param[["threads"]], outfile = "")
        sink()
        doParallel::registerDoParallel(cl)
    } else if (param[["parallel_type"]] == "snow") {
        # since snow could be on any machine, set outfile to /dev/null rather
        # than the parent process (unlink fork)
        # (v2) Spawn new subprocess: benefits can be done on any machine,
        # new process will not have unneeded data
        if (param[["verbose"]]) {
            message(paste0(
                "[", SCRIPT_NAME, "]:\t",
                " Running parallel with snow. Any std.out or std.err (apart",
                " from a failed exit) will be lost."
            ))
        }
        cl <- parallel::makeCluster(param[["threads"]], outfile = "/dev/null")
        parallel::clusterExport(cl, c(
            "SCRIPT_NAME", "param", "adata", "rbindlist_dflist"
        ))
        doSNOW::registerDoSNOW(cl)
    } else if (param[["parallel_type"]] == "domc") {
        # (v3) Shared memory process
        # doesn't use persistent workers like the snow package or ths
        # now-derived functions in the parallel package
        doMC::registerDoMC(param[["threads"]])
    } else {
        stop(paste0("[", SCRIPT_NAME, "]:\tERROR invalid parallel_type."))
    }
    # iterate over each row
    #
    # note that if using foreach and combining the results into a dataframe,
    # errors will not stop execution, but will not necessarily be returned to
    # the user. Will show up as:
    # Error in { : task 1 failed - "replacement has 1 row, data has 0"
    #
    df_final <- foreach::foreach(
        var_row = iterators::iter(adata$var, by = "row"), # iterator row
        #var_row = iterators::isplit(adata, adata$$var$position_id), # iterator split
        .combine = rbindlist_df,
        .inorder = FALSE,
        .multicombine = TRUE,
        .errorhandling = "stop"
    ) %dopar% {
        #message(paste0("iter:\t", row.names(var_row), "\n"))

        # Vector to filter the anndata
        filt_row <- row.names(adata$var) == row.names(var_row)

        # Data frame of variables for the regression
        # Assumes adata$obs + adata$layers contains all of the variables
        # needed for the regression
        df_reg <- adata$obs
        for (i in names(adata$layers)) {
           df_reg[[i]] <- adata$layers[[i]][,filt_row]
        }
        # Assume that X is total counts
        df_reg[["total_counts"]] <- adata$X[,filt_row] 
        #df_reg[["total_counts"]] <- df_reg[["counts_unmodified"]] + df_reg[["counts_modified"]]

        # Remove NAs in Y
        df_reg <- subset(df_reg, !is.na(df_reg$total_counts))

        # Start to fill out the return dataframe
        df_return <- var_row
        df_return$model_id <- param[['model_id']]
        df_return$formula <- param[["formula"]]
        df_return$target_variable <- param[["target_variable"]]
        df_return$n <- nrow(df_reg)
        df_return$count_mean <- mean(df_reg[["total_counts"]], na.rm = T)
        df_return$count_median <- median(df_reg[["total_counts"]], na.rm = T)
        df_return$count_max <- max(df_reg[["total_counts"]], na.rm = T)
        df_return$count_min <- min(df_reg[["total_counts"]], na.rm = T)
        df_return$count_pct_grtr_eq_mean <- sum(
            df_reg[["total_counts"]] >= df_return$count_mean, na.rm = T
        ) / df_return$n
        df_return$count_pct_grtr_eq_median <- sum(
            df_reg[["total_counts"]] >= df_return$count_median, na.rm = T
        ) / df_return$n
        df_return$n_zero_samples <- sum(df_reg[["total_counts"]] == 0)
        if (size_factor_column %in% colnames(df_reg)) {
            total_counts_normalized <- df_reg[["total_counts"]] / df_reg[[size_factor_column]]
            df_return$count_normalized_mean <- mean(total_counts_normalized, na.rm = T)
            df_return$count_normalized_median <- median(total_counts_normalized, na.rm = T)
            df_return$count_normalized_max <- max(total_counts_normalized, na.rm = T)
            df_return$count_normalized_min <- min(total_counts_normalized, na.rm = T)
            df_return$count_normalized_pct_grtr_eq_mean <- sum(
                total_counts_normalized >= df_return$count_normalized_mean, na.rm = T
            ) / df_return$n
            df_return$count_normalized_pct_grtr_eq_median <- sum(
                total_counts_normalized >= df_return$count_normalized_median, na.rm = T
            ) / df_return$n
        }
        df_return$model_message <- NA

        # Check to make sure the target variable has more than one unique value
        target_variable_pass <- TRUE
        if (length(unique(df_reg[[param[["target_variable"]]]])) <= 1) {
            #message(paste0("[", SCRIPT_NAME, "]:\tERROR target variable is unique."))
            target_variable_pass <- FALSE
        } 

        # Process any variables that need processing for the regression
        for (f in param[["factor_covariates"]]) {
            if (grepl(pattern = '::', f, fixed=T)) {
                fact_var <- strsplit(f, split='::', fixed=T)[[1]][1]
                fact_ref <- strsplit(f, split='::', fixed=T)[[1]][2]
                df_reg[[fact_var]] <- factor(
                  df_reg[[fact_var]],
                  levels = c(
                    fact_ref,
                    setdiff(unique(df_reg[[fact_var]]), fact_ref)
                  )
                )
            } else {
                df_reg[[f]] <- as.factor(df_reg[[f]])
            }
        }
        for (i in param[["continuous_covariates"]]) {
            if (i != "library_size") {
                df_reg[[i]] <- as.numeric(scale(
                    df_reg[[i]], 
                    center = TRUE, 
                    scale = TRUE
                ))
                if (inverse_normalize_covariates == TRUE) {
                    df_reg[[i]] <- rank_inverse_normal(df_reg[[i]])
                }
            }
        }

        # Fit the model
        model <- NULL
        if (df_return$n > 3 & target_variable_pass) {
            model <- tryCatch({
                model <- regress(df_reg, param)
                model
            }, warning = function(w) {
                df_return$model_message <<- w$message[[1]]
                if (grepl("algorithm did not converge", w$message)) {
                    message("Warning: Model did not converge.")
                    NULL  # Return NULL or a different value if needed
                }
                #stop("Warning: when fitting model: ", w$message, "\n")
                model # Return the model
            }, error = function(e) {
                df_return$model_message <<- e$message[[1]]
                # If an error occurs, handle it here
                message("Error: when fitting model: ", e$message, "\n")
                NULL  # Return NULL or any other appropriate value
            })
        } else {
            if (df_return$n <= 3) {
                df_return$model_message <- "Model not fit: n<=3."
            } else {
                df_return$model_message <- "Model not fit: target variable is unique."
            }
        }

        cols_to_add <- c(
            "p_value",
            "estimate",
            "std_error",
            "test_statistic",
            "degrees_freedom",
            "cooks_distance_max",
            "cooks_distance_pvalue_min",
            "cooks_distance_pvalue_n_less_0pt05",
            "cooks_distance_pvalue_n_less_0pt01",
            "cooks_distance_pvalue_n_less_0pt01",
            "dispersion",
            #"dispersion_pvalue",
            "vif_target_variable",
            "vif_max", 
            "p_value_empirical",
            "p_value_empirical_distribution",
            "n_permutations"
        )
        for (c in cols_to_add) {
            df_return[[c]] <- NA
        }
        if (!is.null(model)) {
            coef_table <- summary(model)$coefficients
            target_row <- grep(
                pattern = param[["target_variable"]],
                x = rownames(coef_table),
                fixed = TRUE
            )
            if (length(target_row) > 1) {
                message(target_row)
                message(coef_table)
                message(coef_table[target_row,])
                stop(paste0("[", SCRIPT_NAME, "]:\tERROR target variable has more than two levels."))
            } else if (length(target_row) == 1) {
                results <- coef_table[target_row,]

                df_return[["p_value"]] <- results[[
                    names(results)[grepl(
                        pattern="Pr(",
                        x=names(results),
                        fixed=T
                    )]
                ]]
                df_return[["estimate"]] <- results[["Estimate"]]
                df_return[["std_error"]] <- results[["Std. Error"]]
                df_return[["test_statistic"]] <- results[[param[["test_statistic"]]]]
                df_return[["degrees_freedom"]] <- df.residual(model) # summary(model)$df[2]
                # Add in cooks distance for outlier filtering
                cooks_d <- tryCatch({
                    cooks_d <- cooks.distance(model)
                    cooks_d
                }, error = function(e) {
                    #message("Error: when calculating cooks: ", e$message, "\n")
                    NULL  # Return NULL or any other appropriate value
                })
                if (!is.null(cooks_d)) {
                    df_return[["cooks_distance_max"]] <- max(cooks_d)
                    # Cooks distance follows an F-distribution... use that to calcualte pvalues
                    # as in statsmodels
                    cooks_d_pvalues <- 1 - pf(
                        cooks_d, 
                        length(coef_table[, 1]), 
                        df.residual(model)
                    )
                    df_return[["cooks_distance_pvalue_min"]] <- min(cooks_d_pvalues)
                    df_return[["cooks_distance_pvalue_n_less_0pt05"]] <- sum(
                        cooks_d_pvalues <= 0.05
                    )
                    df_return[["cooks_distance_pvalue_n_less_0pt01"]] <- sum(
                        cooks_d_pvalues <= 0.01
                    )
                }
            
                df_return[["dispersion"]] <- sum(residuals(model, type = "pearson")^2) / df.residual(model)
                
                # Add overdispersion check
                # dispersion_test <- tryCatch({
                #     sim_fmp <- DHARMa::simulateResiduals(model, refit = TRUE)
                #     over_d_res <- DHARMa::testOverdispersion(sim_fmp)
                #     print(over_d_res)
                #     stop()
                # }, error = function(e) {
                #     message("Error: when performing overdispersion test: ", e$message, "\n")
                #     NULL  # Return NULL or any other appropriate value
                # })
                # if (is.null(dispersion_test)) {
                #     df_return[["dispersion_pvalue"]] <- NA
                # } else {
                #     df_return[["dispersion_pvalue"]] <- dispersion_test
                # }
                # Add in multicollinearity check
                tmp_vif <- tryCatch({
                    car::vif(model)
                }, error = function(e) {
                    #message("Error: when calculating VIF: ", e$message, "\n")
                    NULL  # Return NULL or any other appropriate value
                })
                if (!is.null(tmp_vif)) {
                    df_return[["vif_target_variable"]] <- tmp_vif[[param[["target_variable"]]]]
                    df_return[["vif_max"]] <- max(unlist(tmp_vif))
                }

                if (!is.na(df_return[["p_value"]])) {
                    if (param[["n_permutations"]] > 0) {
                        tmp_list <- run_permutations(
                            df_reg, 
                            param, 
                            df_return, 
                            n_permutations = param[["n_permutations"]],
                            early_stop_threshold = 0 # 0 means no early stopping
                        )
                        df_return[["p_value_empirical"]] <- tmp_list[["p_value_empirical"]]
                        df_return[["p_value_empirical_distribution"]] <- tmp_list[["p_value_empirical_distribution"]]
                        df_return[["n_permutations"]] <- tmp_list[["n_permutations"]]
                    }
                }
            }
        }

        return(df_return)
    } # end foreach loop

    # stop the cluster
    if (param[["parallel_type"]] != "domc") {
        parallel::stopCluster(cl)
    }
    time_end <- Sys.time()
    if (param[["verbose"]]) {
        message(paste0(
            "\nForeach loop execution time", " [", SCRIPT_NAME, "]:\t",
            difftime(time_end, time_init, units = c("hours")),
            " hours."
        ))
    }

    return(df_final)
}


main <- function() {
    # run analysis
    run_time <- system.time(df_results <- command_line_interface())
    message(paste0(
        "Analysis execution time", " [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]]/3600, # proc.time sec to hours
        " hours."
    ))
}


dev <- function() {
    # for profiling
    suppressPackageStartupMessages(library(profvis))

    param <- list()
    param[["anndata_file"]] <- "../data/demo_data/demo_data.h5ad"
    param[["output_file"]] <- "test"
    # param[["chunk_i"]] <- 1
    # param[["chunk_n"]] <- 2
    param[["formula"]] <- "cbind(counts_unmodified,counts_modified) ~ glucose_condition + library_size"
    param[["target_variable"]] <- "glucose_condition"
    param[["continuous_covariates"]] <- "library_size"
    param[["factor_covariates"]] <- "glucose_condition"
    param[["threads"]] <- 4
    param[["verbose"]] <- TRUE
    param[["parallel_type"]] <- "domc"

    base <- paste0("profile-", gsub(" ", "_", Sys.time()))
    # run the analysis to understand memory usage
    prof <- profvis::profvis(
        run_analysis(
            param, 
            fill_y_na_value = 0, 
            inverse_normalize_covariates = FALSE
        )
    )
    saveRDS(prof, paste0(base, ".Rds.gz"), compress = TRUE)
    #prof <- readRDS(paste0(base, ".Rds.gz"))
    #print(prof)

    param[["parallel_type"]] <- "snow"
    # only run if snow
    prof <- snow::snow.time(   
        run_analysis(
            param, 
            fill_y_na_value = 0, 
            inverse_normalize_covariates = FALSE
        )
    )
    pdf(file = paste0(base, ".pdf"), height = 5, width = 6)
        print(plot(prof))
    dev.off()
}

main()
#dev()