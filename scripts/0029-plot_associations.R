#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

#' Perform regression
#'
#' @param df data.frame.
#'     Data frame with regression information
#'     
#' @param p List.
#'     List of parameters
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

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
    optparse::make_option(c("--anndata_file"),
        type = "character",
        help = paste0(
          "Anndata file with the following layers: counts_modified, counts_unmodified",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--results"),
        type = "character",
        default = "results.tsv.gz",
        help = paste0(
          "Dataframe of results. Needs column `p_value`",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--model_id"),
        type = "character",
        default = "model1",
        help = paste0(
          "Model to select.",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--significance_column"),
        type = "character",
        default = "pvalue",
        help = paste0(
          "Significance column to threshold.",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--significance_threshold"),
        type = "double",
        default = 0.05,
        help = paste0(
          "Threshold for significance.",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--sort_column"),
        type = "character",
        default = "pvalue",
        help = paste0(
          "Sorting column.",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--plot_top_n_associations"),
        type = "integer",
        default = 10,
        help = paste0(
          "Top associations to plot.",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--plot_bottom_n_associations"),
        type = "integer",
        default = 10,
        help = paste0(
          "Bottom associations to plot.",
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

    
    optparse::make_option(c("--out_base"),
        type = "character",
        default = "results_corrected",
        help = paste0(
          "Output base.",
          " [default %default]"
        )
    ),

    optparse::make_option(c("--verbose"),
        action = "store_true",
        default = FALSE,
        help = paste0(
          "Verbose",
          " [default %default]"
        )
    )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Plots results from differential gene expression."
  )
)

# a hack to fix a bug in optparse that won't let you use positional args
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
################################################################################

######################## Required Packages #####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(anndataR))
################################################################################

################################ Functions #####################################
read_adata <- function(adat_file, filter) {
    adata <- anndataR::read_h5ad(adat_file, to = "InMemoryAnnData")
    
    # Deal with filter if exists
    if (!is.na(filter) & (filter != "")) {
        obs_ret <- adata$obs %>%
            dplyr::filter(!! rlang::parse_expr(filter)) %>%
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
  
  return(adata)
}


format_dataframe <- function(
    h5,
    feature_id,
    factor_covs,
    cont_covs,
    params = list(),
    lib_size_normalize = F
) {
    filt_row <- row.names(h5$var) == feature_id
    
    # Data frame of variables for the regression
    # Assumes adata$obs + adata$layers contains all of the variables
    df <- h5$obs
    df$sample_id <- rownames(df)
    for (i in names(h5$layers)) {
        df[[i]] <- h5$layers[[i]][, filt_row]
        
        if (lib_size_normalize) {
            df[[i]] <- df[[i]] / df[['library_size']]
        }
    }
    
    #print(df)
    # Assumes the following layers exist: `counts_unmodified` and `counts_modified`
    df <- df %>%
        dplyr::filter(!is.na(counts_unmodified)) %>%
        dplyr::filter(!is.na(counts_modified)) %>%
        as.data.frame(.)

    # Also process covariates for ease of plotting
    for (f in factor_covs) {
        if (grepl(pattern = '::', f, fixed=T)) {
            fact_var <- strsplit(f, split='::', fixed=T)[[1]][1]
            fact_ref <- strsplit(f, split='::', fixed=T)[[1]][2]
            df[[fact_var]] <- factor(
                df[[fact_var]],
                levels = c(
                    fact_ref,
                    setdiff(unique(df[[fact_var]]), fact_ref)
                )
            )
        } else {
            df[[f]] <- as.factor(df[[f]])
        }
    }
    for (i in cont_covs) {
        # don't scale or normalize for plotting
        df[[i]] <- as.numeric(df[[i]])
    }
    
    df$residuals <- NA
    if (length(params) > 0) {
        model <- tryCatch({
            regress(df, params)
        }, error = function(e) {
            NULL  # If the model fails, return NULL
        })
        if (!is.null(model)) {
            df$residuals <- residuals(model)
        }
    }
   
    # format for plotting
    df <- df %>%
        dplyr::mutate(
            count_ratio = counts_modified / (counts_unmodified + counts_modified),
            total_counts = counts_unmodified + counts_modified
        ) %>%
        tidyr::pivot_longer(
            cols = c('counts_unmodified', 'counts_modified'),
            names_to = 'mod_type',
            values_to = 'counts'
        ) %>%
        dplyr::mutate(
            mod_type = factor(
                mod_type,
                levels=c('counts_unmodified', 'counts_modified'),
                labels=c('Unmodified', 'Modified')
            )
        ) %>%
        as.data.frame(.)

    return(df)
}


plot_association__continuous <- function(
    df,
    x_col,
    base_theme,
    y_label
) {
    y_label <- paste(y_label, " ratio modified / total")
    df$x <- df[[x_col]]
    df <- df %>%
      dplyr::group_by(sample_id, x, count_ratio) %>%
      unique(.)
    
    plot <- ggplot2::ggplot(df, ggplot2::aes(
          x = x,
          y = count_ratio,
          size = total_counts
      )) +
      ggplot2::geom_point(size = 5, alpha = 0.75) +
      ggplot2::theme_bw(base_size = base_theme) +
      ggplot2::labs(
          x=x_col,
          y=y_label
      ) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1.5) +  # Add linear regression
      ggplot2::theme(
          title=element_blank(),
          legend.position='bottom',
          legend.direction='horizontal',
          legend.box='horizontal',
          legend.box.margin=margin(t=-30,r=15,b=0,l=0),
          legend.margin=margin(t=-20,r=0,b=0,l=0),
          legend.title=element_text(
            margin=margin(t=0,r=0,b=10,l=0),
            size=35
          ),
          legend.text=element_text(size=30),
          legend.spacing.y=unit(0, 'cm'),
          plot.background=element_blank(),
          axis.title.y=element_text(size=40),
          axis.text.y=element_text(size=35),
          axis.text.x=element_text(size=35),
          axis.title.x=element_text(size=40),
          panel.grid.minor = element_blank()
      )
  return(plot)
}

plot_association__discrete <- function(
    df,
    facet_col,
    base_theme,
    y_label
) {
    df$facet <- df[[facet_col]]
    plot <- ggplot2::ggplot(df, ggplot2::aes(
          x = mod_type,
          y = counts,
          color = mod_type
      )) +
      ggplot2::geom_boxplot(width = 0.25, linewidth = 2, alpha = 0.9) +
      ggplot2::geom_point(size = 5, alpha = 0.75) +
      ggplot2::geom_line(
          mapping = ggplot2::aes(group = sample_id),
          color = 'black',
          size = 2.5,
          alpha = 0.6
      ) +
      ggplot2::theme_bw(base_size = base_theme) +
      ggplot2::labs(
          x='Modification type',
          y=y_label,
          color='Modification type'
      ) +
      ggplot2::scale_color_brewer(palette = 'Dark2') +
      ggplot2::scale_fill_brewer(palette = 'Dark2') +
      ggplot2::facet_wrap(~ facet, nrow = 1) +
      ggplot2::theme(
          title=element_blank(),
          legend.position='none',
          legend.direction='horizontal',
          legend.box='horizontal',
          legend.box.margin=margin(t=-30,r=15,b=0,l=0),
          legend.margin=margin(t=-20,r=0,b=0,l=0),
          legend.title=element_text(
              margin=margin(t=0,r=0,b=10,l=0),
              size=35
          ),
          legend.text=element_text(size=30),
          legend.spacing.y=unit(0, 'cm'),
          plot.background=element_blank(),
          axis.title.y=element_text(size=40),
          axis.text.y=element_text(size=35),
          axis.text.x=element_text(size=35),
          axis.title.x=element_text(size=40),
          panel.grid.minor = element_blank()
      ) +
      ggplot2::guides(
        color=ggplot2::guide_legend(
            title.position="top",
            title.hjust = 0.5,
            override.aes=list(
                shape = 15,
                alpha = 0.9,
                size = 15
            ),
            keywidth = 0
        )
      )
    return(plot)
}


plot_association__discrete_ratio <- function(
    df,
    facet_col,
    base_theme,
    y_label
) {
    df$facet <- df[[facet_col]]
    plot <- ggplot2::ggplot(df, ggplot2::aes(
          x = facet,
          y = count_ratio,
          fill = facet,
          size = total_counts
      )) +
      ggplot2::geom_boxplot(width = 0.25, linewidth = 2, alpha = 0.9) +
      ggplot2::geom_jitter(alpha = 0.75, width = 0.1) +
      #ggplot2::geom_point(size = 5, alpha = 0.75) +
      ggplot2::theme_bw(base_size = base_theme) +
      ggplot2::labs(
          x='Phenotype',
          y="Ratio modified / total",
          size=paste0("Total ", y_label)
      ) +
      ggplot2::scale_color_brewer(palette = 'Dark2', guide="none") +
      ggplot2::scale_fill_brewer(palette = 'Dark2', guide="none") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme(
          title=element_blank(),
        #   legend.position='none',
        #   legend.direction='horizontal',
        #   legend.box='horizontal',
        #   legend.box.margin=margin(t=-30,r=15,b=0,l=0),
        #   legend.margin=margin(t=-20,r=0,b=0,l=0),
        #   legend.title=element_text(
        #       margin=margin(t=0,r=0,b=10,l=0),
        #       size=35
        #   ),
        #   legend.text=element_text(size=30),
        #   legend.spacing.y=unit(0, 'cm'),
          plot.background=element_blank(),
          axis.title.y=element_text(size=40),
          axis.text.y=element_text(size=35),
          axis.text.x=element_text(size=35),
          axis.title.x=element_text(size=40),
          panel.grid.minor = element_blank()
      ) +
      ggplot2::guides(
        color=ggplot2::guide_legend(
            title.position="top",
            title.hjust = 0.5,
            override.aes=list(
                shape = 15,
                alpha = 0.9,
                size = 15
            ),
            keywidth = 0
        )
      )
    return(plot)
}

plot_association__discrete_residuals <- function(
    df,
    facet_col,
    base_theme,
    y_label
) {
    df$facet <- df[[facet_col]]
    plot <- ggplot2::ggplot(df, ggplot2::aes(
          x = facet,
          y = residuals,
          fill = facet,
          size = total_counts
      )) +
      ggplot2::geom_boxplot(width = 0.25, linewidth = 2, alpha = 0.9) +
      ggplot2::geom_jitter(alpha = 0.75, width = 0.1) +
      #ggplot2::geom_point(size = 5, alpha = 0.75) +
      ggplot2::theme_bw(base_size = base_theme) +
      ggplot2::labs(
          x='Phenotype',
          y="Residuals",
          size=paste0("Total ", y_label)
      ) +
      ggplot2::scale_color_brewer(palette = 'Dark2', guide="none") +
      ggplot2::scale_fill_brewer(palette = 'Dark2', guide="none") +
      ggplot2::theme(
          title=element_blank(),
        #   legend.position='none',
        #   legend.direction='horizontal',
        #   legend.box='horizontal',
        #   legend.box.margin=margin(t=-30,r=15,b=0,l=0),
        #   legend.margin=margin(t=-20,r=0,b=0,l=0),
        #   legend.title=element_text(
        #       margin=margin(t=0,r=0,b=10,l=0),
        #       size=35
        #   ),
        #   legend.text=element_text(size=30),
        #   legend.spacing.y=unit(0, 'cm'),
          plot.background=element_blank(),
          axis.title.y=element_text(size=40),
          axis.text.y=element_text(size=35),
          axis.text.x=element_text(size=35),
          axis.title.x=element_text(size=40),
          panel.grid.minor = element_blank()
      ) +
      ggplot2::guides(
        color=ggplot2::guide_legend(
            title.position="top",
            title.hjust = 0.5,
            override.aes=list(
                shape = 15,
                alpha = 0.9,
                size = 15
            ),
            keywidth = 0
        )
      )
    return(plot)
}

plot_association <- function(
    df,
    target_col,
    title,
    base_theme,
    y_label,
    ratio = FALSE, 
    residuals = FALSE
) {
    if (is.numeric(df[[target_col]])) {
        plt <- plot_association__continuous(
            df,
            target_col,
            base_theme, 
            y_label 
        )
    } else {
        if (ratio) {
            plt <- plot_association__discrete_ratio(
                df,
                target_col,
                base_theme,
                y_label
            )
        } else if (residuals) {
            plt <- plot_association__discrete_residuals(
                df,
                target_col,
                base_theme,
                y_label
            )
        } else {
            plt <- plot_association__discrete(
                df,
                target_col,
                base_theme,
                y_label
            )
        }
    }
  
    plt <- plt +
      ggplot2::labs(title = title) +
      ggplot2::theme(title=element_text(size=40))
    return(plt)
}
################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_base <- arguments$options$out_base
theme_size <- 60

# Read in anndata, format covs
adata <- read_adata(
    arguments$options$anndata_file,
    arguments$options$filter
)
targ_var <- arguments$options$target_variable
cont_covs <- strsplit(arguments$options$continuous_covariates, split = ",")[[1]]
factor_covs <- strsplit(arguments$options$factor_covariates, split = ",")[[1]]

# Read in results
rez <- read.csv(
    arguments$options$results,
    sep ="\t",
    header=T
)
rez$thresh_col <- rez[[arguments$options$significance_column]]
rez$sort_col <- rez[[arguments$options$sort_column]]

rez <- rez %>%
    dplyr::filter(model_id == arguments$options$model_id) %>%
    dplyr::filter(!is.na(thresh_col))

if (verbose) {
    cat(sprintf("Plotting results for %s...\n", arguments$options$results))
}

# Top hits first
top_rez <- rez %>%
  dplyr::arrange(sort_col) %>%
  dplyr::filter(thresh_col <= arguments$options$significance_threshold) %>%
  as.data.frame(.)

if (nrow(top_rez) != 0) {
    top_ids <- top_rez$position_id[
        0:min(nrow(top_rez), arguments$options$plot_top_n_associations)
    ]
    
    pdf(
        file = sprintf("%s/top_associations.pdf", output_base),
        height = cm(9),
        width = cm(9)
    )
    for (feat in top_ids) {
        tmp <- subset(top_rez, position_id == feat)
        #print(tmp)

        # Process the formula and drop the target term
        original_formula <- as.formula(tmp$formula)
        left_hand_side <- as.character(original_formula)[2]
        formula_terms <- attr(terms(original_formula), "term.labels")
        remaining_terms <- setdiff(formula_terms, tmp$target_variable)
        new_formula <- if (length(remaining_terms) == 0) {
            as.formula(paste0(left_hand_side, "~ 1"))
        } else {
            as.formula(paste0(left_hand_side, "~ ", paste(remaining_terms, collapse = " + ")))
        }
        params <- list(
            formula = deparse(new_formula),
            distribution = arguments$options$distribution
        )

        # First raw data
        feat_title <- sprintf('%s - raw data', feat)
        df_plt <- format_dataframe(adata, feat, factor_covs, cont_covs, params)
        print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts", ratio = F))
        print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts", ratio = T))
        # Now normalize data
        feat_title <- sprintf('%s - normalized data', feat)
        df_plt <- format_dataframe(adata, feat, factor_covs, cont_covs, params, lib_size_normalize = T)
        print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts (normalized)", ratio = F))
        print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts (normalized)", ratio = T))
        # Now residuals
        feat_title <- sprintf('%s - model residuals', feat)
        print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Residuals", ratio = F, residuals = T))

    }
    dev.off()
}

# Now bottom hits
null_rez <- rez %>%
    dplyr::arrange(desc(sort_col)) %>%
    as.data.frame(.)

null_ids <- null_rez$position_id[
    0:min(nrow(null_rez), arguments$options$plot_bottom_n_associations)
]

pdf(
    file = sprintf("%s/null_associations.pdf", output_base),
    height = cm(9),
    width = cm(9)
)
for (feat in null_ids) {
    tmp <- subset(null_rez, position_id == feat)
    #print(tmp)

    # Process the formula and drop the target term
    original_formula <- as.formula(tmp$formula)
    left_hand_side <- as.character(original_formula)[2]
    formula_terms <- attr(terms(original_formula), "term.labels")
    remaining_terms <- setdiff(formula_terms, tmp$target_variable)
    new_formula <- if (length(remaining_terms) == 0) {
        as.formula(paste0(left_hand_side, "~ 1"))
    } else {
        as.formula(paste0(left_hand_side, "~ ", paste(remaining_terms, collapse = " + ")))
    }
    params <- list(
        formula = deparse(new_formula),
        distribution = arguments$options$distribution
    )

    # First raw data
    feat_title <- sprintf('%s - raw data', feat)
    df_plt <- format_dataframe(adata, feat, factor_covs, cont_covs, params)
    print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts", ratio = F))
    print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts", ratio = T))
    # Now normalize data
    feat_title <- sprintf('%s - normalized data', feat)
    df_plt <- format_dataframe(adata, feat, factor_covs, cont_covs, params, lib_size_normalize = T)
    print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts (normalized)", ratio = F))
    print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Counts (normalized)", ratio = T))
    # Now residuals
    feat_title <- sprintf('%s - model residuals', feat)
    print(plot_association(df_plt, targ_var, feat_title, theme_size, y_label = "Residuals", ratio = F, residuals = T))

}
dev.off()

if (verbose) {
  cat("Done.\n")
}

################################################################################
