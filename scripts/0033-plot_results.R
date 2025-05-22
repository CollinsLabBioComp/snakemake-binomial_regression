#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
    optparse::make_option(c("--results"),
        type = "character",
        default = "results.tsv.gz",
        help = paste0(
          "Dataframe of results. Needs column `p_value`",
          " [default %default]"
        )
    ),
    
    optparse::make_option(c("--grouping_col"),
        type = "character",
        default = "model_id",
        help = paste0(
          "Grouping column to plot over.",
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
    
    optparse::make_option(c("--output_base"),
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
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
################################################################################

################################ Functions #####################################
plot_volcano_plot <- function(
    df,
    fc_col,
    p_val_col,
    facet_var
) {
    df$neg_log10 <- -log10(df[[p_val_col]])
    plot <- ggplot2::ggplot(df, ggplot2::aes_string(
          x = fc_col,
          y = "neg_log10",
          color = "significant",
          alpha = "significant"
      )) +
      ggplot2::geom_point(size = 0.5) +
      ggplot2::labs(
          x = "Estimate",
          y = "-log10(p-value)",
          color = bquote(FDR<0.05)
      ) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = c("#CF9400", "Black")) +
      ggplot2::scale_alpha_manual(values = c(0.75, 0.25)) +
      ggplot2::facet_wrap(as.formula(paste("~", facet_var)), nrow = 1)
  return(plot)
}


plot_ma_plot <- function(
    df,
    mean_expr_col,
    fc_col,
    facet_var
) {
    plot <- ggplot2::ggplot(df, ggplot2::aes_string(
          x=mean_expr_col,
          y=fc_col,
          color="significant",
          alpha = "significant"
      )) +
      ggplot2::geom_point(size = 0.5) +
      ggplot2::scale_x_continuous(
          trans = "log10",
          labels = scales::comma_format()
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
          x="Mean Expression",
          y="Estimate",
          color=bquote(FDR<0.05)
      ) +
      ggplot2::scale_color_manual(values=c("#CF9400", "Black")) +
      ggplot2::scale_alpha_manual(values = c(0.75, 0.25)) +
      ggplot2::facet_wrap(as.formula(paste("~", facet_var)), ncol = 1)
  return(plot)
}

plot_expression_filter <- function(
    df,
    x_column,
    y_column,
    color_column,
    facet_var,
    title_string = "",
    x_log = FALSE
) {
    plot <- ggplot2::ggplot(df, ggplot2::aes_string(
          x=x_column,
          y=y_column,
          color= color_column
      )) +
      ggplot2::geom_point(size = 0.5) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(as.formula(paste("~", facet_var)), ncol = 1)
  if (x_log) {
    plot <- plot + ggplot2::scale_x_log10()
  }
  if (title_string != "") {
    plot <- plot + ggplot2::labs(title = title_string)
  }
  return(plot)
}


################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_base <- arguments$options$output_base

# Read in results
rez <- read.csv(
        arguments$options$results,
        sep ="\t",
        header=T
    ) %>%
    dplyr::filter(!is.na(qvalue_bh)) %>%
    dplyr::mutate(
        significant = factor(
            qvalue_bh < arguments$options$significance_threshold,
            levels = c(TRUE, FALSE),
            labels = c("True", "False")
        ),
        p_value_neglog10 = -log10(p_value),
        estimate_big = abs(estimate) > 3,
        estimate_log = log10(abs(estimate)),
        cooks_distance_pvalue_min_neglog10 = -log10(cooks_distance_pvalue_min)
    ) %>%
    as.data.frame(.)

if (verbose) {
    cat(sprintf("Plotting results for %s...\n", arguments$options$input_file))
}

png(
    file = sprintf("%s-volcano.png", output_base),
    height = 10,
    width = 12,
    units = "cm",
    res = 320
)
print(plot_volcano_plot(
    rez,
    'estimate',
    'p_value',
    arguments$options$grouping_col
))
dev.off()

png(
    file = sprintf("%s-ma-count_mean.png", output_base),
    height = 10,
    width = 12,
    units = "cm",
    res = 320
)
print(plot_ma_plot(
    rez,
    'count_mean',
    'estimate',
    arguments$options$grouping_col
))
dev.off()

png(
    file = sprintf("%s-ma-count_normalized_mean.png", output_base),
    height = 10,
    width = 12,
    units = "cm",
    res = 320
)
print(plot_ma_plot(
    rez,
    'count_normalized_mean',
    'estimate',
    arguments$options$grouping_col
))
dev.off()


pdf(
    file = sprintf("%s-expression_filter.pdf", output_base),
    height = 10,
    width = 12
)
print(plot_expression_filter(
    rez,
    'count_normalized_mean',
    'p_value_neglog10',
    'significant',
    arguments$options$grouping_col,
    title_string = "Normalized count mean vs. -log10(p-value)",
    x_log = TRUE
))
print(plot_expression_filter(
    rez,
    'count_mean',
    'p_value_neglog10',
    'significant',
    arguments$options$grouping_col,
    title_string = "Count mean vs. -log10(p-value)",
    x_log = TRUE
))
print(plot_expression_filter(
    rez,
    'cooks_distance_pvalue_min_neglog10',
    'estimate',
    'count_normalized_mean',
    arguments$options$grouping_col,
    title_string = "Cook's distance p-value min vs. effect size estimate"
))
print(plot_expression_filter(
    rez,
    'count_normalized_mean',
    'count_normalized_pct_grtr_eq_mean',
    'estimate',
    arguments$options$grouping_col,
    title_string = "Normalized count mean vs. percent greater than mean,\ncolored by effect size estimate",
    x_log = TRUE
))
if (nrow(subset(rez, estimate_big)) > 0) {
    print(plot_expression_filter(
    subset(rez, estimate_big),
        'count_normalized_mean', 
        'count_normalized_pct_grtr_eq_mean',
        'estimate',
        arguments$options$grouping_col,
        title_string = "Normalized count mean vs. percent greater than mean,\ncolored by effect size estimate (|estimate| > 3)",
        x_log = TRUE
    ))
}
if (nrow(subset(rez, qvalue_bh < 0.05)) > 0) {
    print(plot_expression_filter(
        subset(rez, qvalue_bh < 0.05),
        'count_normalized_mean',
        'count_normalized_pct_grtr_eq_mean',
        'estimate',
        arguments$options$grouping_col,
        title_string = "Normalized count mean vs. percent greater than mean,\ncolored by effect size estimate (fdr < 0.05)",
        x_log = TRUE
    ))
}
if (!all(is.na(rez$cooks_distance_pvalue_min))) {
    print(plot_expression_filter(
        rez,
        'count_normalized_mean',
        'count_normalized_pct_grtr_eq_mean',
        'cooks_distance_pvalue_min',
        arguments$options$grouping_col,
        title_string = "Normalized count mean vs. percent greater than mean,\ncolored by Cook's distance p-value min",
        x_log = TRUE
    ))
}
dev.off()


if (verbose) {
  cat("Done.\n")
}

################################################################################
