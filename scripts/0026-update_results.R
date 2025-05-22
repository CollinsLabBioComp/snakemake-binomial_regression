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
    
    optparse::make_option(c("--shrink_distribution"),
        type = "character",
        default = "normal",
        help = paste0(
          "Distribution to use for shrinkage.",
          " [default %default]"
        )
    ),
  
    optparse::make_option(c("--output_file"),
        type = "character",
        default = "results_corrected.tsv.gz",
        help = paste0(
          "Outfile.",
          " [default %default]"
        )
    ),
  
    optparse::make_option(c("--verbose"),
        action = "store_true",
        default = TRUE,
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
      "Corrects differential expression results results using Benjamini-Hochberg procedure."
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
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(qvalue))
################################################################################

################################ Functions #####################################
plot_hits <- function(
    df,
    facet_col,
    p_val_column,
    fdr_threshold = 0.05,
    base_theme = 20
) {
    plt_df <- df
    plt_df$facet <- as.factor(plt_df[[facet_col]])
    plt_df$fdr <- plt_df[[p_val_column]]
    plt_df <- plt_df %>%
        dplyr::group_by(facet) %>%
        dplyr::mutate(n_hits = sum(fdr < fdr_threshold)) %>%
        dplyr::select(facet, n_hits) %>%
        unique(.)
  
    plot <- ggplot2::ggplot(plt_df, ggplot2::aes(
          x=facet,
          y=n_hits
      )) +
      ggplot2::geom_col() +
      ggplot2::theme_bw(base_size=base_theme) +
      ggplot2::labs(
          x="Facet",
          y="Number of genes (FDR < 0.05)"
      ) +
      ggplot2::scale_y_continuous(trans="log10") +
      ggplot2::theme(legend.position="none")
      return(plot)
}
################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
output_file <- arguments$options$output_file
theme_size <- 60

# Get arguments
if (verbose) {
    print("Reading in the data...")
}

# Correct p-values
rez <- read.csv(
      arguments$options$results,
      sep = "\t",
      header = T
  ) %>%
  dplyr::mutate(qvalue_bh = p.adjust(p_value, method='BH')) %>%
  as.data.frame(.)

# Get shrunken estimate
ash_rez <- ashr::ash(
  rez$estimate,
  rez$std_error,
  method = "shrink",
  mixcompdist = "normal"
)
rez$estimate_shrunk <- ash_rez$result$PosteriorMean

## Save result
if (verbose) {
    print("Writing DE results...")
}
gz_file <- gzfile(output_file, "w", compression = 9)
write.table(
    x = rez,
    file = gz_file,
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
)
close(gz_file)

if (verbose) {
    print("Done.")
}

## Plot results
if (verbose) {
    print("Plotting hits per cell type...")
}

pdf(
    file = paste0(
        gsub(pattern='.tsv.gz', replacement='', x=output_file, fixed=T),
        "-hits.pdf"
    ),
    height = cm(9),
    width = cm(9)
)
print(plot_hits(
    rez,
    "model_id",
    "qvalue_bh",
    fdr_threshold = 0.05
))
dev.off()

if (verbose) {
  print("Done.")
}

################################################################################
