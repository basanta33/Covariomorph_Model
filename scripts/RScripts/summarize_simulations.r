#!/usr/bin/env Rscript
# Rscript to obtain median of the simulation outputs
# Author: Basanta Khakurel

library(optparse, quietly = TRUE)
library(dplyr)

option_list <- list(
  make_option(c("-d", "--dir"),
    type = "character", default = NULL,
    help = "output directory", metavar = "character"
  ),
  make_option(c("-n", "--n_sims"),
    type = "integer", default = 10,
    help = "Number of simulations", metavar = "number"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

OUTPUT_DIR <- opt$dir
NUM_SIMULATIONS <- opt$n_sims

RESULTS_DIR <- paste0("Results_", OUTPUT_DIR)
dir.create(RESULTS_DIR, showWarnings = FALSE)

results_sims <- data.frame()

for (i in 1:NUM_SIMULATIONS) {
  log_file <- sprintf("%s/morpho_%d_Covarion.log", OUTPUT_DIR, i)

  if (file.exists(log_file)) {
    tmp <- read.delim(log_file)
    result <- data.frame(
      filename = i,
      sd_category = median(tmp$sd_category, na.rm = TRUE),
      global_switch_rate = median(tmp$global_switch_rate, na.rm = TRUE),
      tree_length = median(tmp$tree_length, na.rm = TRUE)
    )

    results_sims <- rbind(results_sims, result)
  } else {
    warning(sprintf("File %s does not exist", log_file))
  }
}

output_file <- paste0(RESULTS_DIR, "/morpho_sims_median.csv")
write.csv(results_sims, output_file, row.names = FALSE)

cat("Results saved to", output_file, "\n")
