#!/usr/bin/Rscript

# Rscript to combine data simulated under covarion model
# the data is simulated under various rate categories which needs to be corrected.
# Only the state labels from the first Q matrix must be used.
# This script does the following:
#           - read a nexus data as a dataframe
#           - combine data simulated under covarion model the data is simulated
#             under various rate categories which needs to be corrected.
#             Only the state labels from the first Q matrix must be used.
#           - Build a cladistic matrix using Claddis
#           - Write the nexus matrix using Claddis such that the state labels appear correctly
# Authors: Basanta Khakurel and Sebastian Höhna
library(optparse, quietly = TRUE)
library(Claddis, quietly = TRUE)

option_list <- list(
  make_option(c("-d", "--dir"),
    type = "character", default = NULL,
    help = "Directory with the nexus files for conversion", metavar = "character"
  ),
  make_option(c("-n", "--n_sims"),
    type = "integer", default = 100,
    help = "Number of simulations", metavar = "number"
  ),
  make_option(c("-s", "--n_states"),
    type = "integer", default = 4,
    help = "Number of original states", metavar = "number"
  ),
  make_option(c("-c", "--n_cats"),
    type = "integer", default = 4,
    help = "Number of rate categories", metavar = "number"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dir)) {
  print_help(opt_parser)
  stop("Syntax error: missing directory argument", call. = FALSE)
}

DIR <- normalizePath(opt$dir)
n_sims <- opt$n_sims
n_states <- opt$n_states
n_cats <- opt$n_cats

# function to change the state labels from expanded form (as RevBayes Simulates) to shrinked form (as you see in morphological data matrices)
correct_states <- function(data, n_states, n_cats) {
  unique_states <- sort(unique(unlist(data)))
  expected_states <- length(unique_states) / n_cats

  if (n_states != expected_states) {
    stop("Mismatch in expected number of states. Check input files.", call. = FALSE)
  }

  mapping <- setNames(as.character(0:(n_states - 1))[rep(1:n_states, times = n_cats)], unique_states)

  return(as.data.frame(lapply(data, function(col) mapping[as.character(col)])))
}


for (i in 1:n_sims) {
  file_path <- file.path(DIR, paste0("morpho_", i, ".nex"))
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  morpho_data <- as.data.frame(read.nexus.data(file_path))
  corrected_data <- correct_states(morpho_data, n_states, n_cats)
  matrix <- as.matrix(t(corrected_data))

  cladistic_matrix <- build_cladistic_matrix(
    matrix,
    header = paste("Covarion simulated data --", n_states, "states", n_cats, "rate categories")
  )
  write_nexus_matrix(cladistic_matrix, file_path)

  # Remove assumptions block if present
  nexus_lines <- readLines(file_path)
  assumption_start <- grep("BEGIN ASSUMPTIONS;", nexus_lines, ignore.case = TRUE)
  if (length(assumption_start) > 0) {
    writeLines(nexus_lines[1:(assumption_start - 1)], file_path)
  }
}
q()
