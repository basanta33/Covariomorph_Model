# Rscript: summarizing covarion outputs
# Author: Basanta Khakurel

# load required libraries
library(optparse, quietly = TRUE)
library(dplyr)

# to obtain command line arguments
option_list <- list(
  make_option(c("-d", "--dir"),
    type = "character", default = NULL,
    help = "directory with the analysed data", metavar = "character"
  ),
  make_option(c("-f", "--file_name"),
    type = "character", default = NULL,
    help = "filename with the list of analysed data", metavar = "character"
  ),
  make_option(c("-m", "--model_name"),
    type = "character", default = "Cov_Zero",
    help = "model name (corresponding to the output dir)", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

LIST_DATA <- opt$file_name
OUTPUT_DIR <- opt$dir
MODEL <- opt$model_name

# create the results directory using the outpur directory
RESULTS_DIR <- paste0("results_", OUTPUT_DIR)
dir.create(RESULTS_DIR, showWarnings = FALSE)

# list of data files containg the output
data_list <- readLines(LIST_DATA)

# function to process the log files
process_file <- function(filename) {
  log_file <- sprintf("%s/%s_%s.log", OUTPUT_DIR, filename, MODEL)
  tmp <- read.delim(log_file)

  result <- data.frame(
    filename = filename,
    sd_category = median(tmp$sd_category, na.rm = TRUE),
    global_switch_rate = median(tmp$global_switch_rate, na.rm = TRUE),
    tree_length = median(tmp$tree_length, na.rm = TRUE)
  )

  # to accomodate for the different number of rates in different analyses
  max_rates <- min(12, sum(grepl("rate_scalars", names(tmp))))

  rate_columns <- sapply(1:max_rates, function(i) median(tmp[[paste0("rate_scalars.", i, ".")]], na.rm = TRUE))
  rate_names <- paste0("rate", 1:max_rates)

  result[rate_names] <- rate_columns
  return(result)
}

results <- do.call(rbind, lapply(data_list, process_file))

results <- results %>% mutate(filename = as.character(filename))

# the results are saved in the results directory with the name of the output directory
write.csv(results, paste0(RESULTS_DIR, "/", OUTPUT_DIR, ".csv"), row.names = FALSE)

q()
