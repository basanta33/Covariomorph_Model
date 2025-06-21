#!/user/bin/env Rscript
# R script to plot the split frequencies of 3 rate varying configurations
# Authors: Basanta Khakurel and Sebastian Höhna (Modified from Alessio Capobianco's script)
# date: 2025-06-10

library(convenience)
library(tools)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(pilot)
library(extrafont)

set_pilot_family(family = "Montserrat")

OUTPUT_DIR <- "Plots/Empirical_Partitioned/CPP-Plot"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ESS <- 200
BURNIN <- 0.25
expDiffSplits <- convenience:::expectedDiffSplits(ess = ESS)

models <- c(
  "Mk",
  "ACRV_4",
  "Covarion_4"
)

model_labels <- c(
  "Mk" = "Mk",
  "ACRV_4" = "ACRV (m=4)",
  "Covarion_4" = "Covariomorph (m=4)"
)

get_file_path <- function(model, dataset_name) {
  base_output_dir <- "output_Empirical_Partitioned"

  if (model == "Mk") {
    return(file.path(base_output_dir, "output_Mk", paste0(dataset_name, "_Mk_0_Cats.trees")))
  } else if (grepl("^ACRV", model)) {
    k_value <- gsub("ACRV_", "", model)
    return(file.path(
      base_output_dir, paste0("output_ACRV_", k_value),
      paste0(dataset_name, "_ACRV_", k_value, "_Cats.trees")
    ))
  } else if (grepl("^Covarion", model)) {
    k_value <- gsub("Covarion_", "", model)
    return(file.path(
      base_output_dir, paste0("output_Covarion_", k_value, "Cats"),
      paste0(dataset_name, "_Covarion_Part_", k_value, "_Cats.trees")
    ))
  }
}

plot_splits <- function(file1, file2, label1, label2, burnin, expDiffSplits) {
  message(paste("Comparing:", label1, "vs", label2))
  comparison_data <- NULL

  if (!file.exists(file1) || !file.exists(file2)) {
    message(paste("Warning: One or both files do not exist:", file1, file2))
    return(NULL)
  }

  tryCatch(
    {
      comparison_data <- checkConvergence(list_files = c(file1, file2), control = makeControl(burnin = burnin))
    },
    error = function(e) {
      message(paste("Error in checkConvergence:", e$message))
      return(NULL)
    }
  )

  if (is.null(comparison_data) || is.null(comparison_data$tree_parameters) || is.null(comparison_data$tree_parameters$freq_per_run)) {
    message("No valid comparison data found")
    return(NULL)
  }

  splits <- as.data.frame(comparison_data$tree_parameters$freq_per_run)
  col_label1 <- make.names(label1)
  col_label2 <- make.names(label2)
  colnames(splits)[1:2] <- c(col_label1, col_label2)
  na.rowsIndex <- which(rowSums(is.na(splits[, c(col_label1, col_label2)])) == 0)

  if (length(na.rowsIndex) == 0) {
    message("No valid splits found after filtering NA values")
    return(NULL)
  }
  splits <- splits[na.rowsIndex, , drop = FALSE]

  converged.splits.idx.rownames <- c()
  diverged.splits.idx.rownames <- c()

  for (split_idx_rowname in rownames(splits)) {
    run1_freq <- splits[split_idx_rowname, col_label1]
    run2_freq <- splits[split_idx_rowname, col_label2]

    average <- (run1_freq + run2_freq) / 2
    diff <- abs(run1_freq - run2_freq)

    target.index <- which.min(abs(expDiffSplits[1, ] - average))

    if (diff > expDiffSplits[2, target.index]) {
      diverged.splits.idx.rownames <- c(diverged.splits.idx.rownames, split_idx_rowname)
    } else {
      converged.splits.idx.rownames <- c(converged.splits.idx.rownames, split_idx_rowname)
    }
  }

  bivariate_plot <- ggplot(data = splits, aes_string(x = col_label1, y = col_label2)) +
    geom_point(data = splits[converged.splits.idx.rownames, ], size = 2, alpha = 0.7, colour = "black", shape = 16) + # Converged = black dots
    geom_point(data = splits[diverged.splits.idx.rownames, ], size = 2, alpha = 0.7, colour = "blue", shape = 17) + # Diverged = blue triangles
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = label1, y = label2) +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    theme_pilot(axis_title_size = 14, axis_text_size = 10, title_size = 16) +
    theme(plot.margin = margin(5, 5, 5, 5))

  return(list(bivariate_plot = bivariate_plot))
}

# function to label the rows
create_y_axis_label <- function(text_label) {
  ggplot() +
    annotate(
      geom = "text", x = 0.5, y = 0.5,
      label = text_label,
      angle = 90,
      size = 7,
      family = "Montserrat",
      fontface = "bold"
    ) +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = -1, b = 0, l = 0, unit = "pt"))
}

# selected datasets
dataset_info <- list(
  list(name = "Dasyatoidea_MarramaEtAl2023_All", label = "Rays"),
  list(name = "Shirai_1996a", label = "Sharks")
)


all_plot_rows <- list()

# loop along the datasets to obtain comparison plots
for (info in dataset_info) {
  dataset_name <- info$name
  y_axis_label <- info$label

  message(paste("\n\n--- Processing dataset:", dataset_name, "---"))

  comparisons <- list()
  for (i in 1:(length(models) - 1)) {
    for (j in (i + 1):length(models)) {
      comparisons[[length(comparisons) + 1]] <- list(
        file1 = get_file_path(models[i], dataset_name),
        file2 = get_file_path(models[j], dataset_name),
        label1 = model_labels[models[i]],
        label2 = model_labels[models[j]]
      )
    }
  }

  plot_results <- lapply(comparisons, function(comp) {
    plot_splits(
      file1 = comp$file1, file2 = comp$file2,
      label1 = comp$label1, label2 = comp$label2,
      burnin = BURNIN,
      expDiffSplits = expDiffSplits
    )
  })

  bivariate_plots_list <- Filter(Negate(is.null), lapply(plot_results, `[[`, "bivariate_plot"))

  if (length(bivariate_plots_list) > 0) {
    y_label_plot <- create_y_axis_label(y_axis_label)
    plots_row <- wrap_plots(bivariate_plots_list, nrow = 1)
    full_row_with_label <- y_label_plot + plots_row + plot_layout(widths = c(1, 20))
    all_plot_rows <- append(all_plot_rows, list(full_row_with_label))
  }
}

# create a final combined plot
if (length(all_plot_rows) > 0) {
  message("\n--- Combining all dataset rows into a final plot ---")

  final_combined_plot <- wrap_plots(all_plot_rows, ncol = 1)

  output_filename <- file.path(OUTPUT_DIR, "clade_posterior_probabilities_combined.pdf")
  ggsave(
    filename = output_filename,
    plot = final_combined_plot,
    width = 12,
    height = 6,
    dpi = 450,
    device = cairo_pdf,
    limitsize = FALSE
  )
  message(paste("\nScript completed successfully. Final plot saved to:", output_filename))
} else {
  message("\nNo valid plots were generated to create a final plkot.")
}
