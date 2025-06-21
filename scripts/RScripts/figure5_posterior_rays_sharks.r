#!/user/bin/env Rscript
# R script to plot the posteriors of sharks and rays datasets.
# Authors: Basanta Khakurel and Sebastian Höhna
# date: 2025-06-16

# Load necessary libraries
library(ggview, quietly = TRUE)
library(gridExtra)
library(grid)
library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggridges, quietly = TRUE)
library(pilot, quietly = TRUE)
library(extrafont)
library(reshape2)

set_pilot_family(family = "Montserrat")

# function to obtain plot for a single dataset
generate_plots <- function(filename, highlight_cat) {
  # Define output directories
  OUTPUT_DIRS <- c(
    "output_Covarion_2Cats", "output_Covarion_3Cats", "output_Covarion_4Cats",
    "output_Covarion_5Cats", "output_Covarion_6Cats", "output_Covarion_8Cats",
    "output_Covarion_10Cats", "output_Covarion_12Cats", "output_Covarion_14Cats",
    "output_Covarion_16Cats", "output_Covarion_18Cats"
  )
  op_mk <- c("output_Mk", "output_ACRV_4", "output_ACRV_8")

  # Initialize lists to store data
  sd_list <- list()
  switch_rate_list <- list()
  tree_length <- list()

  for (output_dir in OUTPUT_DIRS) {
    logfiles <- list.files(file.path("output_Empirical_Partitioned", output_dir), pattern = "\\.log$", full.names = TRUE)
    logfile <- logfiles[grepl(filename, logfiles) & !grepl("run", logfiles)]
    if (length(logfile) > 0) {
      data <- read.delim(logfile)
      sd_list[[output_dir]] <- data$sd_category
      switch_rate_list[[output_dir]] <- data$global_switch_rate
      tree_length[[output_dir]] <- data$tree_length
    }
  }

  for (output_dir in op_mk) {
    logfiles <- list.files(file.path("output_Empirical_Partitioned", output_dir), pattern = "\\.log$", full.names = TRUE)
    logfile <- logfiles[grepl(filename, logfiles) & !grepl("run", logfiles)]
    if (length(logfile) > 0) {
      data <- read.delim(logfile)
      sd_list[[output_dir]] <- NA
      switch_rate_list[[output_dir]] <- NA
      tree_length[[output_dir]] <- data$tree_length
    }
  }

  combined_data <- bind_rows(
    lapply(names(sd_list), function(dir) {
      data.frame(
        sd_category = sd_list[[dir]],
        global_switch_rate = switch_rate_list[[dir]],
        tree_length = tree_length[[dir]],
        output_dir = dir
      )
    })
  )

  combined_data$output_dir <- recode(combined_data$output_dir,
    "output_Covarion_2Cats" = "2", "output_Covarion_3Cats" = "3",
    "output_Covarion_4Cats" = "4", "output_Covarion_5Cats" = "5",
    "output_Covarion_6Cats" = "6", "output_Covarion_8Cats" = "8",
    "output_Covarion_10Cats" = "10", "output_Covarion_12Cats" = "12",
    "output_Covarion_14Cats" = "14", "output_Covarion_16Cats" = "16",
    "output_Covarion_18Cats" = "18", "output_Mk" = "Mk",
    "output_ACRV_4" = "4 ACRV", "output_ACRV_8" = "8 ACRV"
  )
  factor_levels <- c("18", "16", "14", "12", "10", "8", "6", "5", "4", "3", "2", "8 ACRV", "4 ACRV", "Mk")
  combined_data$output_dir <- factor(combined_data$output_dir, levels = factor_levels)

  # this will highlight the one that is chosen
  combined_data <- combined_data %>%
    mutate(is_highlighted = (output_dir == highlight_cat))

  # Plot 1: Standard Deviation
  plot1 <-
    combined_data %>%
    filter(output_dir %in% c("18", "16", "14", "12", "10", "8", "6", "4", "2")) %>%
    ggplot(aes(x = sd_category, y = output_dir, fill = is_highlighted)) +
    stat_binline(scale = 1, alpha = 0.7, bins = 30, draw_baseline = FALSE) +
    labs(title = "Standard Deviation", x = "Standard Deviation", y = "Number of Rate Categories") +
    scale_x_continuous(limits = c(-0.20, 20), expand = c(0, 0), breaks = seq(0, 20, 5)) +
    scale_fill_manual(values = c("TRUE" = "navy", "FALSE" = "grey80")) +
    theme_pilot(axis_title_size = 25, axis_text_size = 18, title_size = 25) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(clip = "off")

  # Plot 2: Switching Rate
  plot2 <-
    combined_data %>%
    filter(output_dir %in% c("18", "16", "14", "12", "10", "8", "6", "4", "2")) %>%
    ggplot(aes(x = global_switch_rate, y = output_dir, fill = is_highlighted)) +
    stat_binline(bins = 30, scale = 1, alpha = 0.7, draw_baseline = FALSE) +
    labs(title = "Switching Rate", x = "Switching Rate", y = "") +
    scale_fill_manual(values = c("TRUE" = "navy", "FALSE" = "grey80")) +
    theme_pilot(axis_title_size = 25, axis_text_size = 18, title_size = 25) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(clip = "off")

  # for tree length deviation from the Mk model
  all_output_dirs <- c(op_mk, OUTPUT_DIRS)
  tree_length_tl <- list()

  for (output_dir in all_output_dirs) {
    logfiles <- list.files(file.path("output_Empirical_Partitioned", output_dir), pattern = "\\.log$", full.names = TRUE)
    logfile <- logfiles[grepl(filename, logfiles) & !grepl("run", logfiles)]
    if (length(logfile) > 0) {
      data <- read.delim(logfile)
      tree_length_tl[[output_dir]] <- data$tree_length
    }
  }

  min_rows <- min(sapply(tree_length_tl, function(x) length(x)))
  tree_length_tl <- lapply(tree_length_tl, function(x) x[seq_len(min_rows)])
  tl_df <- bind_cols(tree_length_tl)

  percent_diff_tl <- sapply(tl_df[-1], function(x) ((x - tl_df$output_Mk) / tl_df$output_Mk) * 100)

  colnames(percent_diff_tl) <- c("4 ACRV", "8 ACRV", "2", "3", "4", "5", "6", "8", "10", "12", "14", "16", "18")
  percent_diff_tl <- melt(percent_diff_tl)
  percent_diff_tl$Var2 <- factor(percent_diff_tl$Var2, levels = c("18", "16", "14", "12", "10", "8", "6", "5", "4", "3", "2", "8 ACRV", "4 ACRV"))

  percent_diff_tl <- percent_diff_tl %>%
    mutate(is_highlighted = (Var2 == highlight_cat))

  # Plot 3: Tree Length Relative to Mk model
  plot3 <-
    percent_diff_tl %>%
    filter(Var2 %in% c("18", "16", "14", "12", "10", "8", "6", "4", "2", "4 ACRV")) %>%
    ggplot() +
    stat_binline(aes(x = value, y = Var2, fill = is_highlighted), bins = 40, scale = 1, alpha = 0.7, draw_baseline = FALSE) +
    scale_x_continuous(limits = c(-200, 300)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
    scale_fill_manual(values = c("TRUE" = "navy", "FALSE" = "grey80")) +
    labs(x = "Percentage difference (tree length)", y = "", title = "Tree Length relative to Mk model") +
    geom_text(
      aes(x = -190, y = Var2, label = Var2),
      hjust = 0, color = "#404040", nudge_y = 0.5, check_overlap = TRUE, family = "Montserrat", size = 18, size.unit = "pt"
    ) +
    theme_pilot(axis_title_size = 25, axis_text_size = 18, title_size = 25) +
    theme(axis.text.y = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

  final_plot <- grid.arrange(
    plot1, plot2, plot3,
    ncol = 3, nrow = 1,
    widths = c(1, 1, 1)
  )

  return(final_plot)
}

# Generate plots
rays_plot <- generate_plots(filename = "Dasyatoidea_MarramaEtAl2023_All", highlight_cat = "12")
sharks_plot <- generate_plots(filename = "Shirai_1996a", highlight_cat = "18")

# labels a) and b) using arrangeGrob with top titles
rays_labeled <- arrangeGrob(rays_plot, top = textGrob("a) Rays Dataset", gp = gpar(fontsize = 28, fontfamily = "Montserrat"), x = 0, hjust = 0))
sharks_labeled <- arrangeGrob(sharks_plot, top = textGrob("b) Sharks Dataset", gp = gpar(fontsize = 28, fontfamily = "Montserrat"), x = 0, hjust = 0))

final_combined_plot <- grid.arrange(rays_labeled, sharks_labeled, nrow = 2)


output_combined_filename <- file.path("Plots", "Empirical_Partitioned", "Rays_Sharks_Posterior.pdf")

ggsave(final_combined_plot, file = output_combined_filename, width = 20, height = 16, dpi = 450, device = cairo_pdf, create.dir = T)

message(paste("\nScript completed successfully. Final plot saved to:", output_combined_filename))
