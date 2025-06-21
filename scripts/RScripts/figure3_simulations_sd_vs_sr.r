#!/user/bin/env Rscript
# Rscript to plot simulation output
# Authors: Basanta Khakurel and Sebastian Höhna
# date: 2025-05-30

library(tidyverse)
library(patchwork)
library(ggbreak)
library(pilot)
library(extrafont)
library(ggridges)
library(viridis)
library(pals)
library(ggview)

set_pilot_family(family = "Montserrat")

H <- abs(log(10) / (qnorm(0.975) - qnorm(0.025)))
true_sd <- 2 * H
num_states <- 4
num_cats <- 4

original_tl <- 10.54336

Models <- c("0.005_SR", "0.05_SR", "0.5_SR", "1_SR", "5_SR", "Mk", "MkG")

dirs <- lapply(Models, function(model) paste0("Results_Output_Simulated_Data_", model, "_", num_states, "_States_", num_cats, "_Cats"))

sims_data_list <- list()

for (model in Models) {
  model_idx <- 1
  file_path <- paste0("Simulations/Results_Output_Simulated_Data_", model, "_", num_states, "_States_", num_cats, "_Cats/morpho_sims_median.csv")
  temp_data <- read.csv(file_path)
  temp_data$Model <- model
  temp_data$pct_diff <- ((temp_data$tree_length - original_tl) / original_tl) * 100

  sims_data_list[[model]] <- temp_data
  model_idx <- model_idx + 1
}

combined_data <- bind_rows(sims_data_list)

true_sr_values <- data.frame(
  Model = factor(c("Mk", "MkG", "0.005_SR", "0.05_SR", "0.5_SR", "5_SR"),
    levels = c("Mk", "MkG", "0.005_SR", "0.05_SR", "0.5_SR", "5_SR")
  ),
  true_rate = c(0, 0, 0.005, 0.05, 0.5, 5)
)

true_sd_values <- data.frame(
  Model = factor(c("Mk", "MkG", "0.005_SR", "0.05_SR", "0.5_SR", "5_SR"),
    levels = c("Mk", "MkG", "0.005_SR", "0.05_SR", "0.5_SR", "5_SR")
  ),
  true_rate = c(NA, true_sd, true_sd, true_sd, true_sd, true_sd)
)

# plot switching rates vs standard deviation
sr_v_sd <- combined_data %>%
  filter(Model %in% c("0.005_SR", "0.05_SR", "0.5_SR", "5_SR", "Mk", "MkG")) %>%
  mutate(Model = factor(Model, levels = c("Mk", "MkG", "0.005_SR", "0.05_SR", "0.5_SR", "5_SR"))) %>%
  ggplot(aes(x = sd_category)) +
  geom_point(aes(y = global_switch_rate), alpha = 0.6, size = 3, shape = 21, color = "black", fill = "black") +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.005, 0.05, 0.5, 5, 50, 500, 5000),
    labels = c("0.005", "0.05", "0.5", "5", "50", "500", "5000")
  ) +
  scale_x_continuous(limits = c(0, NA)) +
  geom_vline(
    data = true_sd_values,
    aes(
      xintercept = true_rate
    ),
    linetype = "dashed",
    color = "red",
    linewidth = 0.9
  ) +
  geom_hline(
    data = true_sr_values,
    aes(yintercept = true_rate),
    linetype = "dashed",
    color = "blue",
    linewidth = 0.9
  ) +
  labs(x = "Standard Deviation (σ)", y = "Switching Rate (δ)") +
  facet_wrap(
    ~Model,
    nrow = 1,
    labeller = as_labeller(c("0.005_SR" = "0.005", "0.05_SR" = "0.05", "0.5_SR" = "0.5", "5_SR" = "5", "Mk" = "0\n(Mk)", "MkG" = "0\n(Mk + ACRV)"))
  ) +
  theme_pilot(facet_title_size = 18, axis_title_size = 18, axis_text_size = 15, axes = "", grid = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "none",
    strip.text = element_text(face = "plain"),
    strip.background = element_blank()
  ) +
  coord_cartesian(clip = "off")

ggsave("Plots/Simulations/sr_vs_sd_sims.pdf", plot = sr_v_sd, device = cairo_pdf, width = 12, height = 4, dpi = 450, limitsize = FALSE, create.dir = TRUE, bg = "white")
