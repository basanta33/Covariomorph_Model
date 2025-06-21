#!/user/bin/env Rscript
# Rscript to plot empirical results
# Authors: Basanta Khakurel and Sebastian Höhna
# date: 2025-05-30

library(ggview, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggridges, quietly = TRUE)
library(pilot, quietly = TRUE)
library(reshape2)
library(ggview)
library(extrafont)

set_pilot_family(family = "Montserrat")

output_dirs <- c(
  "output_Covarion_2Cats", "output_Covarion_4Cats", "output_Covarion_6Cats",
  "output_Covarion_8Cats", "output_Covarion_10Cats", "output_Covarion_12Cats"
)

data_list <- list()

for (dir in output_dirs) {
  file_path <- paste0("Results/output_Partitioned/results_", dir, "/", dir, ".csv")
  temp_data <- read.csv(file_path)
  temp_data$Directory <- dir
  data_list[[dir]] <- temp_data
}

combined_data <- bind_rows(data_list) %>%
  mutate(Directory = factor(Directory,
    levels = output_dirs,
    labels = c("m = 2", "m = 4", "m = 6", "m = 8", "m = 10", "m = 12")
  ))

combined_data <- combined_data %>%
  dplyr::mutate(
    highlight = dplyr::case_when(
      filename == "Shirai_1996a" ~ "Sharks",
      filename == "Dasyatoidea_MarramaEtAl2023_All" ~ "Rays",
      TRUE ~ "other"
    ),
    highlight_legend = dplyr::case_when(
      highlight %in% c("Sharks", "Rays") ~ highlight,
      TRUE ~ NA_character_
    )
  )

p <- ggplot(combined_data, aes(x = sd_category, y = global_switch_rate)) +
  geom_point(aes(color = highlight_legend, fill = highlight_legend, shape = highlight_legend, alpha = ifelse(highlight %in% c("Sharks", "Rays"), 1, 0.8)), size = 3.5) +
  scale_color_manual(
    values = c("Rays" = "red", "Sharks" = "red"), na.value = "black",
    breaks = c("Rays", "Sharks")
  ) +
  scale_fill_manual(
    values = c("Rays" = "thistle", "Sharks" = "cyan"), na.value = "black",
    breaks = c("Rays", "Sharks")
  ) +
  scale_shape_manual(
    values = c("Rays" = 24, "Sharks" = 23), na.value = 21,
    breaks = c("Rays", "Sharks")
  ) +
  scale_alpha_continuous(range = c(0.6, 1), guide = "none") +
  facet_wrap(~Directory, nrow = 2) +
  labs(
    y = "Switching Rate (δ)",
    x = "Standard Deviation (σ)",
    color = "", fill = "", shape = ""
  ) +
  theme_pilot(facet_title_size = 20, axis_title_size = 20, axis_text_size = 18, axes = "", grid = "hv") +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.005, 0.05, 0.5, 5, 50, 500, 5000),
    labels = c("0.005", "0.05", "0.5", "5", "50", "500", "5000"),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  theme(strip.text = element_text(face = "plain"), strip.background = element_blank(), legend.position = "bottom") +
  coord_cartesian(clip = "off")

ggsave(plot = p, file = "Plots/Empirical_Partitioned/sd_sr_empirical.pdf", device = cairo_pdf, width = 9, height = 6, dpi = 450)

q()
