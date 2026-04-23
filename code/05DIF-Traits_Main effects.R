############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Trait distributions for significant main effects
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Import data
dat <- read.csv(input_file, header = TRUE)

# Summary function: mean ± SD
mean_sd <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  c(y = m, ymin = m - s, ymax = m + s)
}

# Colors
col_intensity <- c("Wet" = "#1994E5", "Dry" = "#E7782F")
col_frequency <- c("Low" = "#E5D019", "High" = "#B38AF1")
col_composition <- c("Conspecific" = "#D95F8D", "Heterospecific" = "#33A27F")

# Theme
theme_traits <- theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 12, color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "gray90", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        legend.position = "none")

# Helper: create parsed trait labels
trait_labels <- c("N" = "bold(Plant ~ number ~ per ~ pot)",
                  "H" = "bold(Plant ~ height ~ (cm))",
                  "R_S" = "bold(Root/Shoot)",
                  "SLA" = "bold(Specific ~ leaf ~ area ~ (mm^{2}/mg))",
                  "LDMC" = "bold(Leaf ~ dry ~ matter ~ content ~ (mg/g))",
                  "LNC" = "bold(Leaf ~ nitrogen ~ content ~ '(%)')",
                  "Leaf_CN" = "bold(Leaf ~ C:N)",
                  "Chl" = "bold(Leaf ~ chlorophyll ~ content ~ (mg/m^{2}))",
                  "PSIIeff" = "bold(Photosystem ~ II ~ efficiency)",
                  "Ψosm" = "bold(Osmotic ~ potential ~ (MPa))",
                  "SS" = "bold(Stomatal ~ size ~ (µm))",
                  "SD" = "bold(Stomatal ~ density ~ (no./mm^{2}))",
                  "RD" = "bold(Root ~ diameter ~ (mm))",
                  "RDMC" = "bold(Root ~ dry ~ matter ~ content ~ (mg/g))",
                  "SRL" = "bold(Specific ~ root ~ length ~ (m/g))",
                  "SRA" = "bold(Specific ~ root ~ area ~ (cm^{2}/g))",
                  "RTD" = "bold(Root ~ tissue ~ density ~ (g/cm^{3}))",
                  "RNC" = "bold(Root ~ nitrogen ~ content ~ '(%)')",
                  "Root_CN" = "bold(Root ~ C:N)",
                  "TB" = "bold(Total ~ biomass ~ per ~ individual ~ (g))")

# Helper: prepare data in long format
prepare_trait_data <- function(data, traits, group_var) {
  data %>%
    pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "value") %>%
    mutate(Trait = factor(Trait, levels = traits),
           Label = recode(as.character(Trait), !!!trait_labels),
           !!group_var := factor(.data[[group_var]], levels = unique(.data[[group_var]])))
}

# Helper: generate violin plot
plot_main_effect <- function(data, x_var, fill_colors, x_label, title = NULL, ncol = 4) {
  ggplot(data, aes(x = .data[[x_var]], y = value, fill = .data[[x_var]])) +
    geom_violin(trim = TRUE, alpha = 0.9) +
    facet_wrap(~ Label, labeller = label_parsed, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = fill_colors) +
    stat_compare_means(method = "t.test", label = "p.signif", label.x.npc = "middle", label.y.npc = "bottom",
                       size = 6, symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                    symbols = c("***", "**", "*", "ns"))) +
    stat_summary(fun.data = mean_sd, color = "white", linewidth = 1, size = 1) +
    labs(x = x_label, y = "Trait value", title = title) +
    theme_traits
}

# Intensity
traits_intensity <- c("N", "H", "R_S", "SLA", "LDMC", "LNC", "Leaf_CN",
                      "Chl", "PSIIeff", "Ψosm", "SS", "SD",
                      "RD", "RDMC", "SRA", "RTD", "RNC", "Root_CN")

dat_intensity <- prepare_trait_data(dat, traits_intensity, "Intensity") %>%
  mutate(Intensity = factor(Intensity, levels = c("Wet", "Dry")))

p_intensity <- plot_main_effect(dat_intensity, x_var = "Intensity",
                                fill_colors = col_intensity,
                                x_label = "Watering intensity")

# Frequency
traits_frequency <- c("RDMC", "SRL", "SRA", "RTD")

dat_frequency <- prepare_trait_data(dat, traits_frequency, "Frequency") %>%
  mutate(Frequency = factor(Frequency, levels = c("Low", "High")))

p_frequency <- plot_main_effect(dat_frequency, x_var = "Frequency",
                                fill_colors = col_frequency,
                                x_label = "Watering frequency")

# Composition
traits_composition <- c("N", "SLA", "LDMC", "RDMC", "RNC", "Root_CN")

dat_composition <- prepare_trait_data(dat, traits_composition, "Composition") %>%
  mutate(Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")))

p_composition <- plot_main_effect(dat_composition, x_var = "Composition",
                                  fill_colors = col_composition,
                                  x_label = "Plant composition")

# Combine panels
p_main_effects <- ggarrange(p_intensity, p_frequency, p_composition, ncol = 1,
                            heights = c(3.8, 1, 1.7), labels = c("(a)", "(b)", "(c)"))
p_main_effects

# Save figure - Figure S4
ggsave(filename = file.path(output_dir, "Traits_main_effects_significant.pdf"), plot = p_main_effects,
       width = 14, height = 20, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "Traits_main_effects_significant.png"), plot = p_main_effects,
       width = 14, height = 20, units = "in", dpi = 300)
