############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Trait responses for significant interactions
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Import data
dat <- read.csv(input_file, header = TRUE)

# Colors
col_intensity <- c("Wet" = "#1994E5", "Dry" = "#E7782F")
col_composition <- c("Conspecific" = "#D95F8D", "Heterospecific" = "#33A27F")

# Theme
theme_interaction <- theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 12, color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "gray90", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        legend.position = "top",
        legend.text = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 12))

# Trait labels
trait_labels <- c("N" = "bold(Plant ~ number ~ per ~ pot)",
                  "TB" = "bold(Total ~ biomass ~ per ~ individual ~ (g))",
                  "R_S" = "bold(Root/Shoot)",
                  "LNC" = "bold(Leaf ~ nitrogen ~ content ~ '(%)')",
                  "Chl" = "bold(Leaf ~ chlorophyll ~ content ~ (mg/m^{2}))",
                  "SS" = "bold(Stomatal ~ size ~ (µm))",
                  "RTD" = "bold(Root ~ tissue ~ density ~ (g/cm^{3}))",
                  "LDMC" = "bold(Leaf ~ dry ~ matter ~ content ~ (mg/g))",
                  "Ψosm" = "bold(Osmotic ~ potential ~ (MPa))")

# Standard error
se_fun <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

# Prepare interaction data
prepare_interaction_data <- function(data, traits, factor_levels = list()) {
  out <- data %>%
    pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "value") %>%
    mutate(Trait = factor(Trait, levels = traits),
           Trait = recode(as.character(Trait), !!!trait_labels))
  
  for (v in names(factor_levels)) {
    out[[v]] <- factor(out[[v]], levels = factor_levels[[v]])
  }
  
  out
}

# Add simple-effect stars
get_star_labels <- function(data, group_var, compare_var) {
  compare_levels <- levels(data[[compare_var]])
  
  pvals <- data %>%
    group_by(Trait, .data[[group_var]]) %>%
    group_modify(~ {
      x <- .x$value[.x[[compare_var]] == compare_levels[1]]
      y <- .x$value[.x[[compare_var]] == compare_levels[2]]
      p <- if (length(na.omit(x)) > 1 && length(na.omit(y)) > 1) {
        t.test(x, y)$p.value
      } else {
        NA_real_
      }
      tibble(p = p)
    }) %>%
    ungroup() %>%
    mutate(label = case_when(is.na(p)  ~ "ns", p < 0.001 ~ "***",
                             p < 0.01  ~ "**", p < 0.05  ~ "*",
                             TRUE      ~ "ns"))
  
  summary_means <- data %>%
    group_by(Trait, .data[[group_var]], .data[[compare_var]]) %>%
    summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Trait, .data[[group_var]]) %>%
    summarise(y_mid = mean(mean, na.rm = TRUE), .groups = "drop")
  
  y_span <- data %>%
    group_by(Trait) %>%
    summarise(y_rng = diff(range(value, na.rm = TRUE)), .groups = "drop")
  
  pvals %>%
    left_join(summary_means, by = c("Trait", group_var)) %>%
    left_join(y_span, by = "Trait") %>%
    mutate(y_star = y_mid + 0.04 * ifelse(is.finite(y_rng), y_rng, 0))
}

# Plot interaction panel
plot_interaction <- function(data, x_var, color_var, colors, title, x_label, legend_label, ncol = 3) {
  summary_data <- data %>%
    group_by(Trait, .data[[x_var]], .data[[color_var]]) %>%
    summarise(mean = mean(value, na.rm = TRUE), se = se_fun(value), .groups = "drop")
  
  ggplot(summary_data, aes(x = .data[[x_var]], y = mean, color = .data[[color_var]], group = .data[[color_var]])) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12) +
    facet_wrap(~ Trait, labeller = label_parsed, scales = "free_y", ncol = ncol) +
    scale_color_manual(values = colors) +
    labs(x = x_label, y = "Trait value", color = legend_label, title = title) +
    theme_interaction
}

# Intensity × Frequency
dat_if <- prepare_interaction_data(dat, traits = c("R_S", "LNC", "Chl"),
                                   factor_levels = list(Intensity = c("Wet", "Dry"),
                                                        Frequency = c("Low", "High")))

p_if <- plot_interaction(dat_if, x_var = "Frequency", color_var = "Intensity",
                         colors = col_intensity, title = "Watering intensity × Watering frequency",
                         x_label = "Watering frequency", legend_label = "Watering intensity")

stars_if <- get_star_labels(dat_if, group_var = "Frequency", compare_var = "Intensity")

p_if <- p_if +
  geom_text(data = stars_if, aes(x = Frequency, y = y_star, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom")

# Intensity × Composition
dat_ic <- prepare_interaction_data(dat, traits = c("N", "SS", "RTD"),
                                   factor_levels = list(Intensity = c("Wet", "Dry"),
                                                        Composition = c("Conspecific", "Heterospecific")))

p_ic <- plot_interaction(dat_ic, x_var = "Composition", color_var = "Intensity",
                         colors = col_intensity, title = "Watering intensity × Plant composition",
                         x_label = "Plant composition", legend_label = "Watering intensity")

stars_ic <- get_star_labels(dat_ic, group_var = "Composition", compare_var = "Intensity")

p_ic <- p_ic +
  geom_text(data = stars_ic, aes(x = Composition, y = y_star, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom")

# Frequency × Composition
dat_fc <- prepare_interaction_data(dat, traits = c("N", "TB", "LDMC", "Ψosm", "RTD"),
                                   factor_levels = list(Frequency = c("Low", "High"),
                                                        Composition = c("Conspecific", "Heterospecific")))

p_fc <- plot_interaction(dat_fc, x_var = "Frequency", color_var = "Composition",
                         colors = col_composition, title = "Watering frequency × Plant composition",
                         x_label = "Watering frequency", legend_label = "Plant composition")

stars_fc <- get_star_labels(dat_fc, group_var = "Frequency", compare_var = "Composition")

p_fc <- p_fc +
  geom_text(data = stars_fc, aes(x = Frequency, y = y_star, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(0.97, 0.05),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.box.background = element_rect(color = "black", linewidth = 0.3),
        legend.key = element_blank())

# Combine panels
p_interactions <- ggarrange(p_if, p_ic, p_fc, ncol = 1, heights = c(1, 1, 1.42),
                            labels = c("(a)", "(b)", "(c)"))
p_interactions

# Save figure - Figure S5
ggsave(filename = file.path(output_dir, "Traits_interactions_lines.pdf"), plot = p_interactions,
       width = 10, height = 12, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "Traits_interactions_lines.png"), plot = p_interactions,
       width = 10, height = 12, units = "in", dpi = 300)

