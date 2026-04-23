############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Standardized model-based trait responses
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(purrr)
library(ggplot2)
library(emmeans)
library(tibble)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Import data
dat <- read.csv(input_file, header = TRUE) %>%
  mutate(Intensity = factor(Intensity, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")))

# Set contrasts for Type III interpretation
options(contrasts = c("contr.sum", "contr.poly"))

# Traits included in the figure
traits <- c("N", "H", "SB", "TB", "R_S",
            "SLA", "LDMC", "LNC", "Leaf_CN",
            "Chl", "PSIIeff", "Ψosm", "SS", "SD",
            "RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root_CN")

# Display labels for traits
trait_label <- function(x) {
  case_match(x, "R_S" ~ "R/S", "Leaf_CN" ~ "Leaf C:N",
             "Root_CN" ~ "Root C:N", .default = x)
}

# Trait categories
trait_category <- tribble(~Trait,      ~Category,
                          "N",         "Growth",
                          "H",         "Growth",
                          "SB",        "Growth",
                          "TB",        "Growth",
                          "R/S",       "Growth",
                          "SLA",       "Leaf",
                          "LDMC",      "Leaf",
                          "LNC",       "Leaf",
                          "Leaf C:N",  "Leaf",
                          "Chl",       "Physiological",
                          "PSIIeff",   "Physiological",
                          "Ψosm",      "Physiological",
                          "SS",        "Physiological",
                          "SD",        "Physiological",
                          "RD",        "Root",
                          "RDMC",      "Root",
                          "SRL",       "Root",
                          "SRA",       "Root",
                          "RTD",       "Root",
                          "RNC",       "Root",
                          "Root C:N",  "Root")

category_levels <- c("Growth", "Leaf", "Physiological", "Root")

trait_levels <- c("N", "H", "SB", "TB", "R/S",
                  "SLA", "LDMC", "LNC", "Leaf C:N",
                  "Chl", "PSIIeff", "Ψosm", "SS", "SD",
                  "RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root C:N")

# Panel titles
panel_levels <- c("Wet \u2190 Intensity \u2192 Dry",
                  "Low \u2190 Frequency \u2192 High",
                  "Conspecific \u2190 Composition \u2192 Heterospecific")

# Panel colors
col_intensity <- c(Wet = "#1994E5", Dry = "#E7782F")
col_frequency <- c(Low = "#E5D019", High = "#B38AF1")
col_composition <- c(Conspecific = "#D95F8D", Heterospecific = "#33A27F")

# Helper: extract contrast estimates and 95% CI
get_contrast_table <- function(emm, contrast_list) {
  contrast(emm, method = contrast_list) %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame()
}

# Extract standardized main-effect contrasts for one trait
get_standardized_contrasts <- function(trait) {
  df <- dat %>%
    select(Intensity, Frequency, Composition, all_of(trait)) %>%
    filter(!is.na(.data[[trait]])) %>%
    mutate(y = as.numeric(scale(.data[[trait]])))
  
  model <- lm(y ~ Intensity * Frequency * Composition, data = df)
  
  contrast_intensity <- get_contrast_table(emmeans(model, ~ Intensity),
                                           list("Dry - Wet" = c(-1, 1))) %>%
    mutate(Panel = panel_levels[1])
  
  contrast_frequency <- get_contrast_table(emmeans(model, ~ Frequency),
                                           list("High - Low" = c(-1, 1))) %>%
    mutate(Panel = panel_levels[2])
  
  contrast_composition <- get_contrast_table(emmeans(model, ~ Composition),
                                             list("Heterospecific - Conspecific" = c(-1, 1))) %>%
    mutate(Panel = panel_levels[3])
  
  bind_rows(contrast_intensity, contrast_frequency, contrast_composition) %>%
    rename(lower = lower.CL, upper = upper.CL) %>%
    mutate(Trait = trait_label(trait),
           sig = ifelse(lower * upper > 0, "sig", "ns")) %>%
    select(Trait, Panel, estimate, lower, upper, sig)
}

# Compile results for all traits
trait_responses <- map_df(traits, get_standardized_contrasts) %>%
  mutate(Trait = factor(Trait, levels = trait_levels),
         Panel = factor(Panel, levels = panel_levels)) %>%
  left_join(trait_category, by = "Trait") %>%
  mutate(Category = factor(Category, levels = category_levels),
         trait_index = match(Trait, trait_levels))

# Background shading for each panel
background_map <- tibble(Panel = factor(panel_levels, levels = panel_levels),
                         left_fill = c(col_intensity["Wet"], col_frequency["Low"], col_composition["Conspecific"]),
                         right_fill = c(col_intensity["Dry"], col_frequency["High"], col_composition["Heterospecific"]))

rect_left <- background_map %>%
  transmute(Panel, xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = left_fill)

rect_right <- background_map %>%
  transmute(Panel, xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = right_fill)

# Plot
p_trait <- ggplot(trait_responses, aes(x = estimate, y = trait_index)) +
  geom_rect(data = rect_left,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_rect(data = rect_right,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            alpha = 0.2, inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_pointrange(aes(xmin = lower, xmax = upper, alpha = sig), linewidth = 0.6) +
  scale_alpha_manual(values = c(sig = 1, ns = 0.25), guide = "none") +
  facet_grid(Category ~ Panel, scales = "free_y", space = "free_y") +
  scale_y_reverse(breaks = seq_along(trait_levels), labels = trait_levels) +
  labs(x = "Standardized trait response (SD units)", y = "Trait") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.line = element_blank(),
        strip.text.x = element_text(face = "bold", size = 12),
        strip.text.y = element_text(face = "bold", size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.spacing.x = unit(1.2, "lines"),
        panel.spacing.y = unit(0.8, "lines"))
p_trait

# Save figure - Figure 1
ggsave(filename = file.path(output_dir, "Traits_standardized_responses.pdf"), plot = p_trait,
       width = 15, height = 11, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "Traits_standardized_responses.png"), plot = p_trait,
       width = 15, height = 11, units = "in", dpi = 300)

