############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Heatmaps of variance explained
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(ggh4x)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file_ifc  <- "data/data_processed/Variance-I*F*C.csv" # sorted results from script 02
input_file_ifb  <- "data/data_processed/Variance-I*F*B.csv" # sorted results from script 02
input_file_ifcb <- "data/data_processed/Variance-I*F*C+B.csv" # sorted results from script 02
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Import data
data_ifc  <- read.csv(input_file_ifc, header = TRUE)
data_ifb  <- read.csv(input_file_ifb, header = TRUE)
data_ifcb <- read.csv(input_file_ifcb, header = TRUE)

# Trait and category order
category_levels <- c("Growth", "Leaf", "Physiological", "Root")

trait_order <- c("N", "H", "SB", "TB", "R/S",
                 "SLA", "LDMC", "LNC", "Leaf C:N",
                 "Chl", "PSIIeff", "Ψosm", "SS", "SD",
                 "RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root C:N")

facet_colors <- c(Growth = "gray90", Leaf = "#C7E9C0",
                  Physiological = "#D9EDF7", Root = "#FCF8E3")

# Plot theme
theme_heatmap <- theme_classic() +
  theme(axis.text = element_text(face = "bold", size = 12, color = "black"),
        axis.title = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(face = "bold", size = 12),
        strip.background = element_blank(),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 10))

# Rename trait labels
rename_traits <- function(x) {
  case_match(x, "R_S" ~ "R/S", "Psi_osm" ~ "Ψosm",
             "Leaf_CN" ~ "Leaf C:N", "Root_CN" ~ "Root C:N", .default = x)
}

# Prepare one dataset for plotting
prepare_heatmap_data <- function(data, factor_labels) {
  df_sig <- data %>%
    mutate(Trait = rename_traits(Trait),
           Factor = dplyr::recode(Factor, !!!factor_labels),
           significant = p < 0.05) %>%
    filter(significant)
  
  avg_row <- df_sig %>%
    group_by(Factor) %>%
    summarise(Variance_Explained = mean(Variance_Explained, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Trait = "Average", Category = "Average")
  
  heat_df <- df_sig %>%
    select(Trait, Category, Factor, Variance_Explained) %>%
    bind_rows(avg_row) %>%
    mutate(label_text = sprintf("%.1f%%", round(Variance_Explained, 1)))
  
  list(avg = heat_df %>%
         filter(Trait == "Average"),
       traits = heat_df %>%
         filter(Trait != "Average"))
}

# Format plotting data
format_heatmap_panels <- function(heat_data, factor_levels) {
  heat_avg <- heat_data$avg %>%
    mutate(Factor = factor(Factor, levels = factor_levels),
           Trait = factor(Trait, levels = "Average"))
  
  heat_traits <- heat_data$traits %>%
    mutate(Factor = factor(Factor, levels = factor_levels),
           Category = factor(Category, levels = category_levels),
           Trait = factor(Trait, levels = trait_order)) %>%
    mutate(Trait = fct_rev(Trait))
  
  list(avg = heat_avg, traits = heat_traits)
}

# Plot average panel
plot_average_panel <- function(data) {
  ggplot(data, aes(x = Factor, y = Trait, fill = Variance_Explained)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label_text, color = Variance_Explained), size = 4) +
    scale_color_gradient(low = "black", high = "white", guide = "none") +
    scale_fill_viridis_c(option = "B", direction = -1, limits = c(0, 100)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_heatmap +
    theme(legend.position = "none",
          strip.text = element_blank())
}

# Plot trait panel with facets
plot_trait_panel <- function(data, title = NULL) {
  ggplot(data, aes(x = Factor, y = Trait, fill = Variance_Explained)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label_text, color = Variance_Explained), size = 4) +
    scale_color_gradient(low = "black", high = "white", guide = "none") +
    scale_fill_viridis_c(option = "B", direction = -1, limits = c(0, 100),
                         name = "Contribution to model-explained variance (%)") +
    scale_x_discrete(position = "top") +
    labs(x = "Factor", y = "Trait", title = title) +
    facet_grid2(rows = vars(Category), scales = "free_y", space = "free_y",
                strip = strip_themed(background_y = elem_list_rect(fill = facet_colors, color = NA),
                                     text_y = elem_list_text(face = "bold", size = 12))) +
    theme_heatmap +
    theme(plot.title = element_text(face = "bold", size = 12, color = "#4A90E2"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5, linetype = "dotted"),
          panel.spacing.y = unit(0.6, "lines"))
}

# Plot trait panel without facets
plot_trait_panel_simple <- function(data, title = NULL) {
  ggplot(data, aes(x = Factor, y = Trait, fill = Variance_Explained)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label_text, color = Variance_Explained), size = 4) +
    scale_color_gradient(low = "black", high = "white", guide = "none") +
    scale_fill_viridis_c(option = "B", direction = -1, limits = c(0, 100),
                         name = "Contribution to model-explained variance (%)") +
    scale_x_discrete(position = "top") +
    labs(x = "Factor", y = "Trait", title = title) +
    theme_heatmap +
    theme(plot.title = element_text(face = "bold", size = 12, color = "#4A90E2"))
}

# Factor labels for each model
factor_labels_ifc <- c("Intensity" = "Intensity (I)",
                       "Frequency" = "Frequency (F)",
                       "Composition" = "Composition (C)",
                       "Intensity:Frequency" = "I × F",
                       "Intensity:Composition" = "I × C",
                       "Frequency:Composition" = "F × C",
                       "Intensity:Frequency:Composition" = "I × F × C")

factor_labels_ifb <- c("Intensity" = "Intensity (I)",
                       "Frequency" = "Frequency (F)",
                       "TB_pot_c" = "Biomass (B)",
                       "Intensity:Frequency" = "I × F",
                       "Intensity:TB_pot_c" = "I × B",
                       "Frequency:TB_pot_c" = "F × B",
                       "Intensity:Frequency:TB_pot_c" = "I × F × B")

factor_labels_ifcb <- c("Intensity" = "Intensity (I)",
                        "Frequency" = "Frequency (F)",
                        "Composition" = "Composition (C)",
                        "TB_pot_c" = "Biomass (B)",
                        "Intensity:Frequency" = "I × F",
                        "Intensity:Composition" = "I × C",
                        "Frequency:Composition" = "F × C",
                        "Intensity:Frequency:Composition" = "I × F × C")

# Factor order for each model
factor_levels_ifc <- c("Intensity (I)", "Frequency (F)", "Composition (C)",
                       "I × F", "I × C", "F × C")

factor_levels_ifb <- c("Frequency (F)", "Biomass (B)",
                       "I × F", "I × B", "F × B", "I × F × B")

factor_levels_ifcb <- c("Intensity (I)", "Frequency (F)", "Composition (C)",
                        "Biomass (B)", "I × F", "F × C")

# Prepare data
heat_ifc <- prepare_heatmap_data(data_ifc, factor_labels_ifc)
heat_ifb <- prepare_heatmap_data(data_ifb, factor_labels_ifb)
heat_ifcb <- prepare_heatmap_data(data_ifcb, factor_labels_ifcb)

panels_ifc <- format_heatmap_panels(heat_ifc, factor_levels_ifc)
panels_ifb <- format_heatmap_panels(heat_ifb, factor_levels_ifb)
panels_ifcb <- format_heatmap_panels(heat_ifcb, factor_levels_ifcb)

# Figure 1: main IFC model
p_avg_ifc <- plot_average_panel(panels_ifc$avg)
p_traits_ifc <- plot_trait_panel(panels_ifc$traits)

p_final_ifc <- ggarrange(p_traits_ifc, p_avg_ifc, ncol = 1, heights = c(12.5, 1))
p_final_ifc

# Save
ggsave(filename = file.path(output_dir, "Heatmap_variance_explained_IFC.pdf"), plot = p_final_ifc,
       width = 10, height = 11, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "Heatmap_variance_explained_IFC.png"), plot = p_final_ifc,
       width = 10, height = 11, units = "in", dpi = 300)


# Figure S6: alternative models
p_traits_ifb <- plot_trait_panel_simple(panels_ifb$traits, title = "I × F × B model")

p_avg_ifcb <- plot_average_panel(panels_ifcb$avg)
p_traits_ifcb <- plot_trait_panel(panels_ifcb$traits, title = "I × F × C + B model")

p_final_ifcb <- ggarrange(p_traits_ifcb, p_avg_ifcb, ncol = 1, heights = c(11, 1))

p_final_alt <- ggarrange(p_final_ifcb, NULL, p_traits_ifb, ncol = 1,
                         heights = c(4, 0.2, 1), labels = c("(a)", " ", "(b)"))
p_final_alt

# Save
ggsave(filename = file.path(output_dir, "Heatmap_variance_explained_alternative_models.pdf"), plot = p_final_alt,
       width = 10, height = 12, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "Heatmap_variance_explained_alternative_models.png"), plot = p_final_alt,
       width = 10, height = 12, units = "in", dpi = 300)
