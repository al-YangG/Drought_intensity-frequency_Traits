############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Heatmaps of variance explained
# Date: June 2026
############################################################

## Load packages
library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(ggh4x)

setwd("~/Desktop/DIF/Traits/Data availability")

## File paths
input_file_base <- "results/tables/ANOVA_Amount_Frequency_Composition.csv"
input_file_bio  <- "results/tables/ANOVA_Amount_Frequency_Composition_Biomass.csv"
input_file_den  <- "results/tables/ANOVA_Amount_Frequency_Composition_InitialDensity.csv"
input_file_full <- "results/tables/ANOVA_Amount_Frequency_Composition_Biomass_InitialDensity.csv"
input_file_sub  <- "results/tables/ANOVA_Amount_Frequency_BiomassSubstitution.csv"

output_dir <- "results/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Import data
data_base <- read.csv(input_file_base, header = TRUE)
data_bio  <- read.csv(input_file_bio,  header = TRUE)
data_den  <- read.csv(input_file_den,  header = TRUE)
data_full <- read.csv(input_file_full, header = TRUE)
data_sub  <- read.csv(input_file_sub,  header = TRUE)

## Trait and category order
category_levels <- c("Growth", "Leaf", "Physiological", "Root")

trait_order <- c("N", "H", "SB", "TB", "R/S",
                 "SLA", "LDMC", "LNC", "Leaf C:N",
                 "Chl", "PSIIeff", "Ψosm", "SS", "SD",
                 "RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root C:N")

facet_colors <- c(Growth = "gray90", Leaf = "#C7E9C0",
                  Physiological = "#D9EDF7", Root = "#FCF8E3")

## Plot theme
theme_heatmap <- theme_classic() +
  theme(axis.text = element_text(face = "bold", size = 11, color = "black"),
        axis.title = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 10))

## Rename trait labels
rename_traits <- function(x) {
  dplyr::recode(x, "R_S" = "R/S", "Psi_osm" = "Ψosm", "Leaf_CN" = "Leaf C:N",
                "Root_CN" = "Root C:N", .default = x)
}

## Add category column
add_category <- function(df) {
  df %>%
    mutate(Category = case_when(Trait %in% c("N", "H", "SB", "TB", "R_S") ~ "Growth",
                                Trait %in% c("SLA", "LDMC", "LNC", "Leaf_CN") ~ "Leaf",
                                Trait %in% c("Chl", "PSIIeff", "Psi_osm", "SS", "SD") ~ "Physiological",
                                Trait %in% c("RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root_CN") ~ "Root",
                                TRUE ~ NA_character_))
}

## Prepare one dataset for plotting
prepare_heatmap_data <- function(data, factor_labels, include_average = TRUE) {
  
  df_sig <- data %>%
    add_category() %>%
    mutate(Trait = rename_traits(Trait), Factor = dplyr::recode(term, !!!factor_labels),
           significant = p < 0.05) %>%
    filter(significant)
  
  heat_df <- df_sig %>%
    select(Trait, Category, Factor, Variance_Explained)
  
  if (include_average) {
    avg_row <- df_sig %>%
      group_by(Factor) %>%
      summarise(Variance_Explained = mean(Variance_Explained, na.rm = TRUE), .groups = "drop") %>%
      mutate(Trait = "Average", Category = "Average")
    
    heat_df <- bind_rows(heat_df, avg_row)
  }
  
  heat_df %>%
    mutate(label_text = sprintf("%.1f%%", round(Variance_Explained, 1)))
}

## Format plotting data
format_heatmap_panel <- function(data, factor_levels, include_average = TRUE) {
  
  if (include_average) {
    heat_avg <- data %>%
      filter(Trait == "Average") %>%
      mutate(Factor = factor(Factor, levels = factor_levels),
             Trait = factor(Trait, levels = "Average"))
    
    heat_traits <- data %>%
      filter(Trait != "Average") %>%
      mutate(Factor = factor(Factor, levels = factor_levels),
             Category = factor(Category, levels = category_levels),
             Trait = factor(Trait, levels = trait_order)) %>%
      mutate(Trait = fct_rev(Trait))
    
    list(avg = heat_avg, traits = heat_traits)
    
  } else {
    data %>%
      mutate(Factor = factor(Factor, levels = factor_levels),
             Category = factor(Category, levels = category_levels),
             Trait = factor(Trait, levels = trait_order)) %>%
      mutate(Trait = fct_rev(Trait))
  }
}

## Plot average panel
plot_average_panel <- function(data) {
  ggplot(data, aes(x = Factor, y = Trait, fill = Variance_Explained)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label_text, color = Variance_Explained), size = 4) +
    scale_color_gradient(low = "black", high = "white", guide = "none") +
    scale_fill_viridis_c(option = "B", direction = -1, limits = c(0, 100)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_heatmap +
    theme(legend.position = "none", strip.text = element_blank())
}

## Plot trait panel with facets
plot_trait_panel <- function(data, title = NULL, show_legend = TRUE) {
  ggplot(data, aes(x = Factor, y = Trait, fill = Variance_Explained)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label_text, color = Variance_Explained), size = 3.5) +
    scale_color_gradient(low = "black", high = "white", guide = "none") +
    scale_fill_viridis_c(option = "B", direction = -1, limits = c(0, 100),
                         name = "Contribution to model-explained variance (%)") +
    scale_x_discrete(position = "top") +
    labs(x = "Factor", y = "Trait", title = title) +
    facet_grid2(rows = vars(Category), scales = "free_y", space = "free_y",
                strip = strip_themed(background_y = elem_list_rect(fill = facet_colors, color = NA),
                                     text_y = elem_list_text(face = "bold", size = 11))) +
    theme_heatmap +
    theme(plot.title = element_text(face = "bold", size = 12, color = "#4A90E2"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4, linetype = "dotted"),
          panel.spacing.y = unit(0.5, "lines"),
          legend.position = ifelse(show_legend, "bottom", "none"))
}


## Factor labels and orders

## Baseline model: A × F × C
factor_labels_base <- c("Amount" = "Amount (A)",
                        "Frequency" = "Frequency (F)",
                        "Composition" = "Composition (C)",
                        "Amount:Frequency" = "A × F",
                        "Amount:Composition" = "A × C",
                        "Frequency:Composition" = "F × C",
                        "Amount:Frequency:Composition" = "A × F × C")

factor_levels_base <- c("Amount (A)", "Frequency (F)", "Composition (C)",
                        "A × F", "A × C", "F × C", "A × F × C")

## Biomass-adjusted model
factor_labels_bio <- c("Amount" = "Amount (A)",
                       "Frequency" = "Frequency (F)",
                       "Composition" = "Composition (C)",
                       "TB_pot_c" = "Biomass (B)",
                       "Amount:Frequency" = "A × F",
                       "Amount:Composition" = "A × C",
                       "Frequency:Composition" = "F × C",
                       "Amount:Frequency:Composition" = "A × F × C")

factor_levels_bio <- c("Amount (A)", "Frequency (F)", "Composition (C)",
                       "Biomass (B)", "A × F", "A × C", "F × C", "A × F × C")

## Density-adjusted model
factor_labels_den <- c("Amount" = "Amount (A)",
                       "Frequency" = "Frequency (F)",
                       "Composition" = "Composition (C)",
                       "Bromus_N_2022_c" = "Density (D)",
                       "Amount:Frequency" = "A × F",
                       "Amount:Composition" = "A × C",
                       "Frequency:Composition" = "F × C",
                       "Amount:Frequency:Composition" = "A × F × C")

factor_levels_den <- c("Amount (A)", "Frequency (F)", "Composition (C)",
                       "Density (D)", "A × F", "A × C", "F × C", "A × F × C")

## Fully adjusted model
factor_labels_full <- c("Amount" = "Amount (A)",
                        "Frequency" = "Frequency (F)",
                        "Composition" = "Composition (C)",
                        "TB_pot_c" = "Biomass",
                        "Bromus_N_2022_c" = "Density",
                        "Amount:Frequency" = "A × F",
                        "Amount:Composition" = "A × C",
                        "Frequency:Composition" = "F × C",
                        "Amount:Frequency:Composition" = "A × F × C")

factor_levels_full <- c("Amount (A)", "Frequency (F)", "Composition (C)",
                        "Biomass", "Density", "A × F", "A × C", "F × C", "A × F × C")

## Biomass substitution model
factor_labels_sub <- c("Amount" = "Amount (A)",
                       "Frequency" = "Frequency (F)",
                       "TB_pot_c" = "Biomass (B)",
                       "Amount:Frequency" = "A × F",
                       "Amount:TB_pot_c" = "A × B",
                       "Frequency:TB_pot_c" = "F × B",
                       "Amount:Frequency:TB_pot_c" = "A × F × B")

factor_levels_sub <- c("Amount (A)", "Frequency (F)", "Biomass (B)",
                       "A × F", "A × B", "F × B", "A × F × B")

## Figure 1: baseline model

heat_base <- prepare_heatmap_data(data_base, factor_labels_base, include_average = TRUE)
panels_base <- format_heatmap_panel(heat_base, factor_levels_base, include_average = TRUE)

p_avg_base <- plot_average_panel(panels_base$avg)
p_traits_base <- plot_trait_panel(panels_base$traits, title = NULL, show_legend = TRUE)

p_final_base <- ggarrange(p_traits_base, p_avg_base, ncol = 1, heights = c(12.5, 1),
                          common.legend = TRUE, legend = "bottom")

p_final_base

ggsave(filename = file.path(output_dir, "Figure1_Heatmap_variance_explained_Baseline.pdf"),
       plot = p_final_base, width = 10, height = 11, units = "in", device = cairo_pdf)

ggsave(filename = file.path(output_dir, "Figure1_Heatmap_variance_explained_Baseline.png"),
       plot = p_final_base, width = 10, height = 11, units = "in", dpi = 300, bg = "white")


ggsave(filename = file.path(output_dir, "Figure1_Heatmap_variance_explained.pdf"),
       plot = p_traits_base, width = 10, height = 12, units = "in", device = cairo_pdf)

ggsave(filename = file.path(output_dir, "Figure1_Heatmap_variance_explained.png"),
       plot = p_traits_base, width = 10, height = 12, units = "in", dpi = 300, bg = "white")


## Simplified figure for presentation: Amount, Frequency, and A × F only
factor_levels_simple <- c("Amount (A)", "Frequency (F)", "A × F")

heat_simple <- heat_base %>%
  filter(Factor %in% factor_levels_simple)

panels_simple <- format_heatmap_panel(heat_simple, factor_levels = factor_levels_simple, include_average = TRUE)

p_traits_simple <- plot_trait_panel(panels_simple$traits, title = NULL, show_legend = TRUE)
p_traits_simple

ggsave(file.path(output_dir, "Heatmap_Amount_Frequency_Interaction_slide.png"), p_traits_simple,
       width = 7, height = 10, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Heatmap_Amount_Frequency_Interaction_slide.pdf"), p_traits_simple,
       width = 7, height = 10, units = "in", device = cairo_pdf)


## Full trait names for the simplified slide figure

trait_full_labels <- c(
  "N"        = "Plant number per pot",
  "H"        = "Plant height",
  "SB"       = "Shoot biomass",
  "TB"       = "Total biomass",
  "R/S"      = "Root-to-shoot ratio",
  "SLA"      = "Specific leaf area",
  "LDMC"     = "Leaf dry matter content",
  "LNC"      = "Leaf nitrogen content",
  "Leaf C:N" = "Leaf carbon-to-nitrogen ratio",
  "Chl"      = "Chlorophyll content",
  "PSIIeff"  = "Efficiency of Photosystem II",
  "Ψosm"     = "Osmotic potential",
  "SS"       = "Stomatal size",
  "SD"       = "Stomatal density",
  "RD"       = "Root diameter",
  "RDMC"     = "Root dry matter content",
  "SRL"      = "Specific root length",
  "SRA"      = "Specific root area",
  "RTD"      = "Root tissue density",
  "RNC"      = "Root nitrogen content",
  "Root C:N" = "Root carbon-to-nitrogen ratio"
)

p_traits_simple <- plot_trait_panel(
  panels_simple$traits,
  title = NULL,
  show_legend = TRUE
) +
  scale_y_discrete(labels = trait_full_labels)

p_traits_simple


ggsave(file.path(output_dir, "Heatmap_Amount_Frequency_Interaction_slide-Sept21.png"), p_traits_simple,
       width = 9, height = 10, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Heatmap_Amount_Frequency_Interaction_slide-Sept21.pdf"), p_traits_simple,
       width = 9, height = 10, units = "in", device = cairo_pdf)



## Figure S5: four alternative models
heat_bio  <- prepare_heatmap_data(data_bio,  factor_labels_bio,  include_average = FALSE)
heat_den  <- prepare_heatmap_data(data_den,  factor_labels_den,  include_average = FALSE)
heat_full <- prepare_heatmap_data(data_full, factor_labels_full, include_average = FALSE)
heat_sub  <- prepare_heatmap_data(data_sub,  factor_labels_sub,  include_average = FALSE)

panel_bio  <- format_heatmap_panel(heat_bio,  factor_levels_bio,  include_average = FALSE)
panel_den  <- format_heatmap_panel(heat_den,  factor_levels_den,  include_average = FALSE)
panel_full <- format_heatmap_panel(heat_full, factor_levels_full, include_average = FALSE)
panel_sub  <- format_heatmap_panel(heat_sub,  factor_levels_sub,  include_average = FALSE)

p_bio <- plot_trait_panel(panel_bio, title = "Biomass-adjusted (A × F × C + B)", show_legend = FALSE)

p_den <- plot_trait_panel(panel_den, title = "Density-adjusted (A × F × C + D)", show_legend = FALSE)

p_full <- plot_trait_panel(panel_full, title = "Fully adjusted (A × F × C + B + D)", show_legend = FALSE)

p_sub <- plot_trait_panel(panel_sub, title = "Biomass substitution (A × F × B)", show_legend = TRUE)

p1 <- ggarrange(p_den, p_sub, ncol = 1, common.legend = TRUE, legend = "bottom",
                labels = c("(b)", ("(d)")), heights = c(1, 0.3))
p1

p2 <- ggarrange(p_bio, p_full, ncol = 1, labels = c("(a)", ("(c)")))
p2


p_alt_models <- ggarrange(p2, p1, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

p_alt_models

ggsave(filename = file.path(output_dir, "FigureS5_Heatmap_variance_explained_alternative_models.pdf"),
       plot = p_alt_models, width = 22, height = 14, units = "in", device = cairo_pdf)

ggsave(filename = file.path(output_dir, "FigureS5_Heatmap_variance_explained_alternative_models.png"),
       plot = p_alt_models, width = 22, height = 14, units = "in", dpi = 300, bg = "white")
