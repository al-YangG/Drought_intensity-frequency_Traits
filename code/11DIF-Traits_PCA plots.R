############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: PCA plots for all traits and trait categories
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(grid)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Import data
dat <- read.csv(input_file, header = TRUE) %>%
  mutate(Intensity = factor(Intensity, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         Treatment = factor(paste(Intensity, Frequency, sep = "-"),
                            levels = c("Dry-High", "Dry-Low", "Wet-High", "Wet-Low")))

# ---------------------------------------------------------
# 1. Trait matrices
# ---------------------------------------------------------

traits_all <- dat %>%
  select(N, H, SB, TB, R_S, SLA, LDMC, LNC, Leaf_CN,
         Chl, PSIIeff, Ψosm, SS, SD, RD, RDMC, SRL, RTD, RNC, Root_CN) %>%
  rename("R/S" = R_S, "Leaf C:N" = Leaf_CN, "Root C:N" = Root_CN)

traits_growth <- dat %>%
  select(N, H, SB, TB, R_S) %>%
  rename("R/S" = R_S)

traits_leaf <- dat %>%
  select(SLA, LDMC, LNC, Leaf_CN) %>%
  rename("Leaf C:N" = Leaf_CN)

traits_physiological <- dat %>%
  select(Chl, PSIIeff, Ψosm, SS, SD)

traits_root <- dat %>%
  select(RD, RDMC, SRL, RTD, RNC, Root_CN) %>%
  rename("Root C:N" = Root_CN)

# ---------------------------------------------------------
# 2. Run PCA
# ---------------------------------------------------------

pca_all <- prcomp(traits_all, center = TRUE, scale. = TRUE)
pca_growth <- prcomp(traits_growth, center = TRUE, scale. = TRUE)
pca_leaf <- prcomp(traits_leaf, center = TRUE, scale. = TRUE)
pca_physiological <- prcomp(traits_physiological, center = TRUE, scale. = TRUE)
pca_root <- prcomp(traits_root, center = TRUE, scale. = TRUE)

# ---------------------------------------------------------
# 3. Helpers
# ---------------------------------------------------------

prepare_pca_plot_data <- function(pca_obj, data, category = NULL, axes = c(1, 2)) {
  stopifnot(length(axes) == 2)
  
  ind <- factoextra::get_pca_ind(pca_obj)
  var <- factoextra::get_pca_var(pca_obj)
  
  ind_df <- as.data.frame(ind$coord)[, paste0("Dim.", axes), drop = FALSE]
  names(ind_df) <- c("Dim.1", "Dim.2")
  
  ind_df <- ind_df %>%
    mutate(Treatment = data$Treatment, Intensity = data$Intensity,
           Frequency = data$Frequency, Composition = data$Composition,
           HullGrp = interaction(Intensity, Frequency, Composition, sep = " × ", drop = TRUE))
  
  var_df <- as.data.frame(var$coord)[, paste0("Dim.", axes), drop = FALSE]
  names(var_df) <- c("Dim.1", "Dim.2")
  var_df$Var <- rownames(var$coord)
  
  range_x_ind <- diff(range(ind_df$Dim.1, na.rm = TRUE))
  range_y_ind <- diff(range(ind_df$Dim.2, na.rm = TRUE))
  range_x_var <- diff(range(var_df$Dim.1, na.rm = TRUE))
  range_y_var <- diff(range(var_df$Dim.2, na.rm = TRUE))
  
  scale_factor <- min(range_x_ind / range_x_var, range_y_ind / range_y_var)
  
  var_df_scaled <- var_df
  var_df_scaled$Dim.1 <- var_df_scaled$Dim.1 * scale_factor
  var_df_scaled$Dim.2 <- var_df_scaled$Dim.2 * scale_factor
  
  eig <- factoextra::get_eigenvalue(pca_obj)
  pc_labels <- c(sprintf("PC%d %.2f%%", axes[1], eig[axes[1], "variance.percent"]),
                 sprintf("PC%d %.2f%%", axes[2], eig[axes[2], "variance.percent"]))
  
  hulls <- ind_df %>%
    group_by(HullGrp, Treatment, Composition) %>%
    filter(n() >= 3) %>%
    slice(chull(Dim.1, Dim.2)) %>%
    ungroup()
  
  list(ind_df = ind_df, var_df_scaled = var_df_scaled, hulls = hulls,
       pc_labels = pc_labels, category = category)
}

make_pca_plot <- function(plot_data, colors, shape_values, category_label = NULL,
                          show_legend = TRUE) {
  p <- ggplot() +
    geom_point(data = plot_data$ind_df, size = 3, alpha = 0.5,
               aes(Dim.1, Dim.2, fill = Treatment, shape = Composition)) +
    geom_segment(data = plot_data$var_df_scaled,
                 aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 linewidth = 0.6, color = "black") +
    geom_polygon(data = plot_data$hulls,
                 aes(Dim.1, Dim.2, group = HullGrp, color = Treatment, linetype = Composition),
                 fill = NA, linewidth = 0.5, alpha = 0.9, show.legend = FALSE) +
    ggrepel::geom_text_repel(data = plot_data$var_df_scaled, aes(Dim.1, Dim.2, label = Var),
                             size = 4, fontface = "bold") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    coord_equal() +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shape_values) +
    labs(x = plot_data$pc_labels[1], y = plot_data$pc_labels[2],
         fill = "Watering treatment", shape = "Plant composition") +
    theme_classic() +
    theme(title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10),
          axis.ticks = element_line(linewidth = 0.3),
          axis.title = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 11),
          legend.text = element_text(face = "bold", size = 11),
          legend.position = if (show_legend) "bottom" else "none",
          legend.box = "vertical",
          legend.direction = "horizontal",
          legend.box.just = "center",
          aspect.ratio = 0.85,
          axis.line = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  
  if (!is.null(category_label)) {
    p <- p +
      annotate("text", x = -Inf, y = Inf, label = category_label, hjust = -0.1, vjust = 1.5,
               fontface = "bold.italic", size = 4, color = "#087C08")
  }
  
  p +
    guides(fill = guide_legend(title = "Watering treatment", order = 1, byrow = TRUE, keyheight = unit(0.6, "lines"),
                               override.aes = list(shape = 21, color = "black", alpha = 0.5)),
           shape = guide_legend(title = "Plant composition", order = 2, byrow = TRUE, keyheight = unit(0.6, "lines"),
                                override.aes = list(fill = "white", color = "black", alpha = 1)))
}

# ---------------------------------------------------------
# 4. Styling
# ---------------------------------------------------------

treatment_colors <- c("Wet-Low" = "#1994E5", "Wet-High" = "#B38AF1", "Dry-Low" = "#E5D019", "Dry-High" = "#E7782F")

shape_comp <- c(Conspecific = 21, Heterospecific = 24)

# ---------------------------------------------------------
# 5. Overall PCA plot
# ---------------------------------------------------------

plot_all <- prepare_pca_plot_data(pca_all, dat, category = "All")

p_all <- make_pca_plot(plot_data = plot_all, colors = treatment_colors, shape_values = shape_comp,
                       category_label = NULL, show_legend = TRUE)
p_all

# Save - Figure 5
ggsave(filename = file.path(output_dir, "PCA_all_traits.pdf"), plot = p_all,
       width = 10, height = 9, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "PCA_all_traits.png"), plot = p_all,
       width = 10, height = 9, units = "in", dpi = 300)


# ---------------------------------------------------------
# 6. Category-specific PCA plots
# ---------------------------------------------------------

plot_growth <- prepare_pca_plot_data(pca_growth, dat, category = "Growth")
plot_leaf <- prepare_pca_plot_data(pca_leaf, dat, category = "Leaf")
plot_physiological <- prepare_pca_plot_data(pca_physiological, dat, category = "Physiological")
plot_root <- prepare_pca_plot_data(pca_root, dat, category = "Root")

p_growth <- make_pca_plot(plot_data = plot_growth, colors = treatment_colors,
                          shape_values = shape_comp, category_label = "Growth traits",
                          show_legend = TRUE)

p_leaf <- make_pca_plot(plot_data = plot_leaf, colors = treatment_colors,
                        shape_values = shape_comp, category_label = "Leaf traits",
                        show_legend = TRUE)

p_physiological <- make_pca_plot(plot_data = plot_physiological, colors = treatment_colors,
                                 shape_values = shape_comp, category_label = "Physiological traits",
                                 show_legend = TRUE)

p_root <- make_pca_plot(plot_data = plot_root, colors = treatment_colors,
                        shape_values = shape_comp, category_label = "Root traits",
                        show_legend = TRUE)

p_categories <- ggarrange(p_growth, p_leaf, p_physiological, p_root, ncol = 2, nrow = 2,
                          common.legend = TRUE, legend = "bottom",
                          labels = c("(a)", "(b)", "(c)", "(d)"))
p_categories

# Save - Figure S8
ggsave(filename = file.path(output_dir, "PCA_trait_categories.pdf"), plot = p_categories,
       width = 10, height = 10, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "PCA_trait_categories.png"), plot = p_categories,
       width = 10, height = 10, units = "in", dpi = 300, bg = "white")

