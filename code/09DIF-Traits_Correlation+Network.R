############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Trait correlations and correlation networks
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(tidyr)
library(Hmisc)
library(ggplot2)
library(ggnewscale)
library(tibble)
library(purrr)
library(tidygraph)
library(ggraph)
library(ggpubr)
library(igraph)
library(ggtext)
library(scales)
library(grid)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
figure_dir <- "results/figures"
table_dir <- "results/tables"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Import data
dat <- read.csv(input_file, header = TRUE)

# ---------------------------------------------------------
# 1. Trait definitions
# ---------------------------------------------------------

trait_order <- c("N", "H", "SB", "TB", "R/S",
                 "SLA", "LDMC", "LNC", "Leaf C:N",
                 "Chl", "PSIIeff", "Ψosm", "SS", "SD",
                 "RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root C:N")

trait_category <- tribble(~Trait,       ~Category,
                          "N",          "Growth",
                          "H",          "Growth",
                          "SB",         "Growth",
                          "TB",         "Growth",
                          "R/S",        "Growth",
                          "SLA",        "Leaf",
                          "LDMC",       "Leaf",
                          "LNC",        "Leaf",
                          "Leaf C:N",   "Leaf",
                          "Chl",        "Physiological",
                          "PSIIeff",    "Physiological",
                          "Ψosm",       "Physiological",
                          "SS",         "Physiological",
                          "SD",         "Physiological",
                          "RD",         "Root",
                          "RDMC",       "Root",
                          "SRL",        "Root",
                          "SRA",        "Root",
                          "RTD",        "Root",
                          "RNC",        "Root",
                          "Root C:N",   "Root")

category_colors <- c(Growth = "#E5E5E5", Leaf = "#C7E9C0",
                     Physiological = "#D9EDF7", Root = "#FCF8E3")

block_sizes <- c(Growth = 5, Leaf = 4, Physiological = 5, Root = 7)

# Helper: select and rename traits
get_trait_data <- function(data) {
  data %>%
    select(N, H, SB, TB, R_S, SLA, LDMC, LNC, Leaf_CN,
           Chl, PSIIeff, Ψosm, SS, SD, RD, RDMC, SRL, SRA, RTD, RNC, Root_CN) %>%
    rename("R/S" = R_S, "Leaf C:N" = Leaf_CN, "Root C:N" = Root_CN)
}

# Helper: z-score standardization
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# ---------------------------------------------------------
# 2. Correlation matrix
# ---------------------------------------------------------

traits_all <- get_trait_data(dat)

traits_std <- traits_all %>%
  mutate(across(everything(), z_score))

# Pearson correlations and p-values
cor_res <- rcorr(as.matrix(traits_std), type = "pearson")

flatten_corr_matrix <- function(cor_mat, p_mat) {
  upper_idx <- upper.tri(cor_mat)
  tibble(row = rownames(cor_mat)[row(cor_mat)[upper_idx]],
         column = colnames(cor_mat)[col(cor_mat)[upper_idx]],
         cor = cor_mat[upper_idx], p = p_mat[upper_idx])
}

cor_table <- flatten_corr_matrix(cor_res$r, cor_res$P)

high_cor <- cor_table %>%
  filter(abs(cor) > 0.7, p < 0.05) %>%
  arrange(desc(abs(cor)))

write.csv(high_cor, file.path(table_dir, "Trait_high_correlations.csv"), row.names = FALSE)

# Correlation matrix for plotting
traits_complete <- traits_all %>%
  drop_na()

cor_mat <- cor(traits_complete, method = "pearson")
p_mat <- cor_res$P

cor_mat <- cor_mat[trait_order, trait_order]
p_mat <- p_mat[trait_order, trait_order]

cor_long <- as.data.frame(as.table(cor_mat)) %>%
  rename(Trait1 = Var1, Trait2 = Var2, r = Freq) %>%
  mutate(p = as.vector(p_mat), i = match(Trait1, trait_order), j = match(Trait2, trait_order)) %>%
  filter(i < j) %>%
  mutate(r_plot = ifelse(p < 0.05, r, NA_real_), x = j, y = length(trait_order) - i + 1)

n_traits <- length(trait_order)

top_strip <- trait_category %>%
  mutate(j = match(Trait, trait_order), x = j, y = n_traits + 1)

left_strip <- trait_category %>%
  mutate(i = match(Trait, trait_order), x = 0, y = n_traits - i + 1)

cuts <- cumsum(block_sizes)
cuts_x <- cuts + 0.5
cuts_x <- cuts_x[cuts_x > 0.5 & cuts_x < (n_traits + 0.5)]

cuts_y <- (n_traits - cuts) + 0.5
cuts_y <- cuts_y[cuts_y > 0.5 & cuts_y < (n_traits + 0.5)]

p_catcorr <- ggplot() +
  geom_tile(data = cor_long, aes(x = x, y = y, fill = r_plot), color = "white") +
  geom_text(data = cor_long, aes(x = x, y = y, label = ifelse(!is.na(r_plot), sprintf("%.2f", r), "")), size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1), na.value = "white", name = "Pearson r") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = top_strip, aes(x = x, y = y, fill = Category), height = 0.9) +
  geom_tile(data = left_strip, aes(x = x, y = y, fill = Category), width = 0.9) +
  scale_fill_manual(values = category_colors, guide = "none") +
  geom_segment(data = data.frame(x = cuts_x),
               aes(x = x, xend = x, y = 0.5, yend = n_traits + 0.5),
               linetype = "dotted", linewidth = 0.4, inherit.aes = FALSE) +
  geom_segment(data = data.frame(y = cuts_y),
               aes(x = 0.5, xend = n_traits + 0.5, y = y, yend = y),
               linetype = "dotted", linewidth = 0.4, inherit.aes = FALSE) +
  scale_x_continuous(breaks = 1:n_traits, labels = trait_order, position = "top",
                     expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_continuous(breaks = 1:n_traits, labels = rev(trait_order),
                     expand = expansion(mult = c(0.12, 0.02))) +
  coord_fixed(clip = "off") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.12, 0.15),
        legend.justification = c(0, 0),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey40", linewidth = 0.4))

starts <- c(1, head(cumsum(block_sizes) + 1, -1))
ends <- cumsum(block_sizes)
centers <- (starts + ends) / 2

category_labels <- data.frame(Category = names(block_sizes), x = centers, y = centers)

p_catcorr <- p_catcorr +
  annotate("text", x = category_labels$x, y = n_traits + 1,
           label = category_labels$Category, fontface = "bold", size = 3.5) +
  annotate("text", x = 0, y = n_traits - category_labels$y + 1,
           label = category_labels$Category, fontface = "bold", size = 3.5, angle = 90)
p_catcorr

# Save - Figure S7
ggsave(filename = file.path(figure_dir, "Trait_correlation_matrix.pdf"), plot = p_catcorr,
       width = 10, height = 10, units = "in", device = cairo_pdf)
ggsave(filename = file.path(figure_dir, "Trait_correlation_matrix.png"), plot = p_catcorr,
       width = 10, height = 10, units = "in", dpi = 300)


# ---------------------------------------------------------
# 3. Trait correlation networks
# ---------------------------------------------------------

build_trait_network <- function(data_subset, cor_thresh = 0.5, p_thresh = 0.01, seed = 42) {
  traits <- get_trait_data(data_subset)
  
  traits_std <- traits %>%
    mutate(across(everything(), z_score))
  
  cor_res <- rcorr(as.matrix(traits_std), type = "pearson")
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  
  sig_mat <- (abs(cor_mat) > cor_thresh) & (p_mat < p_thresh)
  
  edges <- which(sig_mat, arr.ind = TRUE)
  
  edges_df <- tibble(from = rownames(cor_mat)[edges[, 1]],
                     to = colnames(cor_mat)[edges[, 2]],
                     correlation = cor_mat[edges]) %>%
    filter(from < to) %>%
    mutate(weight = abs(correlation))
  
  graph_tbl <- tbl_graph(nodes = tibble(name = colnames(cor_mat)),
                         edges = edges_df, directed = FALSE) %>%
    activate(nodes) %>%
    left_join(trait_category, by = c("name" = "Trait")) %>%
    mutate(Category = factor(Category, levels = c("Growth", "Leaf", "Physiological", "Root")),
           degree = centrality_degree(),
           betweenness = centrality_betweenness(weights = 1 / weight),
           module = as.factor(group_louvain(weights = weight)))
  
  set.seed(seed)
  
  plot_obj <- ggraph(graph_tbl, layout = "circle") +
    geom_edge_arc(aes(width = weight, color = ifelse(correlation > 0, "Positive", "Negative")), strength = 0.3) +
    geom_node_point(aes(size = degree, fill = Category), shape = 21, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold",
                   point.padding = unit(0.4, "lines"), box.padding = unit(0.4, "lines")) +
    scale_edge_width(range = c(0.2, 2), name = "Correlation strength") +
    scale_edge_color_manual(values = c(Positive = "#E31A1C", Negative = "#2C7FB8"),
                            name = "Correlation sign") +
    scale_size_continuous(range = c(2, 8), breaks = c(2, 4, 6), name = "Degree centrality") +
    scale_fill_manual(values = category_colors, name = "Trait category") +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 4)),
           edge_color = guide_legend(order = 2, override.aes = list(edge_width = 1)),
           edge_width = guide_legend(order = 3), size = guide_legend(order = 4)) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "black"),
          plot.subtitle = element_text(face = "bold", size = 10, hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 10))
  
  list(graph = graph_tbl, plot = plot_obj)
}

# Subsets
subsets <- list(Wet = dat %>% filter(Intensity == "Wet"),
                Dry = dat %>% filter(Intensity == "Dry"),
                Low = dat %>% filter(Frequency == "Low"),
                High = dat %>% filter(Frequency == "High"),
                Conspecific = dat %>% filter(Composition == "Conspecific"),
                Heterospecific = dat %>% filter(Composition == "Heterospecific"),
                All = dat)

networks <- list(Wet = build_trait_network(subsets$Wet),
                 Dry = build_trait_network(subsets$Dry),
                 Low = build_trait_network(subsets$Low),
                 High = build_trait_network(subsets$High),
                 Conspecific = build_trait_network(subsets$Conspecific),
                 Heterospecific = build_trait_network(subsets$Heterospecific),
                 All = build_trait_network(subsets$All))

# Panel colors
col_intensity <- c(Wet = "#1994E5", Dry = "#E7782F")
col_frequency <- c(Low = "#E5D019", High = "#B38AF1")
col_composition <- c(Conspecific = "#D95F8D", Heterospecific = "#33A27F")

format_network_plot <- function(plot_obj, title, subtitle, fill_color) {
  plot_obj +
    labs(title = title, subtitle = subtitle) +
    theme(plot.title = ggtext::element_textbox_simple(fill = alpha(fill_color, 0.2), color = "black",
                                                      box.color = NA, halign = 0.5,
                                                      padding = margin(4, 6, 4, 6),
                                                      size = 12, face = "bold"))
}

p_wet <- format_network_plot(networks$Wet$plot, "Watering intensity – Wet",
                             paste0("(n = ", nrow(subsets$Wet), ")"), col_intensity["Wet"])

p_dry <- format_network_plot(networks$Dry$plot, "Watering intensity – Dry",
                             paste0("(n = ", nrow(subsets$Dry), ")"), col_intensity["Dry"])

p_low <- format_network_plot(networks$Low$plot, "Watering frequency – Low",
                             paste0("(n = ", nrow(subsets$Low), ")"), col_frequency["Low"])

p_high <- format_network_plot(networks$High$plot, "Watering frequency – High",
                              paste0("(n = ", nrow(subsets$High), ")"), col_frequency["High"])

p_conspecific <- format_network_plot(networks$Conspecific$plot, "Plant composition – Conspecific",
                                     paste0("(n = ", nrow(subsets$Conspecific), ")"),
                                     col_composition["Conspecific"])

p_heterospecific <- format_network_plot(networks$Heterospecific$plot, "Plant composition – Heterospecific",
                                        paste0("(n = ", nrow(subsets$Heterospecific), ")"),
                                        col_composition["Heterospecific"])

legend_net <- get_legend(p_wet)

net_i <- ggarrange(p_wet + theme(legend.position = "none"), NULL,
                   p_dry + theme(legend.position = "none"), ncol = 3,
                   widths = c(4, 1, 4), labels = c("(a)", "", "(b)"))

net_f <- ggarrange(p_low + theme(legend.position = "none"), legend_net,
                   p_high + theme(legend.position = "none"), ncol = 3,
                   widths = c(4, 1, 4), labels = c("(c)", "", "(d)"))

net_c <- ggarrange(p_conspecific + theme(legend.position = "none"), NULL,
                   p_heterospecific + theme(legend.position = "none"),
                   ncol = 3, widths = c(4, 1, 4), labels = c("(e)", "", "(f)"))

p_net <- ggarrange(net_i, NULL, net_f, NULL, net_c, ncol = 1,
                   heights = c(6, 1, 6, 1, 6))
p_net

# Save - Figure 4
ggsave(filename = file.path(figure_dir, "Trait_networks.pdf"), plot = p_net,
       width = 13, height = 17, units = "in", device = cairo_pdf)
ggsave(filename = file.path(figure_dir, "Trait_networks.png"), plot = p_net,
       width = 13, height = 17, units = "in", dpi = 300, bg = "white")


# ---------------------------------------------------------
# 4. Network metrics
# ---------------------------------------------------------

compute_network_metrics <- function(graph_tbl) {
  g <- igraph::as.igraph(graph_tbl)
  
  clustering <- igraph::cluster_louvain(g)
  membership <- igraph::membership(clustering)
  
  tibble(Edge_density = igraph::edge_density(g),
         Avg_path_length = tryCatch(igraph::mean_distance(g, directed = FALSE), error = function(e) NA_real_),
         Diameter = tryCatch(igraph::diameter(g), error = function(e) NA_real_),
         Avg_clust_coeff = igraph::transitivity(g, type = "average"),
         Modularity = igraph::modularity(clustering),
         Num_modules = length(unique(membership)),
         Num_trait_relationships = igraph::gsize(g))
}

network_metrics <- map_dfr(list(Dry = networks$Dry$graph,
                                Wet = networks$Wet$graph,
                                Low = networks$Low$graph,
                                High = networks$High$graph,
                                Conspecific = networks$Conspecific$graph,
                                Heterospecific = networks$Heterospecific$graph,
                                All = networks$All$graph),
                           compute_network_metrics, .id = "Treatment")

# Save - Table S11
write.csv(network_metrics, file.path(table_dir, "Trait_network_metrics.csv"), row.names = FALSE)
