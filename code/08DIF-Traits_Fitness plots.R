############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Fitness trait plots
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(car)
library(emmeans)
library(glmmTMB)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(ggrepel)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 1. Import and prepare data
# ---------------------------------------------------------

dat <- read.csv(input_file, header = TRUE) %>%
  mutate(Intensity = factor(Intensity, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")))

options(contrasts = c("contr.sum", "contr.poly"))

# Move proportions into (0, 1) for beta regression
squeeze01 <- function(x) {
  n <- sum(!is.na(x))
  (x * (n - 1) + 0.5) / n
}

dat <- dat %>%
  mutate(SER_beta = ifelse(SER <= 0 | SER >= 1, squeeze01(SER), SER),
         SSR1_beta = ifelse(SSR1 <= 0 | SSR1 >= 1, squeeze01(SSR1), SSR1),
         SSR2_beta = ifelse(SSR2 <= 0 | SSR2 >= 1, squeeze01(SSR2), SSR2))

# Helper for emmeans output
fix_emm_names <- function(df) {
  if ("asymp.LCL" %in% names(df)) names(df)[names(df) == "asymp.LCL"] <- "lower.CL"
  if ("asymp.UCL" %in% names(df)) names(df)[names(df) == "asymp.UCL"] <- "upper.CL"
  if ("emmean" %in% names(df)) names(df)[names(df) == "emmean"] <- "response"
  df
}

# P-value helpers
p_to_text <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
}

p_to_star <- function(p) {
  case_when(is.na(p) ~ NA_character_, p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", TRUE ~ "ns")
}

# Colors
col_intensity <- c(Wet = "#1994E5", Dry = "#E7782F")
col_frequency <- c(Low = "#E5D019", High = "#B38AF1")
col_composition <- c(Conspecific = "#D95F8D", Heterospecific = "#33A27F")
col_treatment <- c("#E7782F", "#E5D019", "#1994E5", "#B38AF1")

# ---------------------------------------------------------
# 2. Fit fitness models
# ---------------------------------------------------------

mods <- list(SER = glmmTMB(SER_beta ~ Intensity * Frequency * Composition,
                           data = dat, family = beta_family(link = "logit")),
             SSR1 = glmmTMB(SSR1_beta ~ Intensity * Frequency * Composition,
                            data = dat, family = beta_family(link = "logit")),
             SSR2 = glmmTMB(SSR2_beta ~ Intensity * Frequency * Composition,
                            data = dat, family = beta_family(link = "logit")),
             RP1 = lm(RP1 ~ Intensity * Frequency * Composition, data = dat),
             RP2 = lm(RP2 ~ Intensity * Frequency * Composition, data = dat))

# ---------------------------------------------------------
# 3. Type III tests and marginal means
# ---------------------------------------------------------

extract_type3 <- function(model, metric) {
  anova_table <- car::Anova(model, type = 3)
  p_col <- intersect(c("Pr(>Chisq)", "Pr(>F)"), colnames(anova_table))[1]
  
  tibble(Metric = metric, Effect = rownames(anova_table),
         p_value = unname(anova_table[, p_col]),
         p_label = p_to_text(p_value),
         p_star = p_to_star(p_value))
}

type3_tbl <- imap_dfr(mods, extract_type3)

effects_to_extract <- c("Intensity", "Frequency", "Composition",
                        "Intensity:Frequency", "Intensity:Composition", "Frequency:Composition")

effect_to_formula <- function(effect) {
  if (grepl(":", effect)) {
    vars <- strsplit(effect, ":", fixed = TRUE)[[1]]
    as.formula(paste0("~ ", vars[1], " * ", vars[2]))
  } else {
    as.formula(paste0("~ ", effect))
  }
}

extract_emm <- function(model, metric, effect) {
  emm_obj <- emmeans(model, effect_to_formula(effect), type = "response")
  
  emm_df <- summary(emm_obj, infer = TRUE) %>%
    as.data.frame() %>%
    fix_emm_names() %>%
    mutate(Metric = metric, Effect = effect)
  
  if (!grepl(":", effect)) {
    ct_df <- pairs(emm_obj, adjust = "sidak") %>%
      summary(infer = TRUE) %>%
      as.data.frame() %>%
      mutate(Metric = metric, Effect = effect, ContrastType = "pairwise")
  } else {
    vars <- strsplit(effect, ":", fixed = TRUE)[[1]]
    ct_df <- emmeans(model, as.formula(paste0("~ ", vars[1], " | ", vars[2])), type = "response") %>%
      pairs(adjust = "sidak") %>%
      summary(infer = TRUE) %>%
      as.data.frame() %>%
      mutate(Metric = metric, Effect = effect, ContrastType = paste0(vars[1], " | ", vars[2]))
  }
  
  ct_df <- ct_df %>%
    mutate(p_label = p_to_text(p.value), p_star = p_to_star(p.value))
  
  list(emm = emm_df, contrast = ct_df)
}

res_list <- imap(mods, \(m, metric) {
  map(effects_to_extract, \(eff) extract_emm(m, metric, eff))
})

emm_tbl <- map_dfr(res_list, \(x) map_dfr(x, "emm"))
ct_tbl <- map_dfr(res_list, \(x) map_dfr(x, "contrast"))

sig_effects <- type3_tbl %>%
  filter(Effect != "(Intercept)", !is.na(p_value), p_value < 0.05) %>%
  select(Metric, Effect)

emm_sig <- emm_tbl %>%
  semi_join(sig_effects, by = c("Metric", "Effect"))

# ---------------------------------------------------------
# 4. Plot significant main effects
# ---------------------------------------------------------

main_df <- emm_sig %>%
  filter(!grepl(":", Effect)) %>%
  mutate(X = case_when(Effect == "Intensity" ~ as.character(Intensity),
                       Effect == "Frequency" ~ as.character(Frequency),
                       Effect == "Composition" ~ as.character(Composition),
                       TRUE ~ NA_character_),
    X = factor(X, levels = c("Wet", "Dry", "Low", "High", "Conspecific", "Heterospecific")),
    MetricOrd = factor(Metric, levels = c("SER", "SSR1", "SSR2", "RP1", "RP2")),
    MetricMath = case_when(Metric == "SER" ~ "SER", Metric == "SSR1" ~ "SSR[1]",
                           Metric == "SSR2" ~ "SSR[2]", Metric == "RP1" ~ "RP[1]",
                           Metric == "RP2" ~ "RP[2]", TRUE ~ Metric),
    PanelMath = case_when(Effect == "Intensity" ~ paste0("bold(", MetricMath, " ~ '-' ~ 'Intensity')"),
                          Effect == "Frequency" ~ paste0("bold(", MetricMath, " ~ '-' ~ 'Frequency')"),
                          Effect == "Composition" ~ paste0("bold(", MetricMath, " ~ '-' ~ 'Composition')"))) %>%
  left_join(type3_tbl %>% select(Metric, Effect, p_value), by = c("Metric", "Effect")) %>%
  mutate(sig_lab = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                             p_value < 0.05 ~ "*", TRUE ~ ""))

panel_levels <- main_df %>%
  distinct(PanelMath, MetricOrd, Effect) %>%
  arrange(MetricOrd, Effect) %>%
  pull(PanelMath)

main_df <- main_df %>%
  mutate(PanelMath = factor(PanelMath, levels = panel_levels))

stars_df <- main_df %>%
  group_by(PanelMath) %>%
  summarise(sig_lab = first(sig_lab), x_mid = levels(droplevels(X))[2],
            y_pos = min(lower.CL, na.rm = TRUE) - 0.06 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
            .groups = "drop")

p_main <- ggplot(main_df, aes(x = X, y = response)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, color = X), linewidth = 0.6) +
  geom_text(data = stars_df, aes(x = x_mid, y = y_pos, label = sig_lab),
            inherit.aes = FALSE, fontface = "bold", size = 4, nudge_x = -0.5) +
  scale_color_manual( values = c(col_intensity, col_frequency, col_composition), guide = "none") +
  facet_wrap(~ PanelMath, scales = "free", labeller = label_parsed, ncol = 3) +
  labs(x = NULL, y = "Model-predicted value") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "gray95", color = "gray95"),
        axis.text.x = element_text(face = "bold", size = 11, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", size = 11, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank()) +
  coord_cartesian(clip = "off")
p_main

# ---------------------------------------------------------
# 5. Plot significant interactions
# ---------------------------------------------------------
# Helper: extract EMMs, simple-effect stars, and interaction p-labels
build_interaction_panel <- function(model, panel_id, emm_formula, simple_formula,
                                    x_var, color_var, x_levels, color_levels,
                                    anova_effect, test_col) {
  emm_df <- summary(emmeans(model, emm_formula, type = "response"), infer = TRUE) %>%
    as.data.frame() %>%
    fix_emm_names() %>%
    rename(response = any_of("response")) %>%
    mutate(x = factor(.data[[x_var]], levels = x_levels),
           group = factor(.data[[color_var]], levels = color_levels)) %>%
    transmute(panel = panel_id, x, response, lower.CL, upper.CL, group)
  
  ct_df <- pairs(emmeans(model, simple_formula, type = "response")) %>%
    summary() %>%
    as.data.frame() %>%
    mutate(x = factor(.data[[x_var]], levels = x_levels), star = p_to_star(p.value)) %>%
    select(x, star)
  
  star_df <- emm_df %>%
    group_by(panel, x) %>%
    summarise(y = max(upper.CL, na.rm = TRUE) + 0.06 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
              .groups = "drop") %>%
    left_join(ct_df, by = "x")
  
  interaction_p <- Anova(model, type = 3)[anova_effect, test_col]
  
  short_effect <- dplyr::recode(anova_effect,
                                "Intensity:Composition" = "I × C",
                                "Intensity:Frequency" = "I × F",
                                "Frequency:Composition" = "F × C",
                                .default = gsub(":", " × ", anova_effect))
  
  interaction_lab <- paste0(short_effect, ": p ", p_to_text(interaction_p))
  
  lab_df <- emm_df %>%
    group_by(panel) %>%
    summarise(x = 1.5, y = max(upper.CL, na.rm = TRUE) + 0.18 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(label = interaction_lab)
  
  list(emm = emm_df, stars = star_df, label = lab_df)
}

# Build all four panels
panel_ser_ic <- build_interaction_panel(model = mods$SER,
                                        panel_id = "bold(SER ~ '-' ~ Intensity %*% Composition)",
                                        emm_formula = ~ Intensity * Composition,
                                        simple_formula = ~ Intensity | Composition,
                                        x_var = "Composition", color_var = "Intensity",
                                        x_levels = c("Conspecific", "Heterospecific"), color_levels = c("Wet", "Dry"),
                                        anova_effect = "Intensity:Composition", test_col = "Pr(>Chisq)")

panel_rp1_if <- build_interaction_panel(model = mods$RP1,
                                        panel_id = "bold(RP[1] ~ '-' ~ Intensity %*% Frequency)",
                                        emm_formula = ~ Intensity * Frequency,
                                        simple_formula = ~ Intensity | Frequency,
                                        x_var = "Frequency", color_var = "Intensity",
                                        x_levels = c("Low", "High"), color_levels = c("Wet", "Dry"),
                                        anova_effect = "Intensity:Frequency", test_col = "Pr(>F)")

panel_rp1_ic <- build_interaction_panel(model = mods$RP1,
                                        panel_id = "bold(RP[1] ~ '-' ~ Intensity %*% Composition)",
                                        emm_formula = ~ Intensity * Composition,
                                        simple_formula = ~ Intensity | Composition,
                                        x_var = "Composition", color_var = "Intensity",
                                        x_levels = c("Conspecific", "Heterospecific"), color_levels = c("Wet", "Dry"),
                                        anova_effect = "Intensity:Composition", test_col = "Pr(>F)")

panel_rp2_fc <- build_interaction_panel(model = mods$RP2,
                                        panel_id = "bold(RP[2] ~ '-' ~ Frequency %*% Composition)",
                                        emm_formula = ~ Frequency * Composition,
                                        simple_formula = ~ Composition | Frequency,
                                        x_var = "Frequency", color_var = "Composition",
                                        x_levels = c("Low", "High"), color_levels = c("Conspecific", "Heterospecific"),
                                        anova_effect = "Frequency:Composition", test_col = "Pr(>F)")


# Combine data for one faceted plot
emm_int <- bind_rows(panel_ser_ic$emm, panel_rp1_if$emm, panel_rp1_ic$emm, panel_rp2_fc$emm) %>%
  mutate(panel = factor(panel, levels = c("bold(SER ~ '-' ~ Intensity %*% Composition)",
                                          "bold(RP[1] ~ '-' ~ Intensity %*% Frequency)",
                                          "bold(RP[1] ~ '-' ~ Intensity %*% Composition)",
                                          "bold(RP[2] ~ '-' ~ Frequency %*% Composition)")),
         x = factor(as.character(x), levels = c("Conspecific", "Heterospecific", "Low", "High")),
         group = factor(as.character(group), levels = c("Wet", "Dry", "Conspecific", "Heterospecific")))

stars_int <- bind_rows(panel_ser_ic$stars %>% mutate(panel = "bold(SER ~ '-' ~ Intensity %*% Composition)"),
                       panel_rp1_if$stars %>% mutate(panel = "bold(RP[1] ~ '-' ~ Intensity %*% Frequency)"),
                       panel_rp1_ic$stars %>% mutate(panel = "bold(RP[1] ~ '-' ~ Intensity %*% Composition)"),
                       panel_rp2_fc$stars %>% mutate(panel = "bold(RP[2] ~ '-' ~ Frequency %*% Composition)")) %>%
  mutate(panel = factor(panel, levels = levels(emm_int$panel)),
         x = factor(as.character(x), levels = levels(emm_int$x)))

labs_int <- bind_rows(panel_ser_ic$label %>% mutate(panel = "bold(SER ~ '-' ~ Intensity %*% Composition)"),
                      panel_rp1_if$label %>% mutate(panel = "bold(RP[1] ~ '-' ~ Intensity %*% Frequency)"),
                      panel_rp1_ic$label %>% mutate(panel = "bold(RP[1] ~ '-' ~ Intensity %*% Composition)"),
                      panel_rp2_fc$label %>% mutate(panel = "bold(RP[2] ~ '-' ~ Frequency %*% Composition)")) %>%
  mutate(panel = factor(panel, levels = levels(emm_int$panel)))

# One shared color palette and one shared legend
pal_interaction <- c(Wet = "#1994E5", Dry = "#E7782F",
                     Conspecific = "#D95F8D", Heterospecific = "#33A27F")

theme_interaction <- theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "gray95", color = "gray95"),
        axis.text.x = element_text(face = "bold", size = 11, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", size = 11, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 11),
        axis.line = element_blank())

p_interactions <- ggplot(emm_int, aes(x = x, y = response, color = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.6) +
  geom_text(data = stars_int, aes(x = x, y = y, label = star), inherit.aes = FALSE,
            fontface = "bold", size = 4, vjust = 4.5) +
  geom_text(data = labs_int, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            fontface = "bold", size = 4, hjust = 1) +
  facet_wrap(~ panel, ncol = 2, scales = "free", labeller = label_parsed) +
  scale_color_manual(values = pal_interaction, breaks = c("Wet", "Dry", "Conspecific", "Heterospecific")) +
  labs(x = NULL, y = "Model-predicted value", color = NULL) +
  coord_cartesian(clip = "off") +
  theme_interaction
p_interactions

# ---------------------------------------------------------
# 6. PCA of fitness traits
# ---------------------------------------------------------
dat <- read.csv(input_file, header = TRUE) %>%
  mutate(Intensity = factor(Intensity, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         Treatment = factor(paste(Intensity, Frequency, sep = "-"),
                            levels = c("Dry-High", "Dry-Low", "Wet-High", "Wet-Low")))
fitness_data <- dat %>%
  select(SER, SSR1, SSR2, SPtotal, RP1, RP2)

pca_result <- prcomp(fitness_data, center = TRUE, scale. = TRUE)

prepare_pca_plot_data <- function(pca_obj, groups = NULL, axes = c(1, 2)) {
  ind <- factoextra::get_pca_ind(pca_obj)
  var <- factoextra::get_pca_var(pca_obj)
  
  ind_df <- as.data.frame(ind$coord)[, paste0("Dim.", axes), drop = FALSE]
  names(ind_df) <- c("Dim.1", "Dim.2")
  
  if (!is.null(groups)) {
    ind_df$Group <- factor(groups)
  }
  
  var_df <- as.data.frame(var$coord)[, paste0("Dim.", axes), drop = FALSE]
  names(var_df) <- c("Dim.1", "Dim.2")
  var_df$Var <- rownames(var$coord)
  
  scale_factor <- min(diff(range(ind_df$Dim.1)) / diff(range(var_df$Dim.1)),
                      diff(range(ind_df$Dim.2)) / diff(range(var_df$Dim.2)))
  
  var_df$Dim.1 <- var_df$Dim.1 * scale_factor
  var_df$Dim.2 <- var_df$Dim.2 * scale_factor
  
  eig <- factoextra::get_eigenvalue(pca_obj)
  pc_labels <- c(sprintf("PC%d %.2f%%", axes[1], eig[axes[1], "variance.percent"]),
                 sprintf("PC%d %.2f%%", axes[2], eig[axes[2], "variance.percent"]))
  
  list(ind_df = ind_df, var_df = var_df, pc_labels = pc_labels)
}

plot_data <- prepare_pca_plot_data(pca_result, groups = dat$Treatment)

plot_data$var_df$Label <- case_when(
  plot_data$var_df$Var == "SER" ~ "bold(SER)",
  plot_data$var_df$Var == "SSR1" ~ "bold(SSR['1'])",
  plot_data$var_df$Var == "SSR2" ~ "bold(SSR['2'])",
  plot_data$var_df$Var == "RP1" ~ "bold(RP['1'])",
  plot_data$var_df$Var == "RP2" ~ "bold(RP['2'])",
  plot_data$var_df$Var == "SPtotal" ~ "bold(SP[total])",
  TRUE ~ plot_data$var_df$Var)

plot_data$ind_df <- plot_data$ind_df %>%
  mutate(Intensity = factor(dat$Intensity, levels = c("Wet", "Dry")),
         Frequency = factor(dat$Frequency, levels = c("Low", "High")),
         Composition = factor(dat$Composition, levels = c("Conspecific", "Heterospecific")),
         Group = factor(dat$Treatment),
         HullGrp = interaction(Intensity, Frequency, Composition, sep = " × ", drop = TRUE))

shape_comp <- c(Conspecific = 21, Heterospecific = 24)

hulls <- plot_data$ind_df %>%
  group_by(HullGrp, Group, Composition) %>%
  filter(n() >= 3) %>%
  slice(chull(Dim.1, Dim.2)) %>%
  ungroup()

theme_pca <- theme_classic() +
  theme(axis.text = element_text(face = "bold", color = "black", size = 10),
        axis.title = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(face = "bold", size = 11),
        legend.position = "right",
        aspect.ratio = 0.85,
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

p_pca <- ggplot() +
  geom_point(data = plot_data$ind_df, aes(Dim.1, Dim.2, fill = Group, shape = Composition),
             size = 3, alpha = 0.5) +
  geom_segment(data = plot_data$var_df, aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.6, color = "black") +
  geom_polygon(data = hulls, aes(Dim.1, Dim.2, group = HullGrp, color = Group, linetype = Composition),
               fill = NA, linewidth = 0.5, alpha = 0.9, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = plot_data$var_df, aes(Dim.1, Dim.2, label = Label),
                           parse = TRUE, size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  coord_equal() +
  scale_fill_manual(values = col_treatment) +
  scale_color_manual(values = col_treatment) +
  scale_shape_manual(values = shape_comp) +
  labs(x = plot_data$pc_labels[1], y = plot_data$pc_labels[2],
       fill = "Watering treatment", shape = "Plant composition") +
  theme_pca +
  guides(fill = guide_legend(title = "Watering treatment", order = 1, byrow = TRUE, keyheight = unit(1.5, "lines"),
                             override.aes = list(shape = 21, color = "black", alpha = 0.5)),
         shape = guide_legend(title = "Plant composition", order = 2, byrow = TRUE, keyheight = unit(1.5, "lines"),
                              override.aes = list(fill = "white", color = "black", alpha = 1)))
p_pca

# ---------------------------------------------------------
# 7. Combine and save
# ---------------------------------------------------------
p_right <- ggarrange(p_interactions, p_pca, labels = c("(b)", "(c)"),
                     heights = c(3, 2), nrow = 2)

p_final <- ggarrange(p_main, p_right, labels = c("(a)", " "),
                     widths = c(1.2, 1), nrow = 1)
p_final

# Save - Figure 3
ggsave(filename = file.path(output_dir, "Fitness_all.pdf"), plot = p_final,
       width = 16, height = 10, units = "in", device = cairo_pdf)
ggsave(filename = file.path(output_dir, "Fitness_all.png"), plot = p_final,
       width = 16, height = 10, units = "in", dpi = 300, bg = "white")

