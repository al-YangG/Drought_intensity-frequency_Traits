############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: PCA and RDA of functional traits
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(tibble)
library(purrr)
library(vegan)
library(factoextra)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/tables"

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

# Functional traits used in overall PCA/RDA
# SRA and RB excluded to reduce redundancy
traits_all <- dat %>%
  select(N, H, SB, TB, R_S, SLA, LDMC, LNC, Leaf_CN,
         Chl, PSIIeff, Ψosm, SS, SD, RD, RDMC, SRL, RTD, RNC, Root_CN) %>%
  rename("R/S" = R_S, "Leaf C:N" = Leaf_CN, "Root C:N" = Root_CN, "Phi_osm" = Ψosm)

# Trait subsets
traits_growth <- dat %>%
  select(N, H, SB, TB, R_S) %>%
  rename("R/S" = R_S)

traits_leaf <- dat %>%
  select(SLA, LDMC, LNC, Leaf_CN) %>%
  rename("Leaf C:N" = Leaf_CN)

traits_physiological <- dat %>%
  select(Chl, PSIIeff, Ψosm, SS, SD) %>%
  rename("Phi_osm" = Ψosm)

traits_root <- dat %>%
  select(RD, RDMC, SRL, RTD, RNC, Root_CN) %>%
  rename("Root C:N" = Root_CN)

traits_fitness <- dat %>%
  select(SER, SSR1, SSR2, SPtotal, RP1, RP2)

# Explanatory variables for RDA
design_matrix <- dat %>%
  select(Intensity, Frequency, Composition) %>%
  droplevels()

# ---------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------

# Extract PCA variable loadings, contributions, and cos2
pca_variable_table <- function(pca_obj, n_pc = 2) {
  var_info <- get_pca_var(pca_obj)
  
  out <- tibble(Trait = rownames(var_info$coord))
  
  for (i in seq_len(n_pc)) {
    out[[paste0("PC", i, "_loading")]] <- var_info$coord[, i]
    out[[paste0("PC", i, "_contrib")]] <- var_info$contrib[, i]
    out[[paste0("PC", i, "_cos2")]] <- var_info$cos2[, i]
  }
  
  contrib_cols <- grep("_contrib$", names(out), value = TRUE)
  out$Total_contrib <- rowSums(out[, contrib_cols, drop = FALSE])
  
  out %>%
    arrange(desc(Total_contrib))
}

# Extract PCA sample and trait scores
extract_pca_scores <- function(pca_obj, data, category, n_pcs = 2) {
  sample_scores <- as.data.frame(pca_obj$x[, 1:n_pcs, drop = FALSE])
  colnames(sample_scores) <- paste0("PC", seq_len(n_pcs))
  sample_scores$Category <- category
  sample_scores <- cbind(data[, 1:7], sample_scores)
  
  trait_scores <- as.data.frame(pca_obj$rotation[, 1:n_pcs, drop = FALSE])
  colnames(trait_scores) <- paste0("PC", seq_len(n_pcs))
  trait_scores$Trait <- rownames(pca_obj$rotation)
  trait_scores$Category <- category
  
  list(samples = sample_scores, traits = trait_scores)
}

# Run PCA
run_pca <- function(trait_matrix) {
  prcomp(trait_matrix, center = TRUE, scale. = TRUE)
}

# Run RDA and return tidy permutation table
run_rda <- function(response_matrix, explanatory_matrix, permutations = 999) {
  rda_model <- rda(response_matrix ~ Intensity * Frequency * Composition,
                   data = explanatory_matrix,scale = TRUE)
  
  anova_terms <- as.data.frame(anova(rda_model, by = "term", permutations = permutations))
  
  anova_terms$Term <- rownames(anova_terms)
  
  total_inertia <- rda_model$tot.chi
  
  rda_table <- anova_terms %>%
    mutate(Percent_total_variance = 100 * Variance / total_inertia) %>%
    select(Term, Df, Variance, Percent_total_variance, F, `Pr(>F)`) %>%
    rename(df = Df, Variance_explained = Variance, F_value = F, p_value = `Pr(>F)`)
  
  list(model = rda_model, table = rda_table,
       overall_test = anova(rda_model, permutations = permutations),
       adj_r2 = RsquareAdj(rda_model))
}

# ---------------------------------------------------------
# 3. Overall PCA
# ---------------------------------------------------------

pca_all <- run_pca(traits_all)

# Save sample and trait scores
pca_all_scores <- extract_pca_scores(pca_all, dat, category = "All", n_pcs = 2)

write.csv(pca_all_scores$samples, file.path(output_dir, "PCA_all_sample_scores.csv"), row.names = FALSE)

write.csv(pca_all_scores$traits, file.path(output_dir, "PCA_all_trait_scores.csv"), row.names = FALSE)

# Save PCA variable summary
pca_all_table <- pca_variable_table(pca_all, n_pc = 2)

write.csv(pca_all_table, file.path(output_dir, "PCA_all_trait_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------
# 4. Overall RDA
# ---------------------------------------------------------

set.seed(210)

rda_all <- run_rda(traits_all, design_matrix)

write.csv(rda_all$table, file.path(output_dir, "RDA_all_traits.csv"), row.names = FALSE)

# Save overall test and adjusted R2
rda_all_overall <- tibble(Model = "All traits", F_value = rda_all$overall_test$F[1],
                          p_value = rda_all$overall_test$`Pr(>F)`[1],
                          Adjusted_R2 = rda_all$adj_r2$adj.r.squared)

write.csv(rda_all_overall, file.path(output_dir, "RDA_all_traits_overall_test.csv"), row.names = FALSE)

# ---------------------------------------------------------
# 5. Category-specific PCA
# ---------------------------------------------------------

pca_growth <- run_pca(traits_growth)
pca_leaf <- run_pca(traits_leaf)
pca_physiological <- run_pca(traits_physiological)
pca_root <- run_pca(traits_root)
pca_fitness <- run_pca(traits_fitness)

# Extract scores
scores_growth <- extract_pca_scores(pca_growth, dat, "Growth", 2)
scores_leaf <- extract_pca_scores(pca_leaf, dat, "Leaf", 2)
scores_physiological <- extract_pca_scores(pca_physiological, dat, "Physiological", 2)
scores_root <- extract_pca_scores(pca_root, dat, "Root", 2)
scores_fitness <- extract_pca_scores(pca_fitness, dat, "Fitness", 2)

# Combine and save sample scores
pca_category_samples <- bind_rows(scores_growth$samples,
                                  scores_leaf$samples,
                                  scores_physiological$samples,
                                  scores_root$samples,
                                  scores_fitness$samples)

write.csv(pca_category_samples, file.path(output_dir, "PCA_category_sample_scores.csv"), row.names = FALSE)

# Combine and save trait scores
pca_category_traits <- bind_rows(scores_growth$traits,
                                 scores_leaf$traits,
                                 scores_physiological$traits,
                                 scores_root$traits,
                                 scores_fitness$traits)

write.csv(pca_category_traits, file.path(output_dir, "PCA_category_trait_scores.csv"), row.names = FALSE)

# Save category-specific PCA variable summaries
pca_category_summary <- bind_rows(pca_variable_table(pca_growth, 2) %>% mutate(Category = "Growth"),
                                  pca_variable_table(pca_leaf, 2) %>% mutate(Category = "Leaf"),
                                  pca_variable_table(pca_physiological, 2) %>% mutate(Category = "Physiological"),
                                  pca_variable_table(pca_root, 2) %>% mutate(Category = "Root"),
                                  pca_variable_table(pca_fitness, 2) %>% mutate(Category = "Fitness")) %>%
  select(Category, everything())

write.csv(pca_category_summary, file.path(output_dir, "PCA_category_trait_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------
# 6. Category-specific RDA
# ---------------------------------------------------------

set.seed(123)

rda_growth <- run_rda(traits_growth, design_matrix)
rda_leaf <- run_rda(traits_leaf, design_matrix)
rda_physiological <- run_rda(traits_physiological, design_matrix)
rda_root <- run_rda(traits_root, design_matrix)
rda_fitness <- run_rda(traits_fitness, design_matrix)

rda_category_table <- bind_rows(rda_growth$table %>% mutate(Category = "Growth"),
                                rda_leaf$table %>% mutate(Category = "Leaf"),
                                rda_physiological$table %>% mutate(Category = "Physiological"),
                                rda_root$table %>% mutate(Category = "Root"),
                                rda_fitness$table %>% mutate(Category = "Fitness")) %>%
  select(Category, everything())

write.csv(rda_category_table, file.path(output_dir, "RDA_category_traits.csv"), row.names = FALSE)

# Save category-specific overall RDA tests
rda_category_overall <- bind_rows(tibble(Category = "Growth",
                                         F_value = rda_growth$overall_test$F[1],
                                         p_value = rda_growth$overall_test$`Pr(>F)`[1],
                                         Adjusted_R2 = rda_growth$adj_r2$adj.r.squared),
                                  tibble(Category = "Leaf",
                                         F_value = rda_leaf$overall_test$F[1],
                                         p_value = rda_leaf$overall_test$`Pr(>F)`[1],
                                         Adjusted_R2 = rda_leaf$adj_r2$adj.r.squared),
                                  tibble(Category = "Physiological",
                                         F_value = rda_physiological$overall_test$F[1],
                                         p_value = rda_physiological$`Pr(>F)`[1],
                                         Adjusted_R2 = rda_physiological$adj_r2$adj.r.squared),
                                  tibble(Category = "Root",
                                         F_value = rda_root$overall_test$F[1],
                                         p_value = rda_root$overall_test$`Pr(>F)`[1],
                                         Adjusted_R2 = rda_root$adj_r2$adj.r.squared),
                                  tibble(Category = "Fitness",
                                         F_value = rda_fitness$overall_test$F[1],
                                         p_value = rda_fitness$overall_test$`Pr(>F)`[1],
                                         Adjusted_R2 = rda_fitness$adj_r2$adj.r.squared))

# Save - Table S12
write.csv(rda_category_overall, file.path(output_dir, "RDA_category_overall_tests.csv"), row.names = FALSE)
