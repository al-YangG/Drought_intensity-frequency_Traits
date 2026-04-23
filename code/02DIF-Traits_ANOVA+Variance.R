############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Type III ANOVA for functional traits
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(car)
library(broom)
library(purrr)
library(effectsize)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/tables"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Import and prepare data
dat <- read.csv(input_file, header = TRUE) %>%
  rename(Psi_osm = `Ψosm`) %>%
  mutate(Intensity = factor(Intensity, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         TB_pot_c = TB_pot - mean(TB_pot, na.rm = TRUE))

# Set contrasts for Type III ANOVA
options(contrasts = c("contr.sum", "contr.poly"))

# Traits analysed with linear models
traits <- c("N", "H", "SB", "RB", "TB", "R_S",
            "SLA", "LDMC", "LCC", "LNC", "Leaf_CN",
            "Chl", "PSIIeff", "Psi_osm", "SS", "SD",
            "RD", "RDMC", "SRL", "SRA", "RTD", "RCC", "RNC", "Root_CN")

# Function to fit model and return tidy ANOVA table with effect sizes
run_trait_anova <- function(trait, formula_rhs, data = dat) {
  model_formula <- as.formula(paste(trait, "~", formula_rhs))
  model <- lm(model_formula, data = data)
  anova_table <- car::Anova(model, type = 3)
  
  results <- broom::tidy(anova_table) %>%
    rename(Sum_Sq = sumsq, DF = df, F = statistic, p = p.value) %>%
    filter(term != "(Intercept)")
  
  total_ss <- sum(results$Sum_Sq, na.rm = TRUE)
  
  results <- results %>%
    mutate(Variance_Explained = 100 * Sum_Sq / total_ss)
  
  partial_eta2 <- effectsize::eta_squared(anova_table, partial = TRUE) %>%
    as.data.frame() %>%
    transmute(term = Parameter, Partial_Eta2 = Eta2_partial)
  
  left_join(results, partial_eta2, by = "term") %>%
    mutate(Trait = trait, .before = 1)
}

# Function to run one model across all traits
run_model_set <- function(formula_rhs) {
  map(traits, ~ run_trait_anova(.x, formula_rhs, dat)) %>%
    list_rbind()
}

## 1. Main model: Intensity * Frequency * Composition
results_ifc <- run_model_set("Intensity * Frequency * Composition") %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm"),
         p_adj_fdr = p.adjust(p, method = "fdr")) %>%
  select(Trait, term, DF, Sum_Sq, F, p, p_adj_holm, p_adj_fdr,
         Variance_Explained, Partial_Eta2)

write.csv(results_ifc,file.path(output_dir, "ANOVA_Intensity_Frequency_Composition.csv"),
          row.names = FALSE) # Table S6

## 2. Alternative model: Intensity * Frequency * Composition + TB_pot_c
results_ifcb <- run_model_set("Intensity * Frequency * Composition + TB_pot_c") %>%
  select(Trait, term, DF, Sum_Sq, F, p, Variance_Explained, Partial_Eta2)

write.csv(results_ifcb, file.path(output_dir, "ANOVA_Intensity_Frequency_Composition_Biomass.csv"),
          row.names = FALSE) # Table S8

## 3. Alternative model: Intensity * Frequency * TB_pot_c
results_ifb <- run_model_set("Intensity * Frequency * TB_pot_c") %>%
  select(Trait, term, DF, Sum_Sq, F, p, Variance_Explained, Partial_Eta2)

write.csv(results_ifb, file.path(output_dir, "ANOVA_Intensity_Frequency_Biomass.csv"),
          row.names = FALSE) # Table S9
