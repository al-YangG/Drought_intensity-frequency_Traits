############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Type III ANOVA for functional traits
# Date: June 2026
############################################################

## Load packages
library(dplyr)
library(car)
library(broom)
library(purrr)
library(effectsize)
library(tidyr)
library(stringr)

setwd("~/Desktop/DIF/Traits/Data availability")

## File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/tables"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Import and prepare data
dat <- read.csv(input_file, header = TRUE) %>%
  rename(Psi_osm = `Ψosm`) %>%
  mutate(Amount = factor(Amount, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         TB_pot_c = TB_pot - mean(TB_pot, na.rm = TRUE),
         Bromus_N_2022_c = N2022 - mean(N2022, na.rm = TRUE))

## Set contrasts for Type III ANOVA
options(contrasts = c("contr.sum", "contr.poly"))

## Traits analysed with linear models
traits <- c("N", "H", "SB", "TB", "R_S",
            "SLA", "LDMC", "LNC", "Leaf_CN",
            "Chl", "PSIIeff", "Psi_osm", "SS", "SD",
            "RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root_CN")

## Function to fit model and return tidy ANOVA table with effect sizes
run_trait_anova <- function(trait, formula_rhs, model_name, data = dat) {
  
  model_formula <- as.formula(paste(trait, "~", formula_rhs))
  
  dat_trait <- data %>%
    select(all_of(c(trait, all.vars(as.formula(paste("~", formula_rhs)))))) %>%
    na.omit()
  
  model <- lm(model_formula, data = dat_trait)
  anova_table <- car::Anova(model, type = 3)
  
  results <- broom::tidy(anova_table) %>%
    rename(Sum_Sq = sumsq, DF = df, F = statistic, p = p.value) %>%
    filter(term != "(Intercept)")
  
  total_ss <- sum(results$Sum_Sq, na.rm = TRUE)
  
  results <- results %>%
    mutate(Variance_Explained = 100 * Sum_Sq / total_ss,
           N_obs = nrow(dat_trait), Model = model_name)
  
  partial_eta2 <- effectsize::eta_squared(anova_table, partial = TRUE) %>%
    as.data.frame() %>%
    transmute(term = Parameter, Partial_Eta2 = Eta2_partial)
  
  left_join(results, partial_eta2, by = "term") %>%
    mutate(Trait = trait, .before = 1) %>%
    relocate(Model, .after = Trait)
}

## Function to run one model across all traits
run_model_set <- function(formula_rhs, model_name, trait_vector = traits, data = dat) {
  map(trait_vector, ~ run_trait_anova(.x, formula_rhs, model_name, data)) %>%
    list_rbind()
}

## Model formulas and short names
model_info <- tibble(Model = c("Baseline", "Biomass-adjusted", "Density-adjusted",
                               "Fully adjusted", "Biomass substitution"),
  Formula = c("Amount * Frequency * Composition",
              "Amount * Frequency * Composition + TB_pot_c",
              "Amount * Frequency * Composition + Bromus_N_2022_c",
              "Amount * Frequency * Composition + TB_pot_c + Bromus_N_2022_c",
              "Amount * Frequency * TB_pot_c"))

## Run all five models
results_baseline <- run_model_set(formula_rhs = "Amount * Frequency * Composition",
                                  model_name = "Baseline", trait_vector = traits) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm"),
         p_adj_fdr = p.adjust(p, method = "fdr")) %>%
  select(Trait, Model, term, DF, Sum_Sq, F, p, p_adj_holm, p_adj_fdr, Variance_Explained, Partial_Eta2, N_obs)

write.csv(results_baseline, file.path(output_dir, "ANOVA_Amount_Frequency_Composition.csv"), row.names = FALSE)

results_biomass <- run_model_set(formula_rhs = "Amount * Frequency * Composition + TB_pot_c",
                                 model_name = "Biomass-adjusted", trait_vector = traits) %>%
  select(Trait, Model, term, DF, Sum_Sq, F, p, Variance_Explained, Partial_Eta2, N_obs)

write.csv(results_biomass, file.path(output_dir, "ANOVA_Amount_Frequency_Composition_Biomass.csv"), row.names = FALSE)

results_density <- run_model_set(formula_rhs = "Amount * Frequency * Composition + Bromus_N_2022_c",
                                 model_name = "Density-adjusted", trait_vector = traits) %>%
  select(Trait, Model, term, DF, Sum_Sq, F, p, Variance_Explained, Partial_Eta2, N_obs)

write.csv(results_density, file.path(output_dir, "ANOVA_Amount_Frequency_Composition_InitialDensity.csv"), row.names = FALSE)

results_full <- run_model_set(formula_rhs = "Amount * Frequency * Composition + TB_pot_c + Bromus_N_2022_c",
                              model_name = "Fully adjusted", trait_vector = traits) %>%
  select(Trait, Model, term, DF, Sum_Sq, F, p, Variance_Explained, Partial_Eta2, N_obs)

write.csv(results_full, file.path(output_dir, "ANOVA_Amount_Frequency_Composition_Biomass_InitialDensity.csv"), row.names = FALSE)

results_biomass_substitution <- run_model_set(formula_rhs = "Amount * Frequency * TB_pot_c", model_name = "Biomass substitution", trait_vector = traits) %>%
  select(Trait, Model, term, DF, Sum_Sq, F, p, Variance_Explained, Partial_Eta2, N_obs)

write.csv(results_biomass_substitution, file.path(output_dir, "ANOVA_Amount_Frequency_BiomassSubstitution.csv"), row.names = FALSE)

## Combine all ANOVA results
all_anova_results <- bind_rows(results_baseline, results_biomass, results_density,
                               results_full, results_biomass_substitution)

write.csv(all_anova_results, file.path(output_dir, "ANOVA_AllModels.csv"), row.names = FALSE)


## Summary table: significant predictors by trait and model

## Function to make terms more readable
clean_term <- function(x) {
  x %>%
    str_replace_all("Amount", "Amount") %>%
    str_replace_all("Frequency", "Frequency") %>%
    str_replace_all("Composition", "Composition") %>%
    str_replace_all("TB_pot_c", "Community biomass") %>%
    str_replace_all("Bromus_N_2022_c", "Initial density")
}

## Long-format summary: one row per trait × model
sig_summary_long <- all_anova_results %>%
  mutate(term_clean = clean_term(term), sig = p < 0.05) %>%
  group_by(Category = case_when(Trait %in% c("N", "H", "SB", "TB", "R_S") ~ "Growth",
                                Trait %in% c("SLA", "LDMC", "LNC", "Leaf_CN") ~ "Leaf",
                                Trait %in% c("Chl", "PSIIeff", "Psi_osm", "SS", "SD") ~ "Physiological",
                                Trait %in% c("RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root_CN") ~ "Root",
                                TRUE ~ "Other"), Trait, Model) %>%
  summarise(Significant_predictors = ifelse(any(sig, na.rm = TRUE), paste(term_clean[sig], collapse = "; "), "None"),
            .groups = "drop") %>%
  arrange(Category, Trait, factor(Model, levels = model_info$Model))

write.csv(sig_summary_long, file.path(output_dir, "ANOVA_SignificantPredictors_Long.csv"), row.names = FALSE)

## Wide-format summary: one row per trait, one column per model
sig_summary_wide <- sig_summary_long %>%
  pivot_wider(names_from = Model, values_from = Significant_predictors) %>%
  select(Category, Trait, Baseline, `Biomass-adjusted`, `Density-adjusted`,
         `Fully adjusted`, `Biomass substitution`)

write.csv(sig_summary_wide, file.path(output_dir, "ANOVA_SignificantPredictors_Wide.csv"), row.names = FALSE)
