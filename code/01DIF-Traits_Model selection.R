############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Comparison of alternative models
# Author: Yang
# Date: April 2026
############################################################

## Load package
library(dplyr)

setwd("~/Desktop/DIF/Traits/Data availability")

## File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/tables"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Import and prepare data
dat <- read.csv(input_file, header = TRUE) %>%
  rename(Psi_osm = `Ψosm`) %>%
  mutate(Intensity = factor(Intensity, levels = c("Dry", "Wet")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")),
         TB_pot_c = TB_pot - mean(TB_pot, na.rm = TRUE))

## Traits analysed with linear models
traits <- c("N", "H", "SB", "RB", "TB", "R_S",
            "SLA", "LDMC", "LCC", "LNC", "Leaf_CN",
            "Chl", "PSIIeff", "Psi_osm", "SS", "SD",
            "RD", "RDMC", "SRL", "SRA", "RTD", "RCC", "RNC", "Root_CN")

## Function to compare models for one trait
compare_models <- function(trait) {
  model1 <- lm(as.formula(paste(trait, "~ Intensity * Frequency * Composition")), data = dat)
  model2 <- lm(as.formula(paste(trait, "~ Intensity * Frequency * Composition + TB_pot_c")), data = dat)
  model3 <- lm(as.formula(paste(trait, "~ Intensity * Frequency * TB_pot_c")), data = dat)
  
  tibble(Trait = trait,
         Model = c("NoBiomass", "WithBiomass", "BiomassOnly"),
         AIC = c(AIC(model1), AIC(model2), AIC(model3)),
         BIC = c(BIC(model1), BIC(model2), BIC(model3)))
}

## Run model comparison for all traits
model_comp <- bind_rows(lapply(traits, compare_models))

## Rank models within each trait
model_ranked <- model_comp %>%
  group_by(Trait) %>%
  arrange(AIC, .by_group = TRUE) %>%
  mutate(Rank_AIC = row_number(),
         Delta_AIC = AIC - min(AIC),
         Rank_BIC = rank(BIC, ties.method = "first")) %>%
  ungroup()

## Extract best model for each trait
best_models <- model_ranked %>%
  filter(Rank_AIC == 1) %>%
  select(Trait, BestModel = Model, AIC, BIC, Delta_AIC)

## Save results - Table S5
write.csv(model_ranked, file.path(output_dir, "Trait_ModelComparison_AICBIC.csv"), row.names = FALSE)
write.csv(best_models, file.path(output_dir, "Trait_BestModels.csv"), row.names = FALSE)

