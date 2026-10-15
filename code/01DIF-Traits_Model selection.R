############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Comparison of alternative models
# Date: June 2026
############################################################

## Load packages
library(dplyr)
library(tidyr)
library(tibble)

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

## Traits analysed with linear models
traits <- c("N", "H", "SB", "TB", "R_S",
            "SLA", "LDMC", "LNC", "Leaf_CN",
            "Chl", "PSIIeff", "Psi_osm", "SS", "SD",
            "RD", "RDMC", "SRL", "SRA", "RTD", "RNC", "Root_CN")

## Function to compare models for one trait
compare_models <- function(trait) {
  
  ## Use the same rows for all candidate models of a given trait
  dat_trait <- dat %>%
    select(all_of(c(trait,"Amount", "Frequency", "Composition",
                    "TB_pot_c", "Bromus_N_2022_c"))) %>%
    na.omit()
  
  model_baseline <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition")),
                       data = dat_trait)
  
  model_biomass <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition + TB_pot_c")),
                      data = dat_trait)
  
  model_density <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition + Bromus_N_2022_c")),
                      data = dat_trait)
  
  model_biomass_density <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition + TB_pot_c + Bromus_N_2022_c")),
                              data = dat_trait)
  
  model_biomass_replacing_composition <- lm(as.formula(paste(trait, "~ Amount * Frequency * TB_pot_c")),
                                            data = dat_trait)
  
  tibble(Trait = trait, Model = c("Baseline", "Community biomass", "Initial density",
                                  "Community biomass + initial density",
                                  "Community biomass replacing composition"),
         AIC = c(AIC(model_baseline), AIC(model_biomass), AIC(model_density),
                 AIC(model_biomass_density), AIC(model_biomass_replacing_composition)),
         N_obs = nrow(dat_trait))
}

## Run model comparison for all traits
model_comp <- bind_rows(lapply(traits, compare_models))

## Rank models within each trait and calculate ΔAIC
model_ranked <- model_comp %>%
  group_by(Trait) %>%
  mutate(Delta_AIC = AIC - min(AIC), Rank_AIC = rank(AIC, ties.method = "first")) %>%
  ungroup()

## Extract best-supported model for each trait
best_model_summary <- model_ranked %>%
  filter(Rank_AIC == 1) %>%
  transmute(Trait, N_obs, Best_AIC_model = Model)

## Create one-row-per-trait ΔAIC summary table
delta_aic_summary <- model_ranked %>%
  select(Trait, Model, Delta_AIC) %>%
  pivot_wider(names_from = Model, values_from = Delta_AIC, names_prefix = "Delta_AIC_")

readable_summary <- best_model_summary %>%
  left_join(delta_aic_summary, by = "Trait") %>%
  select(Trait, N_obs, Best_AIC_model, `Delta_AIC_Baseline`,
         `Delta_AIC_Community biomass`, `Delta_AIC_Initial density`,
         `Delta_AIC_Community biomass + initial density`,
         `Delta_AIC_Community biomass replacing composition`) %>%
  arrange(Trait)

## Optional: count how many traits were best supported by each model
best_model_counts <- readable_summary %>%
  count(Best_AIC_model, name = "Number_of_traits") %>%
  arrange(desc(Number_of_traits))

## Save outputs
write.csv(readable_summary, file.path(output_dir, "Trait_ModelComparison_AICSummary.csv"), row.names = FALSE)

write.csv(best_model_counts, file.path(output_dir, "Trait_BestModel_Counts.csv"), row.names = FALSE)



## Additional sensitivity analysis:
## Initial Bromus density and initial Inula density

## - Conspecific pots are assigned Inula density = 0
## - Heterospecific pots keep observed Inula density
## - Heterospecific pots with no established/surviving Inula are also 0
##
## This treats Inula density as realised heterospecific-neighbour density.

dat_inula <- dat %>%
  mutate(Inula_N_2022 = case_when(Composition == "Conspecific" ~ 0,
                                  Composition == "Heterospecific" & is.na(N2022_Inula) ~ 0,
                                  Composition == "Heterospecific" ~ N2022_Inula),
         Inula_N_2022_c = Inula_N_2022 - mean(Inula_N_2022, na.rm = TRUE))


## Function to compare baseline, Bromus density, Inula density,
## and Bromus + Inula density models for one trait
compare_bromus_inula_density_models <- function(trait) {
  
  dat_trait <- dat_inula %>%
    select(all_of(c(trait, "Amount", "Frequency", "Composition",
                    "Bromus_N_2022_c", "Inula_N_2022_c"))) %>%
    na.omit() %>%
    droplevels()
  
  model_baseline <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition")),
                       data = dat_trait)
  
  model_bromus_density <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition + Bromus_N_2022_c")),
                             data = dat_trait)
  
  model_inula_density <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition + Inula_N_2022_c")),
                            data = dat_trait)
  
  model_bromus_inula_density <- lm(as.formula(paste(trait, "~ Amount * Frequency * Composition + Bromus_N_2022_c + Inula_N_2022_c")),
                                   data = dat_trait)
  
  tibble(Trait = trait, Model = c("Baseline", "Initial Bromus density", "Initial Inula density",
                                  "Initial Bromus + Inula density"),
         AIC = c(AIC(model_baseline), AIC(model_bromus_density),
                 AIC(model_inula_density), AIC(model_bromus_inula_density)),
         N_obs = nrow(dat_trait))
}


## Run model comparison for all traits
density_model_comp <- bind_rows(lapply(traits, compare_bromus_inula_density_models))

## Rank models within each trait and calculate ΔAIC
density_model_ranked <- density_model_comp %>%
  group_by(Trait) %>%
  mutate(Delta_AIC = AIC - min(AIC, na.rm = TRUE),
         Rank_AIC = rank(AIC, ties.method = "first", na.last = "keep")) %>%
  ungroup()

## Extract best-supported model for each trait
density_best_model_summary <- density_model_ranked %>%
  filter(Rank_AIC == 1) %>%
  transmute(Trait, N_obs, Best_AIC_model = Model)

## Create one-row-per-trait ΔAIC summary table
density_delta_aic_summary <- density_model_ranked %>%
  select(Trait, Model, Delta_AIC) %>%
  pivot_wider(names_from = Model, values_from = Delta_AIC, names_prefix = "Delta_AIC_")

density_readable_summary <- density_best_model_summary %>%
  left_join(density_delta_aic_summary, by = "Trait") %>%
  select(Trait, N_obs, Best_AIC_model, `Delta_AIC_Baseline`,
         `Delta_AIC_Initial Bromus density`,
         `Delta_AIC_Initial Inula density`,
         `Delta_AIC_Initial Bromus + Inula density`) %>%
  arrange(Trait)


## Count how many traits were best supported by each model
density_best_model_counts <- density_readable_summary %>%
  count(Best_AIC_model, name = "Number_of_traits") %>%
  arrange(desc(Number_of_traits))

## Save outputs
write.csv(density_readable_summary, file.path(output_dir, "Trait_Bromus_InulaDensity_AICSummary.csv"), row.names = FALSE)

write.csv(density_best_model_counts, file.path(output_dir, "Trait_Bromus_InulaDensity_BestModel_Counts.csv"), row.names = FALSE)
