############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Fitness traits
# Author: Yang
# Date: April 2026
############################################################

# Load packages
library(dplyr)
library(tibble)
library(purrr)
library(car)
library(glmmTMB)
library(DHARMa)

setwd("~/Desktop/DIF/Traits/Data availability")

# File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/tables"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Import data
dat <- read.csv(input_file, header = TRUE) %>%
  mutate(Intensity = factor(Intensity, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")))

# Set contrasts for Type III tests
options(contrasts = c("contr.sum", "contr.poly"))

# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------

# Smithson & Verkuilen adjustment to move proportions into (0, 1)
squeeze01 <- function(x) {
  n <- sum(!is.na(x))
  (x * (n - 1) + 0.5) / n
}

# Prepare data for one response variable
prepare_data <- function(data, response) {
  df <- data %>%
    filter(!is.na(.data[[response]]))
  
  contrasts(df$Intensity) <- contr.sum
  contrasts(df$Frequency) <- contr.sum
  contrasts(df$Composition) <- contr.sum
  
  df
}

# Check whether all treatment combinations are present
has_full_3way <- function(data) {
  all(xtabs(~ Intensity + Frequency + Composition, data) > 0)
}

# Convert p-values to significance labels
p_to_label <- function(p) {
  case_when(is.na(p) ~ NA_character_, p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", TRUE ~ "ns")
}

# ---------------------------------------------------------
# 1. Beta regression: SER, SSR1, SSR2
# ---------------------------------------------------------

fit_beta_trait <- function(data, response) {
  df <- prepare_data(data, response)
  
  response_beta <- paste0(response, "_beta")
  df[[response_beta]] <- if (any(df[[response]] <= 0 | df[[response]] >= 1, na.rm = TRUE)) {
    squeeze01(df[[response]])
  } else {
    df[[response]]
  }
  
  formula_3way <- as.formula(paste(response_beta, "~ Intensity * Frequency * Composition"))
  formula_2way <- as.formula(paste(response_beta, "~ (Intensity + Frequency + Composition)^2"))
  
  model_formula <- if (has_full_3way(df)) formula_3way else formula_2way
  
  model <- glmmTMB(model_formula, data = df, family = beta_family(link = "logit"))
  
  list(model = model, anova = Anova(model, type = 3), formula_used = deparse(model_formula))
}

beta_traits <- c("SER", "SSR1", "SSR2")

beta_results <- map_dfr(beta_traits, function(trait) {
  fit <- fit_beta_trait(dat, trait)
  
  as.data.frame(fit$anova) %>%
    rownames_to_column("Effect") %>%
    mutate(Trait = trait, Formula = fit$formula_used, .before = 1)
})

# Save - Table S7
write.csv(beta_results, file.path(output_dir, "Fitness_beta_regression.csv"), row.names = FALSE)

# ---------------------------------------------------------
# 2. Tweedie GLMM: SP1, SP2, SPtotal
# ---------------------------------------------------------

fit_seed_trait <- function(data, response, seed = 42) {
  df <- prepare_data(data, response)
  
  if (all(df[[response]] == 0)) {
    return(list(model = NULL, anova = NULL,
                formula_used = paste(response, "~ (Intensity + Frequency + Composition)^2"),
                diagnostics = NULL, model_type = "All zeros - not estimable"))
  }
  
  model_formula <- as.formula(paste(response, "~ (Intensity + Frequency + Composition)^2"))
  
  model <- glmmTMB(model_formula, data = df, family = tweedie(link = "log"))
  
  set.seed(seed)
  sim_res <- simulateResiduals(model, n = 1000)
  
  list(model = model, anova = Anova(model, type = 3, test.statistic = "Chisq"),
       formula_used = deparse(model_formula), diagnostics = sim_res, model_type = "Tweedie GLMM")
}

seed_traits <- c("SP1", "SP2", "SPtotal")

seed_results <- map_dfr(seed_traits, function(trait) {
  fit <- fit_seed_trait(dat, trait)
  
  if (is.null(fit$anova)) {
    tibble(Trait = trait, Model = fit$model_type, Formula = fit$formula_used,
           Effect = NA, Chisq = NA, Df = NA, `Pr(>Chisq)` = NA)
  } else {
    as.data.frame(fit$anova) %>%
      rownames_to_column("Effect") %>%
      mutate(Trait = trait, Model = fit$model_type, Formula = fit$formula_used, .before = 1)
  }
})

# Save - Table S7
write.csv(seed_results, file.path(output_dir, "Fitness_seed_production.csv"), row.names = FALSE)

# ---------------------------------------------------------
# 3. Linear models: RP1, RP2
# ---------------------------------------------------------

fit_ramet_trait <- function(data, response) {
  df <- prepare_data(data, response)
  
  model_formula <- as.formula(paste(response, "~ Intensity * Frequency * Composition"))
  
  model <- lm(model_formula, data = df)
  
  list(model = model, anova = Anova(model, type = 3), formula_used = deparse(model_formula))
}

ramet_traits <- c("RP1", "RP2")

ramet_results <- map_dfr(ramet_traits, function(trait) {
  fit <- fit_ramet_trait(dat, trait)
  
  as.data.frame(fit$anova) %>%
    rownames_to_column("Effect") %>%
    mutate(Trait = trait, Formula = fit$formula_used, .before = 1)
})

# Save - Table S7
write.csv(ramet_results, file.path(output_dir, "Fitness_ramet_production.csv"), row.names = FALSE)
