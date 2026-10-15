############################################################
# Project: Drought intensity and frequency (DIF) - Traits
# Analysis: Inula responses in heterospecific pots
# Date: June 2026
############################################################

## Load packages
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(car)
library(glmmTMB)
library(DHARMa)

setwd("~/Desktop/DIF/Traits/Data availability")

## File paths
input_file <- "data/DIF-Traits_data.csv"
output_dir <- "results/tables"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Import and prepare data
dat <- read.csv(input_file, header = TRUE) %>%
  mutate(Amount = factor(Amount, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")),
         Composition = factor(Composition, levels = c("Conspecific", "Heterospecific")))

## Set contrasts for Type III tests
options(contrasts = c("contr.sum", "contr.poly"))

## Keep only heterospecific communities, where Inula was planted
inula_dat <- dat %>%
  filter(Composition == "Heterospecific") %>%
  droplevels()

## Check treatment replication
table(inula_dat$Amount, inula_dat$Frequency)

## 1. Linear models: Inula performance and functional traits
## For N_Inula, SB_Inula, and TB_Inula, zero values are retained, where no Inula survived.
##
## For R_S_Inula, SLA_Inula, LDMC_Inula, Chl_Inula, and
## PSIIeff_Inula, analyses are based on surviving Inula only.
## The interaction term is omitted because the Dry-High treatment
## has no surviving Inula trait measurements.

inula_lm_traits <- c("N_Inula", "SB_Inula", "TB_Inula",
                     "R_S_Inula", "SLA_Inula", "LDMC_Inula",
                     "Chl_Inula", "PSIIeff_Inula")

## Keep only traits present in the dataset
inula_lm_traits <- inula_lm_traits[inula_lm_traits %in% names(inula_dat)]

## Traits for which the full Amount x Frequency model can be fitted
inula_performance_traits <- c("N_Inula", "SB_Inula", "TB_Inula")
inula_performance_traits <- inula_performance_traits[inula_performance_traits %in% inula_lm_traits]

## Traits measured only for surviving Inula
inula_survivor_traits <- setdiff(inula_lm_traits, inula_performance_traits)


## Full factorial models for Inula performance responses
inula_performance_results <- map_dfr(inula_performance_traits, function(trait) {
  
  df <- inula_dat %>%
    select(all_of(c(trait, "Amount", "Frequency"))) %>%
    na.omit() %>%
    droplevels()
  
  contrasts(df$Amount) <- contr.sum
  contrasts(df$Frequency) <- contr.sum
  
  model_formula <- as.formula(paste(trait, "~ Amount * Frequency"))
  model <- lm(model_formula, data = df)
  
  as.data.frame(Anova(model, type = 3)) %>%
    rownames_to_column("Effect") %>%
    mutate(Trait = trait, N_obs = nrow(df), Formula = deparse(model_formula), .before = 1)
})


## Additive models for functional traits measured among surviving Inula
inula_survivor_results <- map_dfr(inula_survivor_traits, function(trait) {
  
  df <- inula_dat %>%
    select(all_of(c(trait, "Amount", "Frequency"))) %>%
    na.omit() %>%
    droplevels()
  
  contrasts(df$Amount) <- contr.sum
  contrasts(df$Frequency) <- contr.sum
  
  model_formula <- as.formula(paste(trait, "~ Amount + Frequency"))
  model <- lm(model_formula, data = df)
  
  as.data.frame(Anova(model, type = 3)) %>%
    rownames_to_column("Effect") %>%
    mutate(Trait = trait, N_obs = nrow(df), Formula = deparse(model_formula),
      Note = "Functional trait analysed among surviving Inula only; interaction omitted because Dry-High had no surviving Inula trait measurements",
      .before = 1)
})


## Combine and save linear model results
inula_lm_results <- bind_rows(inula_performance_results %>%
                                mutate(Note = NA_character_), inula_survivor_results)
inula_lm_results

write.csv(inula_lm_results, file.path(output_dir, "Inula_Performance_Traits_ANOVA.csv"), row.names = FALSE)

## 2. Beta regression: Inula SER and SSR
## Smithson & Verkuilen adjustment to move proportions into (0, 1)
squeeze01 <- function(x) {
  n <- sum(!is.na(x))
  (x * (n - 1) + 0.5) / n
}

## Prepare data for one proportional response
prepare_inula_beta_data <- function(data, response) {
  
  df <- data %>%
    filter(!is.na(.data[[response]])) %>%
    droplevels()
  
  contrasts(df$Amount) <- contr.sum
  contrasts(df$Frequency) <- contr.sum
  
  df
}

## Check whether all Amount x Frequency treatment combinations are present
has_full_2way <- function(data) {
  all(xtabs(~ Amount + Frequency, data) > 0)
}

## Fit beta regression for one proportional response
fit_inula_beta_trait <- function(data, response) {
  
  df <- prepare_inula_beta_data(data, response)
  
  response_beta <- paste0(response, "_beta")
  
  df[[response_beta]] <- if (any(df[[response]] <= 0 | df[[response]] >= 1, na.rm = TRUE)) {
    squeeze01(df[[response]])
  } else {
    df[[response]]
  }
  
  formula_full <- as.formula(paste(response_beta, "~ Amount * Frequency"))
  formula_additive <- as.formula(paste(response_beta, "~ Amount + Frequency"))
  
  model_formula <- if (has_full_2way(df)) formula_full else formula_additive
  
  model <- glmmTMB(model_formula, data = df, family = beta_family(link = "logit"))
  
  sim_res <- simulateResiduals(model, n = 1000)
  
  list(model = model, anova = Anova(model, type = 3), diagnostics = sim_res,
       formula_used = deparse(model_formula))
}

## Replace names here if needed
inula_beta_traits <- c("SER_Inula", "SSR2_Inula")
inula_beta_traits <- inula_beta_traits[inula_beta_traits %in% names(inula_dat)]

inula_beta_results <- map_dfr(inula_beta_traits, function(trait) {
  
  fit <- fit_inula_beta_trait(inula_dat, trait)
  
  as.data.frame(fit$anova) %>%
    rownames_to_column("Effect") %>%
    mutate(Trait = trait, Formula = fit$formula_used, .before = 1)
})

inula_beta_results

write.csv(inula_beta_results, file.path(output_dir, "Inula_SER_SSR_beta_regression.csv"), row.names = FALSE)

## 3. Model-based plot: key Inula responses
library(emmeans)
library(ggplot2)
library(ggpubr)

## Output directory for figures
figure_dir <- "results/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

## Helper for emmeans output
fix_emm_names <- function(df) {
  if ("asymp.LCL" %in% names(df)) names(df)[names(df) == "asymp.LCL"] <- "lower.CL"
  if ("asymp.UCL" %in% names(df)) names(df)[names(df) == "asymp.UCL"] <- "upper.CL"
  if ("emmean" %in% names(df)) names(df)[names(df) == "emmean"] <- "response"
  df
}

## Helper for significance stars
p_to_stars <- function(p) {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01  ~ "**",
            p < 0.05  ~ "*", TRUE ~ "ns")
}

## Prepare beta-regression variables
inula_dat <- inula_dat %>%
  mutate(SER_Inula_beta = ifelse(SER_Inula <= 0 | SER_Inula >= 1, squeeze01(SER_Inula), SER_Inula),
         SSR2_Inula_beta = ifelse(SSR2_Inula <= 0 | SSR2_Inula >= 1, squeeze01(SSR2_Inula), SSR2_Inula))

## Fit models used for plotting
m_SER_plot <- glmmTMB(SER_Inula_beta ~ Amount * Frequency, data = inula_dat %>% filter(!is.na(SER_Inula_beta)),
                      family = beta_family(link = "logit"))

m_SSR_plot <- glmmTMB(SSR2_Inula_beta ~ Amount * Frequency, data = inula_dat %>% filter(!is.na(SSR2_Inula_beta)),
                      family = beta_family(link = "logit"))

m_N_plot <- lm(N_Inula ~ Amount * Frequency, data = inula_dat %>% filter(!is.na(N_Inula)))

## Extract model-predicted marginal means

## Main Amount effects
emm_N_amount <- summary(emmeans(m_N_plot, ~ Amount), infer = TRUE) %>%
  as.data.frame() %>%
  fix_emm_names() %>%
  mutate(Amount = factor(Amount, levels = c("Wet", "Dry")))

emm_SER_amount <- summary(emmeans(m_SER_plot, ~ Amount, type = "response"), infer = TRUE) %>%
  as.data.frame() %>%
  fix_emm_names() %>%
  mutate(Amount = factor(Amount, levels = c("Wet", "Dry")))

emm_SSR_amount <- summary(emmeans(m_SSR_plot, ~ Amount, type = "response"), infer = TRUE) %>%
  as.data.frame() %>%
  fix_emm_names() %>%
  mutate(Amount = factor(Amount, levels = c("Wet", "Dry")))

## SER Amount × Frequency interaction
emm_SER_interaction <- summary(emmeans(m_SER_plot, ~ Amount * Frequency, type = "response"), infer = TRUE) %>%
  as.data.frame() %>%
  fix_emm_names() %>%
  mutate(Amount = factor(Amount, levels = c("Wet", "Dry")),
         Frequency = factor(Frequency, levels = c("Low", "High")))

## P-values from the saved result tables
## Linear model table: N_Inula
p_N_amount <- inula_lm_results %>%
  filter(Trait == "N_Inula", Effect == "Amount") %>%
  pull(`Pr(>F)`)

## Beta-regression table: SER_Inula and SSR2_Inula
p_SER_amount <- inula_beta_results %>%
  filter(Trait == "SER_Inula", Effect == "Amount") %>%
  pull(`Pr(>Chisq)`)

p_SER_interaction <- inula_beta_results %>%
  filter(Trait == "SER_Inula", Effect == "Amount:Frequency") %>%
  pull(`Pr(>Chisq)`)

p_SSR_amount <- inula_beta_results %>%
  filter(Trait == "SSR2_Inula", Effect == "Amount") %>%
  pull(`Pr(>Chisq)`)

## Check p-values
p_N_amount
p_SER_amount
p_SER_interaction
p_SSR_amount
p_to_stars(p_SER_amount)

## Star positions for Amount panels
star_N <- emm_N_amount %>%
  summarise(x = 1.5, y = max(upper.CL, na.rm = TRUE) + 0.10 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
            sig_lab = p_to_stars(p_N_amount))

star_SER_amount <- emm_SER_amount %>%
  summarise(x = 1.5, y = max(upper.CL, na.rm = TRUE) + 0.10 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
            sig_lab = p_to_stars(p_SER_amount))

star_SSR_amount <- emm_SSR_amount %>%
  summarise(x = 1.5, y = max(upper.CL, na.rm = TRUE) + 0.10 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
            sig_lab = p_to_stars(p_SSR_amount))

## Stars for SER Amount × Frequency panel
star_SER_interaction <- emm_SER_interaction %>%
  summarise(x = 1.5, y = max(upper.CL, na.rm = TRUE) + 0.22 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
            sig_lab = paste0("A × F: ", p_to_stars(p_SER_interaction)))

## Add Wet–Dry comparisons within each Frequency level
SER_amount_within_freq <- emmeans(m_SER_plot, ~ Amount | Frequency, type = "response")

SER_amount_within_freq_contrasts <- summary(pairs(SER_amount_within_freq), infer = TRUE) %>%
  as.data.frame()

stars_SER_within_freq <- emm_SER_interaction %>%
  group_by(Frequency) %>%
  summarise(y = max(upper.CL, na.rm = TRUE) + 0.10 * diff(range(c(lower.CL, upper.CL), na.rm = TRUE)),
            .groups = "drop") %>%
  left_join(SER_amount_within_freq_contrasts %>%
              mutate(Frequency = factor(Frequency, levels = c("Low", "High"))) %>%
              transmute(Frequency, sig_lab = p_to_stars(p.value)),by = "Frequency")

## Colours and theme
col_amount <- c(Wet = "#1994E5", Dry = "#E7782F")

theme_inula <- theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 11, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 11, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size =10))

## Panel 1: N_Inula by Amount
p_N <- ggplot(emm_N_amount, aes(x = Amount, y = response, color = Amount)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.7) +
  geom_text(data = star_N, aes(x = x, y = y, label = sig_lab),
            inherit.aes = FALSE, fontface = "bold", size = 4.5) +
  scale_color_manual(values = col_amount, guide = "none") +
  labs(x = "Water amount", y = "Model-predicted value") +
  ggtitle("Inula plant number") +
  coord_cartesian(clip = "off") +
  theme_inula
p_N

## Panel 2: SER_Inula by Amount
p_SER_amount <- ggplot(emm_SER_amount, aes(x = Amount, y = response, color = Amount)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.7) +
  geom_text(data = star_SER_amount, aes(x = x, y = y, label = sig_lab),
            inherit.aes = FALSE, fontface = "bold", size = 4.5) +
  scale_color_manual(values = col_amount, guide = "none") +
  labs(x = "Water amount", y = "Model-predicted value") +
  ggtitle("Seedling establishment rate") +
  coord_cartesian(clip = "off") +
  theme_inula
p_SER_amount


## Panel 3: SSR2_Inula by Amount
p_SSR <- ggplot(emm_SSR_amount, aes(x = Amount, y = response, color = Amount)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.7) +
  geom_text(data = star_SSR_amount, aes(x = x, y = y, label = sig_lab),
            inherit.aes = FALSE, fontface = "bold", size = 4.5) +
  scale_color_manual(values = col_amount, guide = "none") +
  labs(x = "Water amount", y = "Model-predicted value") +
  ggtitle("Seedling survival rate") +
  coord_cartesian(clip = "off") +
  theme_inula
p_SSR


## Panel 4: SER_Inula Amount × Frequency
## Coloured by Amount
p_SER_interaction <- ggplot(emm_SER_interaction,
                            aes(x = Frequency, y = response, color = Amount, group = Amount)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.12)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.7,
                  position = position_dodge(width = 0.12)) +
  geom_text(data = stars_SER_within_freq, aes(x = Frequency, y = y, label = sig_lab),
            inherit.aes = FALSE, fontface = "bold", size = 4.5) +
  geom_text(data = star_SER_interaction, aes(x = x, y = y, label = sig_lab),
            inherit.aes = FALSE, fontface = "bold", size = 4) +
  scale_color_manual(values = col_amount, name = "Water amount") +
  labs(x = "Watering frequency", y = "Model-predicted value") +
  ggtitle("Seedling establishment rate") +
  coord_cartesian(clip = "off") +
  theme_inula +
  theme(legend.position = c(0.78, 0.78),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.key = element_blank())
p_SER_interaction


## Combine panels
p_inula_model <- ggarrange(p_N, p_SER_amount, p_SSR, p_SER_interaction, ncol = 2, nrow = 2,
                           labels = c("(a)", "(b)", "(c)", "(d)"))
p_inula_model

## Save figure
ggsave(filename = file.path(figure_dir, "FigureS6_Inula_model_predicted_responses.pdf"),
       plot = p_inula_model, width = 8, height = 7, units = "in", device = cairo_pdf)

ggsave(filename = file.path(figure_dir, "FigureS6_Inula_model_predicted_responses.png"),
       plot = p_inula_model, width = 8, height = 7, units = "in", dpi = 300, bg = "white")

