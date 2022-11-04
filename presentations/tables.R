

library(magrittr)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggmap)
library(ggstance)
library(gridExtra)
library(cowplot)
library(rgeos)
library(sp)
library(knitr)
library(modelr)
library(car)
library(rstan)
library(latex2exp)
library(ggthemes)

library(econometr)
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))

fit_version <- 66

quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

output_basepath = str_glue("temp-data/output_dist_fit{fit_version}")

canva_palette_vibrant <- "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))

load(file.path("data", "takeup_village_pot_dist.RData"))

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data


## Fit Loading

load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))

dist_fit_data %<>%
  left_join(tribble(
    ~ fit_type,        ~ model_color,
      "fit",           "black", 
      "prior-predict", "darkgrey",
  ), by = "fit_type") %>%
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS")))

delta <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))

belief_data = dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  pull(beliefs_results) %>%
  first() 


cols_we_want = c(
    "est_takeup_level",
    "beliefs_results",
    "wtp_results",
    "est_takeup_te",
    "est_takeup_dist_te"
)


dist_fit_data %>%
    filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
    select(est_takeup_te) %>%
    unnest() %>%
    filter(
        mu_assigned_treatment_left == assigned_treatment_left,
        mu_assigned_treatment_right == assigned_treatment_right, 
        assigned_treatment_right == "control", 
        is.na(assigned_dist_group_right))

dist_fit_data %>%
    filter(fit_type == "fit") %>%
    select(any_of(cols_we_want)) %>%
    unnest(cols = c(est_takeup_level))





incentive_te = dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  mutate(
    est_takeup_te =
      map_if(est_takeup_te, fct_match(model_type, "structural"),
             filter, mu_assigned_treatment_left == assigned_treatment_left, mu_assigned_treatment_right == assigned_treatment_right) %>%
        map(filter,
            (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) | (assigned_dist_group_left == assigned_dist_group_right),
            assigned_treatment_left != assigned_treatment_right,
            fct_match(assigned_treatment_right, c("control")),
            fct_match(assigned_treatment_left, "bracelet") | !fct_match(assigned_treatment_right, "calendar")),
    model_color = canva_pal(canva_palette_vibrant)(n())
  ) %>% 
  select(model, model_name, est_takeup_te, fit_type, model_color) %>% 
  mutate(
    est_takeup_te = map(
      est_takeup_te,
      mutate,
      assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
        fct_relabel(str_to_title) %>% 
        fct_relevel("Combined"),
      assigned_treatment_left = str_to_title(assigned_treatment_left)
    )) %>%
    unnest(est_takeup_te) %>%
    filter(assigned_dist_group_left == "Combined" ) %>%
    mutate(combined_incentive_te = round(mean_est, 3)) %>%
    select(treatment = assigned_treatment_left, combined_incentive_te )

signalling_te = dist_fit_data %>% 
  filter(fit_type == "fit") %>%
  select(model, model_name, est_takeup_te, fit_type) %>% 
  mutate(
    est_takeup_te = map(
      est_takeup_te,
      filter,
      (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) | (assigned_dist_group_left == assigned_dist_group_right),
      across(c(assigned_treatment_left, assigned_treatment_right), fct_match, "control"),
      !is.na(mu_assigned_treatment_left),
      fct_match(mu_assigned_treatment_left, "bracelet") | !fct_match(mu_assigned_treatment_right, "calendar"),
      fct_match(mu_assigned_treatment_right, "control"),
    ) %>% 
      map(
        mutate,
        assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
          fct_relabel(str_to_title) %>% 
          fct_relevel("Combined"),
        mu_assigned_treatment_left = str_to_title(mu_assigned_treatment_left),
      ),
    model_color = canva_pal(canva_palette_vibrant)(n())) %>%
    unnest(est_takeup_te) %>%
    filter(assigned_dist_group_left == "Combined") %>%
    select(mu_assigned_treatment_left, mean_est) %>%
    mutate(combined_signalling_te = round(mean_est, 3)) %>%
    select(treatment = mu_assigned_treatment_left, combined_signalling_te)

private_te = dist_fit_data %>% 
  filter(fit_type == "fit") %>%
  select(model, model_name, est_takeup_te, fit_type) %>% 
  mutate(
    est_takeup_te = map(
      est_takeup_te,
      filter,
      (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) | (assigned_dist_group_left == assigned_dist_group_right),
      !is.na(mu_assigned_treatment_left),
      !is.na(mu_assigned_treatment_right),
      across(c(mu_assigned_treatment_left, mu_assigned_treatment_right), fct_match, "control"),
      fct_match(assigned_treatment_left, "bracelet") | !fct_match(assigned_treatment_right, "calendar"),
      fct_match(assigned_treatment_right, "control"),
    ) %>% 
      map(
        mutate,
        assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
          fct_relabel(str_to_title) %>% 
          fct_relevel("Combined"),
        assigned_treatment_left = str_to_title(assigned_treatment_left),
      ),
    model_color = canva_pal(canva_palette_vibrant)(n())) %>%
    unnest(est_takeup_te) %>%
    filter(assigned_dist_group_left == "Combined")  %>%
    select(assigned_treatment_left, mean_est) %>%
    mutate(combined_private_te = round(mean_est, 3)) %>%
    select(treatment = assigned_treatment_left, combined_private_te)


structural_combined_te = inner_join(
    incentive_te,
    private_te,
    by = "treatment"
) %>%
    left_join(signalling_te, by = "treatment")



structural_combined_te

structural_combined_te %>%
    write_csv(
        file.path(output_basepath, 
        str_glue("structural-combined-te.csv"))
    )


