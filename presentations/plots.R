
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



fit_version <- 66

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
  ), by = "fit_type")

delta <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))

belief_data = dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  pull(beliefs_results) %>%
  first() 


#### Separate Beliefs ####
belief_data$ate_knows %>%
  filter(assigned_treatment_left != "control") %>%
  mutate(assigned_treatment_left = as.character(assigned_treatment_left)) %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 1,
    crossbar_width = 0.4) +
    labs(subtitle = "Treatment Effects") +
    NULL

ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-belief-te-1ord-plots.png")), width = 7.5, height = 5.0, dpi = 500)

belief_data$ate_knows %>%
  filter(assigned_treatment_left != "control") %>%
  mutate(assigned_treatment_left = as.character(assigned_treatment_left)) %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 2,
    crossbar_width = 0.4) +
    labs(subtitle = "Treatment Effects") +
    NULL

ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-belief-te-2ord-plots.png")), width = 7.5, height = 5.0, dpi = 500)



#### Rate of Change ####
dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST", "STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  transmute(
    model, model_name,
    y_rate_of_change_diff = map2(cluster_roc_diff, stan_data, ~ {
      mutate(.x,
        roc_distance = roc_distance / 1000,
        across(starts_with("per_"), divide_by, sd(.y$analysis_data$cluster.dist.to.pot)), 
        across(starts_with("per_"), multiply_by, 1000),
        across(starts_with("per_"), multiply_by, 100)
      ) 
  })) %>%
  unnest(y_rate_of_change_diff) %>% 
  filter(fct_match(assigned_treatment, c("ink", "bracelet"))) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = assigned_treatment)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.25) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.25) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Difference in Rate of Change",
    x = "Distance to Treatment [km]", y = "Difference in Rate of Change [pp/km]",
    caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
  ) +
#   facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title)) +
  guides( colour = "none") +
  theme_bw() +
  theme(legend.position = "bottom") +
  # coord_cartesian(ylim = c(0, 0.15)) + 
  NULL 
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-diff-roc-by-treat.png")), width = 7.5, height = 5.0, dpi = 500)

dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST", "STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  transmute(
    model, model_name,
    y_rate_of_change_diff = map2(cluster_roc_diff, stan_data, ~ {
      mutate(.x,
        roc_distance = roc_distance / 1000,
        across(starts_with("per_"), divide_by, sd(.y$analysis_data$cluster.dist.to.pot)), 
        across(starts_with("per_"), multiply_by, 1000),
        across(starts_with("per_"), multiply_by, 100)
      ) 
  })) %>%
  unnest(y_rate_of_change_diff) %>% 
  filter(fct_match(assigned_treatment, c("ink", "bracelet"))) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = model_name)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = model_name), alpha = 0.25) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = model_name), alpha = 0.25) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Difference in Rate of Change",
    x = "Distance to Treatment [km]", y = "Difference in Rate of Change [pp/km]",
    caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
  ) +
  facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title)) +
  guides(fill = "none", colour = "none") +
  theme_bw() +
  theme(legend.position = "bottom") +
  # coord_cartesian(ylim = c(0, 0.15)) + 
  NULL 
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-diff-roc-facet-treat.png")), width = 7.5, height = 5.0, dpi = 500)


#### Difference Rep Returns by Dist ####
dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>%
  select(model_name, cluster_rep_return_dist) %>%
  unnest(cluster_rep_return_dist) %>% 
  mutate(
    roc_distance = roc_distance / 1000,
    across(starts_with("per_"), divide_by, 1000)
  ) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = model_name)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = model_name), alpha = 0.4) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = model_name), alpha = 0.4) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Valuation of Reputational Returns in Terms of Distance",
    x = "Distance to Treatment [km]", y = "Distance Value [km]",
  ) +
  facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title), nrow = 1) +
  theme(legend.position = "top") +
  NULL

dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>%
  select(model_name, cluster_rep_return_dist) %>%
  unnest(cluster_rep_return_dist) %>% 
  mutate(
    roc_distance = roc_distance / 1000,
    across(starts_with("per_"), divide_by, 1000)
  ) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = assigned_treatment)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.4) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.4) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Valuation of Reputational Returns in Terms of Distance",
    x = "Distance to Treatment [km]", y = "Distance Value [km]",
  ) +
#   facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title), nrow = 1) +
  theme(legend.position = "bottom") +
  NULL
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-rep-returns-dist-by-treat.png")), width = 7.5, height = 5.0, dpi = 500)

dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>%
  select(model_name, cluster_rep_return_dist) %>%
  unnest(cluster_rep_return_dist) %>% 
  mutate(
    roc_distance = roc_distance / 1000,
    across(starts_with("per_"), divide_by, 1000)
  ) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = assigned_treatment)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.4) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.4) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Valuation of Reputational Returns in Terms of Distance",
    x = "Distance to Treatment [km]", y = "Distance Value [km]",
  ) +
  facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title), nrow = 2) +
  theme_bw() +
  theme(legend.position = "top") + 
  guides(fill = "none", colour = "none") +
  NULL
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-rep-returns-dist-facet-treat.png")), width = 7.5, height = 5.0, dpi = 500)
