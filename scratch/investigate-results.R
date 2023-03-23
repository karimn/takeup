library(magrittr)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggmap)
library(ggstance)
library(gridExtra)
library(cowplot)
library(sf)
library(knitr)
library(modelr)
library(car)
library(rstan)
library(latex2exp)
library(ggthemes)
library(cmdstanr)
library(posterior)
library(tidybayes)

library(econometr)

load("temp-data/processed_dist_fit80_lite.RData")



grab_cal_draws = function(x){
  fit_obj = as_cmdstan_fit(
    str_glue("data/stan_analysis_data/dist_fit{x}_STRUCTURAL_LINEAR_U_SHOCKS-1.csv"))

  prefer_cal_draws = gather_rvars(
    fit_obj,
    prob_prefer_calendar[index]
  )
  return(prefer_cal_draws)
}

fit_obj = as_cmdstan_fit(
  str_glue("data/stan_analysis_data/dist_fit80_STRUCTURAL_LINEAR_U_SHOCKS-1.csv"))

fit_obj 

p_hat_draws = gather_rvars(
  fit_obj,
  base_mu_rep, 
  centered_cluster_beta_1ord[j,k],
  centered_cluster_dist_beta_1ord[j,k]
)

rm(fit_obj)
gc()

subset_p_hat_draws = p_hat_draws %>%
  filter(is.na(j) | j == 1)


treat_levels = c("control", "ink", "calendar", "bracelet")

wide_p_hat_draws = subset_p_hat_draws %>%
  select(-j)  %>%
  mutate(base_mu_rep = .value[.variable == "base_mu_rep"]) %>%
  mutate(control_beta = .value[.variable == "centered_cluster_beta_1ord" & k == 1]) %>% 
  # mutate(control_dist_beta = .value[.variable == "centered_cluster_dist_beta_1ord" & k == 1]) %>% 
  filter(.variable != "base_mu_rep") %>%
  pivot_wider(
    names_from = .variable, values_from = .value
  )  %>%
  mutate(
    treatment = factor(treat_levels, levels = treat_levels)
  )

inv_logit = function(x){1/(1+exp(-x))}

calculate_p_hat = function(wide_data, dist) {
  wide_data %>%
    mutate(not_control = k != 1) %>%
    mutate(dist = dist) %>%
    mutate(
      p_hat = inv_logit(not_control*control_beta + centered_cluster_beta_1ord + centered_cluster_dist_beta_1ord*dist)
      ) %>%
    select(treatment, dist, p_hat, everything())
}


sd_of_dist = read_rds("temp-data/sd_of_dist.rds")


distances = seq(from = 0, to = 2500, length.out = 20) / sd_of_dist


p_hat_df = map_dfr(
  distances,
  ~calculate_p_hat(
    wide_p_hat_draws, 
    .x
  )
)



p_hat_df %>%
  mutate(
    p_hat = mean(p_hat), 
    distance = dist*sd_of_dist
  ) %>%
  ggplot(aes(
    x = distance, 
    y = p_hat, 
    colour = treatment
  )) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Distance", 
    y = "P hat", 
    title = "Estimated 1st Order Belief Proportion vs Distance", 
    colour = ""
  ) +
  ggthemes::scale_color_canva(
    palette = "Primary colors with a vibrant twist", 
    labels = as_labeller(str_to_title)
  )

ggsave("temp-plots/anne-check/pr-1ord-continuous-distance.png", 
width = 8, 
height = 6, 
dpi = 500)


# prefer_cal_draws_77 = grab_cal_draws(77)
# prefer_cal_draws_75 = grab_cal_draws(75)
# prefer_cal_draws_71 = grab_cal_draws(71)
prefer_cal_draws = grab_cal_draws(80)

prefer_cal_draws %>%
  slice(1) %>%
  unnest_rvars() %>%
  print(n = 50)

stop()


stop()


if (interactive()) {
  params = lst(
    fit_version = 75,
    input_path = "data/stan_analysis_data",
    output_path = "temp-data", 
    models = "STRUCTURAL_LINEAR_U_SHOCKS"
  )
}

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path("dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))



fit_version <- params$fit_version



default_top_levels = c("Bracelet", "Combined")

# 66 ed fit
# 60 Karim fit
# 62 also Karim fit


model_fit_by = if_else(fit_version %in% c(60, 62), "Karim", "Ed")


models_we_want = c(
  params$models
)


quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)


output_basepath = str_glue("temp-data/output_dist_fit{fit_version}")
dir.create(output_basepath, showWarnings = FALSE)
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
load(file.path("temp-data", str_interp("processed_dist_fit71_lite.RData")))

dist_fit_data %<>%
  left_join(tribble(
    ~ fit_type,        ~ model_color,
      "fit",           "black", 
      "prior-predict", "darkgrey",
  ), by = "fit_type")


dist_fit_data = dist_fit_data %>%
  filter(model %in% c("REDUCED_FORM_NO_RESTRICT", "STRUCTURAL_LINEAR_U_SHOCKS"))

rf_analysis_data <- dist_fit_data %>% 
  filter(
    fct_match(model_type, "reduced form"),
    fct_match(fit_type, "fit"),
  ) %$% 
  stan_data[[1]]$analysis_data 
delta <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))

stan_data = (dist_fit_data %>%
  filter(fct_match(model_type, "structural")) %>%
  select(stan_data) %>%
  pull())[[1]]




belief_data = dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  pull(beliefs_results) %>%
  first() 



dist_fit_data

#### Separate Beliefs ####
belief_te_1ord = belief_data$ate_knows %>%
  filter(assigned_treatment_left != "control") %>%
  mutate(assigned_treatment_left = fct_drop(assigned_treatment_left)) %>%
  mutate(assigned_treatment_left = fct_relabel(assigned_treatment_left, str_to_title)) %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 1,
    crossbar_width = 0.4) +
  labs(title = "") +
  theme_minimal() + 
  theme(legend.position = "bottom") + 
    NULL


iter_data = (dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  pull(beliefs_results) %>%
  first())$prob_knows %>%
  head(1) %>%
  unnest(iter_data)




full_data_env = new.env()

with_env = function(f, e = parent.frame()) {
    stopifnot(is.function(f))
    environment(f) = e
    f
}

load_full_data_function = function(){
    load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}.RData")))
    return(dist_fit_data)
}


full_dist_fit_data = with_env(load_full_data_function, full_data_env)() %>%
  filter(model %in% c("REDUCED_FORM_NO_RESTRICT", "STRUCTURAL_LINEAR_U_SHOCKS"))

full_iter_data = (full_dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  pull(beliefs_results) %>%
  first())$prob_knows %>%
  head(1) %>%
  unnest(iter_data)


full_iter_data %>%
    ggplot() +
    geom_histogram(
        aes(x = iter_est)) 
    


fit_obj_75 = as_cmdstan_fit("data/stan_analysis_data/dist_fit75_STRUCTURAL_LINEAR_U_SHOCKS-1.csv")
fit_obj_71 = as_cmdstan_fit("data/stan_analysis_data/dist_fit71_STRUCTURAL_LINEAR_U_SHOCKS-1.csv")


cf_71 = gather_rvars(fit_obj_71, cluster_w_cutoff[a,b,c])
cf_75 = gather_rvars(fit_obj_75, cluster_w_cutoff[a,b,c])


cf_71 %>%
  head(1) %>%
  unnest_rvars()

cf_75 %>%
  head(1) %>%
  unnest_rvars()

full_iter_data %>%
    filter(is.na(iter_est))

full_iter_data %>%
    filter(iter_est == 0)

iter_data  %>%
  filter(is.na(iter_est))
