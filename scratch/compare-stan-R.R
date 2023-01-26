
library(tidyverse)
library(sf)
library(broom)
library(posterior)
library(tidyverse)
library(tidybayes)
library(broom)
library(rstan)
library(sf)
library(nleqslv)
library(cmdstanr)
library(econometr)
library(furrr)

source(file.path("optim", "optim-functions.R"))
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("stan_models", "stan_owen_t.R"))

stan_owen_t = exposed_stan_func$stan_owen_t

set.seed(19484)


theme_set(theme_minimal())

struct_param_draws = read_csv(
    file.path(
        "data/stan_analysis_data", 
        str_interp(
            "param_posterior_draws_dist_fit71_STRUCTURAL_LINEAR_U_SHOCKS.csv"
        )
    )
)







max_draw = max(struct_param_draws$.draw)
draw_treat_grid = expand.grid(
    draw = sample(1:max_draw, size = 100),
    treatment = c("control", "calendar", "bracelet", "ink") 
)


dist_data = read_rds(
    file.path(
        "optim/data/",
        "full-experiment.rds"
    )
)

sd_of_dist = dist_data$sd_of_dist

pred_functions = map2(
    draw_treat_grid$draw,
    draw_treat_grid$treatment,
    ~extract_params(
        param_draws = struct_param_draws,
        private_benefit_treatment = .y,
        visibility_treatment = .y,
        draw_id = .x,
        dist_sd = sd_of_dist,
        j_id = 1,
        rep_cutoff = Inf,
        dist_cutoff = Inf, 
        bounds = c(-Inf, Inf),
        mu_rep_type = 0
    ) %>% find_pred_takeup()
)


library(furrr)
plan(multicore, workers = 12)
dist_df = tibble(
    dist = seq(from = 0, to = 3500, by = 100)
)
pred_dist_df = future_imap_dfr(
    pred_functions,
    ~{
        dist_df %>%
            mutate( 
                as_tibble(.x(dist)),
                draw = draw_treat_grid[.y, "draw"],
                treatment = draw_treat_grid[.y, "treatment"]
            )
    },
    .options = furrr_options(
        seed = TRUE
    ),
    .progress = TRUE
)



data_to_save = lst(
    param_grid = struct_param_draws,
    dist_sd = sd_of_dist, 
    pred_df = pred_dist_df
)


saveRDS(
    data_to_save,
    "temp-data/r-predicted-data.rds"
)

extract_and_harmonise_variables = function(stan_df){
    v_star_df = stan_df %>%
        select(cluster_w_cutoff) %>%
        unnest(cols = c(cluster_w_cutoff))  %>%
        select(roc_distance, assigned_treatment, mean_est, per_0.5) %>%
        rename(
            dist = roc_distance,
            treatment = assigned_treatment, 
            mean_est = mean_est, 
            median_est = per_0.5
        ) %>%
        mutate(variable = "v_star")

    pred_takeup_df = stan_df %>%
        select(obs_cluster_takeup_level) %>%
        unnest(cols = c(obs_cluster_takeup_level)) %>%
        select(
            dist = assigned_dist_obs, 
            treatment = assigned_treatment_obs, 
            median_est = per_0.5
        ) %>% 
        mutate(variable = "actual_takeup")


    rep_return_df = stan_df %>%
        select(cluster_rep_return) %>%
        unnest(cols = c(cluster_rep_return)) %>%
        select( 
            dist = roc_distance, 
            mean_est = mean_est, 
            median_est = per_0.5, 
            treatment = assigned_treatment
        ) %>% 
        mutate(variable = "rep_return")

    prop_df = stan_df %>%
        select(cluster_takeup_prop) %>%
        unnest(c(cluster_takeup_prop)) %>%
        select(
            dist = roc_distance, 
            mean_est = mean_est, 
            median_est = per_0.5, 
            treatment = assigned_treatment
        ) %>%
        mutate(variable = "pred_takeup")

    df = bind_rows(
        v_star_df, 
        pred_takeup_df, 
        rep_return_df, 
        prop_df
    )
    return(df)

}


params = list(structural_takeup_version = 71)
dist_fit_data <- str_split(params$structural_takeup_version, ",") %>% 
  pluck(1) %>% 
  as.integer() %>% 
  map_dfr(~ {
    temp_env <- new.env()
    load(file.path("temp-data", str_glue("processed_dist_fit{.x}.RData")), envir = temp_env)
    
    temp_env$dist_fit_data %>% 
      filter(fct_match(model_type, "structural")) %>% 
      mutate(version = .x)
  }) %>% 
  group_by(model, fit_type, model_type) %>% 
  filter(min_rank(version) == n()) %>% 
  ungroup()



# fit_version = 71
# load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))

stan_fit_df = dist_fit_data %>%
    filter(model == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
    filter(fct_match(fit_type, "fit"))

stan_clean_df = extract_and_harmonise_variables(
    stan_fit_df
)


load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# stick to monitored sms.treatment group
# remove sms.treatment.2
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data %>% 
mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group,
    sms_treatment = factor(sms.treatment.2))

stan_mu_df = dist_fit_data %>% 
  filter(fct_match(model_type, "structural"))  %>%
  filter(fit_type == "fit") %>%
  filter(model == 'STRUCTURAL_LINEAR_U_SHOCKS') %>%
  select(fit_type, model, obs_cluster_mu_rep)   %>%
  unnest(obs_cluster_mu_rep) %>%
  unnest(iter_data)  %>%
  group_by(variable, model) %>%
  summarise(
    median = median(iter_est), 
    mean = mean(iter_est), 
    conf.low = quantile(iter_est, 0.25), 
    conf.high = quantile(iter_est, 0.75)) %>%
  mutate(cluster_id = str_extract(variable, "\\d+") %>% as.numeric) %>%
  arrange(cluster_id) %>%
  left_join(
    analysis_data %>%
        select(cluster_id,assigned_dist_group, assigned.treatment, cluster.dist.to.pot), 
        by = "cluster_id"
  ) 

summ_stan_mu_df = stan_mu_df %>%
    ungroup() %>%
    select(
        dist = cluster.dist.to.pot, 
        treatment = assigned.treatment, 
        median_est = median, 
        mean_est = mean
    ) %>%
    unique() 

summ_df = pred_dist_df %>%
    gather(variable, value, pred_takeup:total_error_sd) %>%
    group_by(
        dist, 
        treatment, 
        variable
    ) %>%
    summarise(
        mean_est = mean(value), 
        median_est = median(value)
    ) %>%
    ungroup() 


comp_df = bind_rows(
    summ_df %>% mutate(type = "ed"), 
    stan_clean_df %>% mutate(type = 'stan') %>%
        filter(variable != "mu_rep"),
    summ_stan_mu_df %>% 
        mutate(type = "stan") %>% 
        mutate(variable = "mu_rep")
)


comp_df %>% 
    filter(dist < 3500) %>%
    filter(variable %in% c("pred_takeup", "v_star", "mu_rep")) %>%
    ggplot(aes(
        x = dist, 
        y = median_est, 
        colour = type
    )) +
    geom_point(size = 2) +
    facet_grid(
        variable~treatment, 
        scales = "free") +
    theme_bw()


ggsave("temp-plots/compare-stan-R-panel.png", width = 10, height = 10, dpi = 500)
