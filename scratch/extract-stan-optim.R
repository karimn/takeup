library(tidyverse)
library(tidybayes)
library(testthat)
library(posterior)
library(cmdstanr)
library(sf)
library(nleqslv)


from_cmdstan_fit = TRUE
load_exposed_functions = FALSE

if (load_exposed_functions) {
    source("scratch/expose-stan-functions.R")
}

source("dist_structural_util.R")
source(file.path("optim", "optim-functions.R"))

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
  

distance_data = read_rds("optim/data/full-experiment.rds")

dist_df = distance_data$long_distance_mat %>%
    filter(dist < 10000) %>%
    mutate(dist_index = 1:n())


if (from_cmdstan_fit) {
    fit_files = str_glue("data/stan_analysis_data/dist_fit86_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-1.csv")

    fit = as_cmdstan_fit(fit_files)

    fit

    library(tidybayes)


    dist_beta_draws = gather_draws(
        fit, 
        dist_beta_v[j], 
        centered_cluster_dist_beta_1ord[j, k],
    )

    rm(fit)
    gc()

    dist_beta_draws %>%
        filter(j == 1 | is.na(j)) %>%
        summarise_draws()

    rm(dist_beta_draws)

    dist_beta_draws %>%
        summarise_draws()


    optim_fit_draws = fit$draws(
        variables = c(
            # "optim_w", 
            # "total_error_sd",
            # "net_private_benefit",
            # "visibility", 
            # "structural_cluster_benefit",
            # "obs_cluster_mu_rep",
            # "sim_delta", 
            "dist_beta"
    
        )
    )
    optim_fit_draws %>%
        saveRDS(
            "temp-data/optim-fit-draws.rds"
        )
    rm(fit)
    gc()
} else {
    optim_fit_draws = read_rds("temp-data/optim-fit-draws.rds")
}
optim_fit_draws[, , "structural_cluster_benefit[3,1]"]
optim_fit_draws[, , "visibility[1,1,1]"]
optim_fit_draws[, , "obs_cluster_mu_rep[1]"]
optim_fit_draws[, , "sim_delta[1]"]

stop()

posterior_param_draws = read_csv(
    file.path(
        "data/stan_analysis_data", 
        str_interp(
            "param_posterior_draws_dist_fit71_STRUCTURAL_LINEAR_U_SHOCKS.csv"
        )
    )
)


posterior_R_predictions = read_rds(
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


dist_fit_data <- str_split("71", ",") %>% 
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

stan_fit_df = dist_fit_data %>%
    filter(model == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
    filter(fct_match(fit_type, "fit"))

stan_fit_draws = extract_and_harmonise_variables(
    stan_fit_df
)



stan_fit_draws


optim_rvars = spread_rvars(
    optim_fit_draws, 
    # optim_w[b_index, mu_index, dist_index], 
    net_private_benefit[b_index, mu_index, dist_index], 
    visibility[b_index, mu_index, dist_index],
    total_error_sd
    # ndraws = 100
)





close_dist_indices = dist_df %>%
    filter(dist < 5000) %>%
    pull(dist_index)

optim_dist_df = dist_df %>%
    left_join(
        optim_rvars %>%
            filter(dist_index %in% close_dist_indices)
    ) %>%
        filter(dist < 5000)


optim_rvars

long_optim_dist_df = optim_dist_df %>%
    unnest_rvars()


long_optim_dist_df  = long_optim_dist_df %>%
    mutate(
        pred_takeup = 1 - pnorm(optim_w / `total_error_sd[1,1]`)
    ) %>%
    ungroup() %>%
    mutate(
        u_sd = sqrt(`total_error_sd[1,1]`^2 - 1)
    )







long_complete_optim_dist_df = long_optim_dist_df %>%
    filter(!is.nan(optim_w)) %>%
    head(1000) 




library(furrr)
plan(multicore, workers = 12)
long_complete_optim_dist_df = long_complete_optim_dist_df %>%
  mutate(
      R_v_star = future_pmap(
        list(
          dist, 
          net_private_benefit,
          visibility, 
          `total_error_sd[1,1]`, 
          u_sd
        ),
        ~find_v_star(
          distance = ..1, 
          b = ..2, 
          mu_rep = ..3, 
          total_error_sd = ..4, 
          u_sd = ..5, 
          bounds = c(-Inf, Inf)
      ),
      .progress = TRUE
  )
)
long_complete_optim_dist_df = long_complete_optim_dist_df %>%
    unnest_wider(R_v_star)


n_vstar_disagree = long_complete_optim_dist_df %>%
    mutate(
        diff = abs(optim_w - v_star)
    ) %>%
    filter(diff > 0.01) %>%
    nrow()

test_that("V star matches on same values", {
    expect_equal(n_vstar_disagree, 0)
})


optim_rvars

long_optim_dist_df %>%
    filter(.draw == 1) %>%
    tail()

treatments = c(
    "control",
    "ink", 
    "calendar", 
    "bracelet"
)

long_optim_dist_df = long_optim_dist_df %>%
    mutate(
        b_treatment = treatments[b_index], 
        mu_treatment = treatments[mu_index]
    )

stop()




summ_optim_df = long_optim_dist_df %>%
    filter(!is.nan(optim_w)) %>%
    group_by(
        b_treatment, 
        mu_treatment, 
        dist_index
    ) %>%
    summarise(
        mean_est = mean(pred_takeup), 
        median_est = median(pred_takeup), 
        dist = unique(dist)
    )


long_var_df = long_optim_dist_df %>%
    filter(!is.na(optim_w)) %>%
    select(
        b_treatment, 
        mu_treatment, 
        dist, 
        optim_w, 
        pred_takeup, 
        net_private_benefit, 
        visibility
    ) %>%
    gather(
        variable, value, 
        optim_w:visibility
    ) %>%
    group_by(
        b_treatment, mu_treatment, dist, variable
    ) %>%
    summarise(
        mean_est = mean(value), 
        median_est = median(value)
    )




cluster_id_dist_crosswalk = analysis_data %>%
    select(
        cluster_id,
        treatment = assigned.treatment, 
        dist = cluster.dist.to.pot
        ) %>%
        unique()


long_optim_dist_df

stan_fit_df %>%
    select(obs_cluster_mu_rep) %>%
    unnest()  %>%
    unnest(iter_data) 


comp_df = bind_rows(
    
    
    stan_fit_df %>%
    select(obs_cluster_mu_rep) %>%
    unnest()  %>%
    unnest(iter_data) %>%
    group_by(variable) %>%
    summarise(
        mean_est = mean(iter_est)
    ) %>%
    mutate(
        cluster_id = str_extract(variable, "\\d+") %>% as.numeric
    )  %>%
    left_join(
        cluster_id_dist_crosswalk
    )  %>% mutate(type = "stanfit"),
long_var_df %>%
    filter(b_treatment == mu_treatment) %>%
    filter(variable == "visibility") %>%
    mutate(type = "gen_quantities") %>%
    mutate(treatment = b_treatment) 
)


comp_df %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = type
    )) +
    geom_point() +
    facet_wrap(~treatment)



long_var_df %>%
    filter(b_treatment == mu_treatment) %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = b_treatment
    )) +
    geom_point() +
    facet_wrap(
        ~variable, 
        scales = "free"
        )


summ_optim_df %>%
    write_csv("temp-data/stan-optim-quant-draws.csv")

summ_optim_df %>%
    filter(mu_treatment == "control") %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = b_treatment 
    )) +
    geom_point() 

summ_optim_df %>%
    filter(b_treatment == "control") %>%
    filter(dist < 3500) %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = mu_treatment 
    )) +
    geom_point() 


summ_optim_df %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = b_treatment
    )) +
    geom_point() +
    facet_wrap(
        ~mu_treatment
    )

optim_distance_rvars = left_join(
    optim_rvars, 
    dist_df, 
    by = "dist_index"
)


optim_distance_rvars
stop()


optim_rvars

optim_rvars %>%
    filter(!is.nan(`total_error_sd[, 1]`))

long_optim_rvars = optim_rvars %>%
    unnest_rvars()

long_optim_rvars

optim_draws_df = as_draws_df(
    optim_fit_draws
)

optim_draws_df[, 1] %>%
    filter(!is.nan(`optim_w[1,1,1]`))

optim_fit_draws
# cf_fit_draws = fit$draws(
#     variables = c("cluster_cf_cutoff", "beta", "structural_cluster_benefit_cost")
# )


beta_rvar = gather_rvars(
    cf_fit_draws,
    structural_cluster_benefit_cost
)

beta_rvar

rm(fit)
gc()
stop()




ed = spread_rvars(
    cf_fit_draws,
    cluster_cf_cutoff[b_index, mu_index, cluster_index]
)
ed

dist_ed = ed %>%
  left_join(
    analysis_data %>%
        select(
            cluster_id,
            dist = cluster.dist.to.pot) %>%
        unique(), 
        by = c("cluster_index" = "cluster_id")
  ) 
dist_ed
summ_dist_ed = dist_ed %>%
    group_by(b_index, mu_index, cluster_index) %>%
    summarise(value = mean(cluster_cf_cutoff, na.rm = TRUE), dist = unique(dist)) 
summ_dist_ed
summ_dist_ed %>%
    tail()



summ_dist_ed %>%
    filter(b_index == 1 | b_index == 5) %>%
    filter((b_index == 1 & dist <= 1250) | (b_index == 5 & dist > 1250)) %>%
    ggplot(aes( 
        x = dist, 
        y = value, 
        colour = factor(mu_index)
    )) +
    geom_point()

cf_fit_draws

summ_ed_df = ed %>%
    group_by(i, j, k) %>%
    summarise(ed = mean(cluster_w_cutoff))


treatments = c(
    "control",
    "ink",
    "calendar",
    "bracelet"
)
dist_seq = seq(0, 5000, 100)

summ_ed_df = summ_ed_df %>%
    mutate(
        treatment = factor(treatments[k], levels = treatments), 
        dist = dist_seq[i]
        )
summ_ed_df %>%
    ggplot(aes(
        x = dist, 
        y = ed, 
        colour = treatment
    )) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = 3000, type = "longdash")

summ_ed_df %>%
    filter(treatment == "bracelet" ) %>%
    ggplot(aes(
        x = dist, 
        y = ed, 
        colour = treatment
    )) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = 3000, type = "longdash") 


stop()
fit_rvars

fit_draws

fit_rvars %>%
    filter(!is.nan(cluster_w_cutoff))

bc_draws = fit$draws(
  variables = c(
    "structural_cluster_benefit_cost", 
    "cluster_dist_cost", 
    "obs_cluster_mu_rep",
    "total_error_sd[1]", 
    "u_sd[1]", 
    "beta", 
    "dist_beta_v")
)


betas = spread_rvars(
    bc_draws, 
    beta[k], 
    dist_beta_v[k]
    )