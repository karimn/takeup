library(tidyverse)
library(tidybayes)
library(posterior)
library(cmdstanr)
library(sf)

fit_file = "data/stan_analysis_data/dist_fit75_STRUCTURAL_LINEAR_U_SHOCKS-1.csv"

fit = as_cmdstan_fit(fit_file)


distance_data = read_rds("optim/data/full-experiment.rds")

dist_df = distance_data$long_distance_mat %>%
    filter(dist < 10000) %>%
    mutate(dist_index = 1:n())

# To find fixed point need:
# benefit_cost
# mu_rep
# total_error_sd
# u_sd


# fit_draws = fit$draws(
#     variables = c(
#         "cluster_w_cutoff",
#         "cluster_rep_return"
#     )
# )


optim_fit_draws = fit$draws(
    variables = c(
        "optim_w", 
        "total_error_sd"
    )
)




rm(fit)
gc()

optim_rvars = spread_rvars(
    optim_fit_draws, 
    optim_w[b_index, mu_index, dist_index], 
    total_error_sd
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

long_optim_dist_df = optim_dist_df %>%
    unnest_rvars()


long_optim_dist_df  = long_optim_dist_df %>%
    mutate(
        pred_takeup = 1 - pnorm(optim_w / `total_error_sd[1,1]`)
    )




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