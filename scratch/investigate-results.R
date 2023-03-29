library(magrittr)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggstance)
library(gridExtra)
library(cowplot)
library(sf)
library(knitr)
library(modelr)
library(car)
library(rstan)
library(ggthemes)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(furrr)


models_we_want = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"


treat_levels = c("control", "ink", "calendar", "bracelet")
inv_logit = function(x){1/(1+exp(-x))}

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


analysis_data <- monitored_nosms_data %>% 
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group)

cluster_treatment_map = distinct(analysis_data, assigned_treatment, assigned_dist_group) %>% 
  arrange(assigned_dist_group, assigned_treatment) # We must arrange by distance first

analysis_data %<>% 
  nest_join(
    endline.know.table.data %>% 
      filter(fct_match(know.table.type, "table.A")),
    by = "KEY.individ", 
    name = "knowledge_data"
  ) %>% 
  mutate(
    map_dfr(knowledge_data, ~ {
      tibble(
        obs_know_person = sum(.x$num.recognized),
        obs_know_person_prop = mean(.x$num.recognized),
        knows_other_dewormed = sum(fct_match(.x$dewormed, c("yes", "no")), na.rm = TRUE),
        knows_other_dewormed_yes = sum(fct_match(.x$dewormed, "yes"), na.rm = TRUE),
        thinks_other_knows = sum(fct_match(.x$second.order, c("yes", "no")), na.rm = TRUE),
        thinks_other_knows_yes = sum(fct_match(.x$second.order, "yes"), na.rm = TRUE),
      )
    }
  ))

beliefs_ate_pairs <- cluster_treatment_map %>% 
  mutate(treatment_id = seq(n())) %>% {
      bind_rows(
        left_join(., filter(., fct_match(assigned_treatment, "control")), by = c("assigned_dist_group"), suffix = c("", "_control")) %>% 
          filter(assigned_treatment != assigned_treatment_control) %>% 
          select(treatment_id, treatment_id_control),
        
        left_join(., filter(., fct_match(assigned_dist_group, "close")), by = c("assigned_treatment"), suffix = c("", "_control")) %>% 
          filter(assigned_dist_group != assigned_dist_group_control) %>% 
          select(treatment_id, treatment_id_control),
      )
} %>%
  arrange(treatment_id, treatment_id_control) 



calculate_belief_latent_predictor = function(beta, 
                                             dist_beta, 
                                             dist, 
                                             control_beta, 
                                             control_dist_beta, 
                                             control) {
    # Don't want to double count control and add it twice
    if (control == FALSE) {
        val = beta + dist*dist_beta + control_beta + control_dist_beta*dist
    } else {
        val = beta + dist*dist_beta 
    }
    return(val)
}

belief_data = analysis_data %>%
  filter(obs_know_person > 0)  %>%
  mutate(
    j = 1:n()
  )


belief_data %>%
  group_by(
    assigned_treatment, 
    assigned_dist_group
  ) %>%
  select(obs_know_person, knows_other_dewormed, assigned_treatment, assigned_dist_group) %>%
  group_by(
    assigned_treatment, 
    assigned_dist_group
  ) %>%
  summarise(
    pr = mean(knows_other_dewormed/obs_know_person)
  )  %>%
  spread( 
    assigned_treatment, 
    pr
  )




load_param_draws = function(fit_version, model, chain) {
  fit_obj = as_cmdstan_fit(
    str_glue("data/stan_analysis_data/dist_fit{fit_version}_{model}-{chain}.csv"))
  draws = spread_rvars(
    fit_obj,
    centered_obs_beta_1ord[j, k],
    centered_obs_dist_beta_1ord[j, k]
  ) %>%
    mutate(chain = chain, model = model, fit_version = fit_version)
  return(draws)
}

load_p_draws = function(fit_version, model, chain) {
  fit_obj = as_cmdstan_fit(
    str_glue("data/stan_analysis_data/dist_fit{fit_version}_{model}-{chain}.csv"))
  draws = spread_rvars(
    fit_obj,
    obs_lin_pred_1ord[j, k],
    obs_prob_1ord[j, k]
  ) %>%
    mutate(chain = chain, model = model, fit_version = fit_version)
  return(draws)
}


belief_data = analysis_data %>%
  filter(obs_know_person > 0)  %>%
  mutate(
    j = 1:n()
  )

raw_ndgt_model_draws = map_dfr(
  1, 
  ~load_p_draws(
    fit_version = 86, 
    model = models_we_want,  
    chain = .x
  )
)

raw_ndgt_model_draws = raw_ndgt_model_draws %>%
  group_by(
    j, k, model, fit_version
  ) %>%
  summarise(
    across(where(is_rvar), ~rvar(c(draws_of(.x))))
  ) %>%
  ungroup()

k_cw_df = tibble(
  k = 1:4, 
  treatment = c("control", "ink", "calendar", "bracelet")
)
k_dist_cw_df = tibble(
  k = 1:8, 
  treatment = rep(treat_levels, 2)
)

ndgt_model_draws = raw_ndgt_model_draws %>%
  left_join(
    belief_data %>%
      select(j, assigned_dist_group), 
    by = "j"
  )



subset_ndgt_model_draws = ndgt_model_draws %>%
  mutate(
    stan_dist = if_else(k <= 4, "close", "far")
  ) %>%
  filter(stan_dist == assigned_dist_group) %>%
  left_join(k_dist_cw_df)

model_draws = map_dfr(
  1,
  ~load_param_draws(
    fit_version = 86,
    model = models_we_want,
    chain = .x
  )
)

# pool over chains


model_draws = model_draws %>%
  group_by(
    j, k, model, fit_version
  ) %>%
  summarise(
    across(where(is_rvar), ~rvar(c(draws_of(.x))))
  ) %>%
  ungroup()






dist_model_draws = model_draws %>%
  left_join(
    belief_data %>%
      select(assigned_dist_group, j, standard_cluster.dist.to.pot)
    ) %>%
  mutate(
    not_control = k != 1, 
    not_close = assigned_dist_group != "close"
  ) %>%
  group_by(j) %>%
  mutate( 
    control_beta = centered_obs_beta_1ord[k == 1],
    control_dist_beta = centered_obs_dist_beta_1ord[k == 1]
  ) %>%
  ungroup() %>%
  mutate(
    ed_pred = centered_obs_beta_1ord + 
              centered_obs_dist_beta_1ord*standard_cluster.dist.to.pot + 
              not_control*(control_beta + control_dist_beta*standard_cluster.dist.to.pot)
              # not_close*(centered_obs_beta_1ord)
  )  



joint_dist_model_draws = dist_model_draws %>%
  left_join(k_cw_df) %>%
  left_join(
    subset_ndgt_model_draws %>% select(j, treatment, contains("1ord")), 
    by = c("j", "treatment")
  )

joint_dist_model_draws %>%
  select(obs_lin_pred_1ord)




far_dist_model_draws = joint_dist_model_draws %>%
  filter(assigned_dist_group == "far") %>%
  select( 
    j,
    k, 
    obs_lin_pred_1ord, 
    ed_pred,
    standard_cluster.dist.to.pot
  )

joint_dist_model_draws %>%
  # filter(assigned_dist_group == "far") %>%
  mutate(
    difference = obs_lin_pred_1ord - ed_pred, 
    d_s = difference/standard_cluster.dist.to.pot
  )  %>%
  group_by(k, assigned_dist_group) %>%
  slice(1:5) %>%
  select(k, d_s) %>%
  print(n = 50)




summ_dist_model_draws = joint_dist_model_draws %>%
  # filter(assigned_dist_group == "close") %>%
  group_by(
    treatment, assigned_dist_group
  ) %>%
  summarise(
    ed_p_hat = mean(rvar_mean(inv_logit(ed_pred))),
    stan_pred_p_hat = mean(rvar_mean(inv_logit(obs_lin_pred_1ord))),
    stan_p_hat = mean(rvar_mean(obs_prob_1ord))
  ) %>%
  arrange(assigned_dist_group) %>%
  left_join(k_cw_df) 

summ_draw_table = summ_dist_model_draws %>%
  gather(
    variable, value, contains("p_hat")
  ) %>%
  mutate(treatment = factor(treatment, levels = treat_levels)) %>%
  pivot_wider(
    id_cols = c(variable, assigned_dist_group), 
    names_from = treatment, 
    values_from = value
  ) %>%
  arrange(assigned_dist_group, variable) %>%
  select(variable, assigned_dist_group, all_of(treat_levels))


raw_summ = belief_data %>%
  group_by(
    assigned_treatment, 
    assigned_dist_group
  ) %>%
  select(obs_know_person, knows_other_dewormed, assigned_treatment, assigned_dist_group) %>%
  group_by(
    assigned_treatment, 
    assigned_dist_group
  ) %>%
  summarise(
    pr = mean(knows_other_dewormed/obs_know_person)
  )  %>%
  mutate(assigned_treatment = factor(assigned_treatment, levels = treat_levels)) %>%
  spread( 
    assigned_treatment, 
    pr
  )

full_stan_summ =  ndgt_model_draws %>%
  group_by(k) %>%
  summarise(
    pr = rvar_mean(obs_prob_1ord)
  ) %>%
  mutate(pr = mean(pr)) %>%
  left_join(
    k_dist_cw_df
  ) %>%
  mutate(assigned_dist_group = if_else(
    k <= 4, "close", "far"
  )) %>%
  mutate(treatment = factor(treatment, levels = treat_levels)) %>%
  select(-k) %>%
  spread(
    treatment, pr
  )


raw_summ
summ_draw_table
full_stan_summ






dist_summ_df = joint_dist_model_draws %>%
  # filter(assigned_dist_group == "close") %>%
  group_by(
    treatment, standard_cluster.dist.to.pot
  ) %>%
  summarise(
    ed_p_hat = mean(rvar_mean(inv_logit(ed_pred))),
    stan_pred_p_hat = mean(rvar_mean(inv_logit(obs_lin_pred_1ord))),
    stan_p_hat = mean(rvar_mean(obs_prob_1ord))
  ) %>%
  arrange(standard_cluster.dist.to.pot) %>%
  left_join(k_cw_df) 


unique_dists = unique(dist_summ_df$standard_cluster.dist.to.pot)


dist_summ_df %>%
  filter(standard_cluster.dist.to.pot == unique_dists[[3]])

dist_summ_df %>%
  select(-k) %>%
  gather(variable, value, -standard_cluster.dist.to.pot, -treatment )  %>%
  filter(variable == "ed_p_hat") %>%
  ggplot(aes(
     x = standard_cluster.dist.to.pot*630,
     y = value, 
     colour = treatment
  )) +
  facet_wrap(~variable, ncol = 1) +
  geom_line(size=1)
  # geom_smooth() 

  ggsave("temp-plots/pr-by-treat-dist.png", width = 8, height = 6, dpi = 500)






rf_glm_fit = belief_data %>%
  select(obs_know_person, knows_other_dewormed, assigned_treatment, assigned_dist_group) %>%
  glm(
    data = .,
    knows_other_dewormed/obs_know_person ~ 0 + assigned_treatment:assigned_dist_group,
    weights = obs_know_person, 
    family = binomial(link = "logit")
  )


rf_cts_glm_fit = belief_data %>%
  select(
    obs_know_person, 
    knows_other_dewormed, 
    assigned_treatment, 
    assigned_dist_group, 
    standard_cluster.dist.to.pot, 
    cluster.dist.to.pot) %>%
  glm(
    data = .,
    knows_other_dewormed/obs_know_person ~ 0 + assigned_treatment*cluster.dist.to.pot,
    weights = obs_know_person, 
    family = binomial(link = "logit")
  )



treatment = c(
  "bracelet", 
  "calendar", 
  "control", 
  "ink"
)


discrete_pred_hat = map_dfr(
  treatment,
  ~predictions(
    rf_glm_fit, 
    newdata = belief_data %>%
    select(
      obs_know_person, 
      knows_other_dewormed, 
      assigned_treatment, 
      assigned_dist_group, 
      standard_cluster.dist.to.pot, 
      cluster.dist.to.pot)  %>%
    mutate(
      assigned_treatment = .x
    )) %>% as_tibble()
) 

cts_pred_hat = map_dfr(
  treatment,
  ~predictions(
    rf_cts_glm_fit, 
    newdata = belief_data %>%
    select(
      obs_know_person, 
      knows_other_dewormed, 
      assigned_treatment, 
      assigned_dist_group, 
      standard_cluster.dist.to.pot, 
      cluster.dist.to.pot)  %>%
    mutate(
      assigned_treatment = .x
    )) %>% as_tibble()
) 


levels_discrete_pred_hat = discrete_pred_hat %>%
  group_by(
    assigned_dist_group, 
    assigned_treatment
  ) %>%
  summarise(
    p_hat = median(estimate)
  ) %>%
  spread(assigned_treatment, p_hat) 

levels_cts_pred_hat = cts_pred_hat %>%
  group_by(
    assigned_dist_group, 
    assigned_treatment
  ) %>%
  summarise(
    p_hat = median(estimate)
  ) %>%
  spread(assigned_treatment, p_hat) 

cts_rough_te = levels_cts_pred_hat %>%
  mutate(across(where(is.numeric), ~.x - control))

discrete_rough_te = levels_discrete_pred_hat %>%
  mutate(across(where(is.numeric), ~.x - control))


cts_rough_te
discrete_rough_te

cts_pred_hat %>%
  group_by(
    assigned_dist_group, 
    assigned_treatment
  ) %>%
  summarise(
    p_hat = mean(estimate)
  ) %>%
  ggplot(aes(
    x = p_hat, 
    y = assigned_treatment, 
    colour = assigned_dist_group
  )) +
  geom_point(size = 2) +
  theme_minimal() +
  theme(legend.position = "bottom")


ggsave("temp-plots/glm-cts-fit.png", width = 8, height= 6, dpi = 500)

stop()


library(marginaleffects)

rf_p_hat = predictions(
  rf_glm_fit, 
  newdata = datagrid(
    assigned_treatment = treat_levels, 
    assigned_dist_group = c("close", "far")
  )
)


rf_plot_df = rf_p_hat %>%
  as_tibble() %>%
    select(
      treatment = assigned_treatment, 
      assigned_dist_group, 
      value = estimate, conf.low, 
      conf.high
      )  %>%
  mutate(name = "frequentist probit") %>%
  mutate(type = "FREQ")



# rf_plot_df = belief_data %>%
#   select(obs_know_person, knows_other_dewormed, assigned_treatment, assigned_dist_group) %>%
#   glm(
#     data = .,
#     knows_other_dewormed/obs_know_person ~ 0 + assigned_treatment:assigned_dist_group,
#     weights = obs_know_person, 
#     family = binomial(link = "logit")
#   ) %>%
#   tidy(conf.int = TRUE) %>%
#   mutate(
#     treatment = str_extract(term, "(?<=treatment).*(?=:)"), 
#     assigned_dist_group = str_extract(term, "(?<=dist_group).*$")
#   ) %>%
#   rename(value = estimate)




ndgt_plot_df = ndgt_model_draws %>%
  group_by(k) %>%
  summarise(
    pr = rvar_mean(obs_prob_1ord, na.rm = TRUE)
  ) %>%
  left_join(k_dist_cw_df) %>%
  mutate(assigned_dist_group = rep(c("close", "far"), each = 4))  %>%
  rename(value = pr) %>%
  median_qi(value, na.rm = TRUE) %>%
  to_broom_names() %>%
  mutate(name = "imputed distance", type = "BAYES")

plot_df = joint_dist_model_draws %>%
  # filter(assigned_dist_group == "close") %>%
  group_by(
    treatment, assigned_dist_group
  ) %>%
  mutate(
    ed_p_hat = inv_logit(ed_pred), 
    stan_pred_p_hat = inv_logit(obs_lin_pred_1ord), 
    stan_p_hat = obs_prob_1ord 
  ) %>%
  select(contains("p_hat")) %>%
  summarise(across(contains("p_hat"), rvar_mean, na.rm = TRUE))  %>%
  pivot_longer(
      contains("p_hat")
  ) %>%
  median_qi(value, na.rm = TRUE) %>%
  to_broom_names() %>%
  mutate(type = "BAYES") %>%
  bind_rows(
    ndgt_plot_df, 
    rf_plot_df
  ) %>%
  mutate(treatment = factor(treatment, levels = treat_levels)) %>%
  mutate( 
    name = factor(name, levels = c("frequentist probit", "imputed distance", "ed_p_hat", "stan_pred_p_hat", "stan_p_hat"))
  )


plot_df %>%
  filter(name != "frequentist probit") %>%
  filter(!(name %in% c("stan_pred_p_hat", "stan_p_hat"))) %>%
  ggplot(aes(
    x = value, 
    xmin = conf.low, 
    xmax = conf.high, 
    y = treatment, 
    colour = name
  )) +
  geom_pointrange(
    position = position_dodge(0.5)
  ) +
  geom_pointrange(
    data = plot_df %>%
      filter(name == "frequentist probit"),
    position = position_dodge(1), 
    colour = "black", 
    size = 1
  ) +
  facet_wrap(
    ~assigned_dist_group, 
    ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggthemes::scale_color_canva(palette = "Primary colors with a vibrant twist") +
  labs(
    caption = "Black indicates frequentist probit.", 
    title = str_glue("1ord Beliefs, {models_we_want}")
  )


ggsave(
  str_glue("temp-plots/{models_we_want}_figuring-out-close-far.png"),
  width = 10,
  height = 10, 
  dpi = 500)


stop()

ndgt_model_draws %>%
  group_by(k) %>%
  summarise(
    pr = rvar_mean(obs_prob_1ord)
  ) %>%
  mutate(pr = mean(pr)) %>%
  left_join(k_dist_cw_df) %>%
  mutate(dist = rep(c("close", "far"), each = 4)) %>%
  select(-k) %>%
  spread(treatment, pr) 

ndgt_model_draws %>%
  left_join(
    belief_data %>%
      mutate(j = 1:n()) %>%
      select(j, assigned_dist_group))  %>%
    left_join(k_dist_cw_df) %>%
  group_by(treatment, assigned_dist_group) %>%
  summarise(
    pr = rvar_mean(obs_prob_1ord)
  ) %>%
  mutate(pr = mean(pr)) %>%
  spread(treatment, pr)




  ## here ed

comp_dist_model_draws %>%
  filter(assigned_dist_group == "far") %>%
  group_by(
    k, assigned_dist_group
  ) %>%
  summarise(
    ed_p_hat = mean(rvar_mean(inv_logit(ed_pred))),
    stan_pred_p_hat = mean(rvar_mean(inv_logit(obs_lin_pred_1ord))),
    stan_p_hat = mean(rvar_mean(obs_prob_1ord))
  ) %>%
  arrange(assigned_dist_group) %>%
  left_join(k_cw_df)

comp_dist_model_draws %>%
  # select(ed_pred, obs_lin_pred_1ord) %>%
  mutate(diff = mean(ed_pred - obs_lin_pred_1ord))  %>%
  filter(abs(diff) > 0.01) %>%
  select(standard_cluster.dist.to.pot) %>%
  unique() %>%
  arrange(standard_cluster.dist.to.pot)





  mutate(
    centered_obs_beta_1ord 
  )

stop()

load_draws = function(fit_version, model, type = "cluster") {
  fit_obj = as_cmdstan_fit(
    str_glue("data/stan_analysis_data/dist_fit{fit_version}_{model}-1.csv"))

  if (type == "cluster") {
    lin_pred_draws = gather_rvars(
      fit_obj,
      base_mu_rep, 
      centered_cluster_beta_1ord[j,k],
      centered_cluster_dist_beta_1ord[j,k]
    ) %>%
      group_by(.variable, k) %>%
      summarise(.value = rvar_mean(.value)) 
  }

  if (type == "obs") {
    lin_pred_draws = gather_rvars(
      fit_obj,
      # base_mu_rep, 
      centered_obs_beta_1ord[j,k],
      centered_obs_dist_beta_1ord[j,k]
    ) %>%
      group_by(.variable, k) %>%
      summarise(.value = rvar_mean(.value)) 
  }

  p_hat_draws = spread_rvars(
    fit_obj, 
    obs_prob_1ord[j, k], 
    obs_lin_pred_1ord[j, k]
  )
  return(lst(
    p_hat_draws,
    lin_pred_draws
  ))
}




# p_hat_71_struct = load_draws(71, "STRUCTURAL_LINEAR_U_SHOCKS")


# p_hat_80_struct_phat = load_draws(80, "structural_linear_u_shocks_phat_mu_rep")

p_hat_71_struct_obs = load_draws(
  82, 
  "STRUCTURAL_LINEAR_U_SHOCKS", 
  type = "obs")



stop()


clean_p_hat_draws = p_hat_71_struct_obs$p_hat_draws %>%
  left_join(
    belief_data %>%
      select(j, dist = cluster.dist.to.pot, dist_std = standard_cluster.dist.to.pot)
  ) %>%
  left_join(
    k_cw_df
  )

clean_p_hat_summary = clean_p_hat_draws %>%
  group_by(treatment, dist_group) %>%
  summarise(
    .value = rvar_mean(obs_prob_1ord, na.rm = TRUE)
  ) 
  # mutate(.value = mean(.value, na.rm = TRUE))



lin_pred_draws = bind_rows(
  # p_hat_71_struct$lin_pred_draws %>% mutate(fit_version = 71, obs_version = "cluster"),
  # p_hat_80_struct_phat$lin_pred_draws %>% mutate(fit_version = 80, obs_version = "cluster"),
  p_hat_71_struct_obs$lin_pred_draws %>% mutate(fit_version = 71, obs_version = "obs"),
  # p_hat_80_struct_phat_obs$lin_pred_draws %>% mutate(fit_version = 80, obs_version = "obs")
) %>%
  mutate(
    .variable = case_when(
      str_detect(.variable, "centered_.{3,7}_beta_1ord") ~ "beta_1ord",
      str_detect(.variable, "centered_.{3,7}_dist_beta_1ord") ~ "dist_beta_1ord",
      TRUE ~ .variable
    )
  )

wide_lin_pred_draws = lin_pred_draws %>%
  filter(.variable != "base_mu_rep") %>%
  group_by(fit_version, obs_version) %>%
  mutate(control_beta = .value[.variable == "beta_1ord" & k == 1]) %>% 
  mutate(
    control_dist_beta = .value[.variable == "dist_beta_1ord" & k == 1]
  ) %>%
  mutate(control = k == 1) %>%
  pivot_wider(
    names_from = .variable, values_from = .value
  )  %>%
  mutate(
    treatment = factor(treat_levels, levels = treat_levels)
  ) 


wide_lin_pred_draws





  ed = calculate_belief_latent_predictor(
    wide_lin_pred_draws[[1, "beta_1ord"]], 
    wide_lin_pred_draws[[1, "dist_beta_1ord"]], 
    unique_dists[[1]], 
    wide_lin_pred_draws[[1, "control_beta"]], 
    wide_lin_pred_draws[[1, "control_dist_beta"]], 
    control = TRUE
  )


wide_lin_pred_draws %>%
  mutate(
    lin_pred = 
  ) %>%
  select( 
    val
  )

ed

clean_p_hat_draws %>%
  filter(dist_std == unique_dists[[1]]) %>%
  filter(k == 1) %>%
  filter(j == 1)





calculate_dist_linear_predictor = function(wide_data, k, dist) {
  calculate_belief_latent_predictor(
    wide_data[[k, "beta_1ord"]], 
    wide_data[[k, "dist_beta_1ord"]], 
    dist, 
    wide_data[[k, "control_beta"]], 
    wide_data[[k, "control_dist_beta"]], 
    control = k == 1
  )
}



# calculate_p_hat = function(wide_data, dist) {
#   wide_data %>%
#     mutate(not_control = k != 1) %>%
#     mutate(dist_std = dist) %>%
#     mutate(
#       p_hat = inv_logit(beta_1ord + dist_beta_1ord*dist)
#       ) %>%
#     select(treatment, dist_std, p_hat, fit_version, everything())
# }


sd_of_dist = read_rds("temp-data/sd_of_dist.rds")


distances = seq(from = 0, to = 2500, length.out = 20) / sd_of_dist


lin_pred_dist_draws = map_dfr(
  1:4, 
  ~ {
      belief_data %>%
        select( 
          assigned_dist_group, dist_std = standard_cluster.dist.to.pot
        ) %>% 
        mutate(k = .x) %>%
        mutate(
          lin_pred = map2(
            dist_std,
            k,
            ~calculate_dist_linear_predictor(
              wide_lin_pred_draws, 
              k = .y,
              dist = .x
            )
          )
        )
  }
) %>%
  unnest(lin_pred)


wide_lin_pred_draws




lin_pred_dist_draws = lin_pred_dist_draws %>%
  group_by(k) %>%
  mutate(j = 1:n()) %>%
  left_join(
    k_cw_df
  )  %>%
  ungroup()

j_idx = 10
lin_pred_dist_draws %>%
  filter(j == j_idx) %>%
  select(assigned_dist_group, dist_std, treatment, lin_pred) %>%
  mutate(p_hat = inv_logit(lin_pred))

clean_p_hat_draws %>%
  filter(j == j_idx)  %>%
  select(dist_group, dist_std, treatment, obs_lin_pred_1ord, obs_prob_1ord)

clean_p_hat_draws %>%
  filter(j == 1) 





clean_p_hat_draws %>%
  slice(10)








lin_pred_dist_draws %>%
  mutate(
    p_hat = inv_logit(lin_pred)
  ) %>%
  group_by(k, assigned_dist_group) %>%
  summarise(
    p_hat = rvar_mean(p_hat), 
    lin_pred = rvar_mean(lin_pred)
  ) %>%
  mutate(
    p_hat = mean(p_hat), 
    lin_pred = mean(lin_pred)
  ) %>%
  left_join(
    k_cw_df
  ) %>%
  arrange(assigned_dist_group, treatment)




lin_pred_dist_draws %>%
  mutate(
    p_hat = inv_logit(lin_pred)
  ) %>%
  group_by(k, assigned_dist_group) %>%
  summarise(
    p_hat = rvar_mean(p_hat), 
    lin_pred = rvar_mean(lin_pred)
  ) %>%
  mutate(
    p_hat = mean(p_hat), 
    lin_pred = mean(lin_pred)
  ) %>%
  left_join(
    k_cw_df
  ) %>%
  write_csv(
    "temp-data/80-lin-pred-p-hat.csv"
  )


comp_lin_pred = bind_rows(
  read_csv(
    "temp-data/80-lin-pred-p-hat.csv"
  ) %>% mutate(fit = 80),
  read_csv(
    "temp-data/71-lin-pred-p-hat.csv"
  ) %>% mutate(fit = 71)
)

comp_lin_pred %>%
  select(-k, -dist_group) %>%
  select(-lin_pred)  %>%
  spread(
    fit, p_hat
  )


stan_data

stop()


clean_p_hat_summary %>%
  summarise(.value = mean(.value, na.rm = TRUE))


pred_p_hat_draws = belief_data %>%
  select( 
    assigned_dist_group, standard_cluster.dist.to.pot
  ) %>%
  mutate(
    p_hat_draws = map(
      standard_cluster.dist.to.pot, 
      ~calculate_p_hat(
        wide_lin_pred_draws,
         .x
      )
    )
  )

clean_pred_p_hat_draws = pred_p_hat_draws %>%
  unnest_wider(p_hat_draws) %>%
  unnest(cols = c(treatment, dist_std, p_hat))


clean_pred_p_hat_summary = clean_pred_p_hat_draws %>%
  group_by(treatment, assigned_dist_group) %>%
  summarise(p_hat = rvar_mean(p_hat, na.rm = TRUE)) 


inner_join(
  clean_p_hat_summary,
  clean_pred_p_hat_summary, 
  by = c("treatment", "dist_group" = "assigned_dist_group")
) %>%
  mutate(
    diff = .value - p_hat, 
    .value = mean(.value, na.rm = TRUE), 
    p_hat = mean(p_hat, na.rm = TRUE)
  ) %>%
  select(treatment, dist_group, .value, p_hat, diff) %>%
  median_qi(diff, na.rm = TRUE)  %>%
  to_broom_names() %>%
  ggplot(aes(
    x = diff, 
    xmin = conf.low, 
    xmax = conf.high,
    y = treatment, 
    colour = dist_group
  )) +
  geom_pointrange(
    position = position_dodge(0.5)
  ) +
  geom_vline(
    xintercept = 0, 
    linetype = "longdash"
  ) +
  theme_minimal()

stop()

clean_p_hat_draws


pred_p_hat_draws = map_dfr(
  belief_data$standard_cluster.dist.to.pot,
  ~calculate_p_hat(
    wide_lin_pred_draws, 
    .x
  )
)


clean_p_hat_draws %>%
  filter(dist_std == unique_dists[[1]]) %>%
  group_by(treatment) %>%
  summarise(
    .value = rvar_mean(.value, na.rm = TRUE)
  ) %>%
  mutate(.value = mean(.value, na.rm = TRUE)) %>%
  mutate(treatment = factor(treatment, levels = treat_levels)) %>%
  arrange(treatment)

pred_p_hat_draws %>%
  filter(
    dist_std == unique_dists[[1]]
  )  %>%
  group_by(treatment)  %>%
  summarise(
    p_hat = rvar_mean(p_hat)
  ) %>%
  mutate(p_hat = mean(p_hat)) %>%
  arrange(treatment)

pred_p_hat_draws %>%
  mutate(
    p_hat = mean(p_hat)
  ) %>%
  ggplot(aes(
    x = dist, 
    y = p_hat, 
    colour = treatment
  )) +
  geom_point()




p_hat_df %>%
  filter(fit_version == 71) %>%
  mutate(
    p_hat = mean(p_hat), 
    distance = dist*sd_of_dist
  ) %>%
  ggplot(aes(
    x = distance, 
    y = p_hat, 
    colour = treatment, 
    linetype = factor(obs_version)
  )) +
  geom_point() +
  geom_line() +
  theme_bw() +
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
  ) +
  facet_grid(~fit_version) +
  geom_hline( 
    yintercept = 1, 
    linetype = "longdash"
  )

p_hat_df %>%
  mutate(
    p_hat = mean(p_hat), 
    distance = dist*sd_of_dist
  ) %>%
  ggplot(aes(
    x = distance, 
    y = p_hat, 
    colour = treatment, 
    linetype = factor(fit_version)
  )) +
  geom_point() +
  geom_line() +
  theme_bw() +
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
  ) +
  facet_grid(obs_version~fit_version) +
  geom_hline( 
    yintercept = 1, 
    linetype = "longdash"
  )


p_hat_df %>%
  mutate(
    p_hat = mean(p_hat), 
    distance = dist*sd_of_dist
  ) %>%
  ggplot(aes(
    x = distance, 
    y = p_hat, 
    colour = obs_version, 
    linetype = factor(fit_version)
  )) +
  geom_point() +
  geom_line() +
  theme_bw() +
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
  ) +
  facet_wrap(~treatment, ncol = 1) +
  geom_hline( 
    yintercept = 1, 
    linetype = "longdash"
  )

ggsave("temp-plots/anne-check/pr-1ord-continuous-distance.png", 
width = 8, 
height = 6, 
dpi = 500)



belief_data = analysis_data %>%
  filter(obs_know_person > 0) 

p_hat_belief_data = belief_data %>%
  select(
    assigned_dist_group, 
    assigned_treatment,
    standard_cluster.dist.to.pot) %>%
  mutate(
    p_hat = map(
      standard_cluster.dist.to.pot, 
      ~calculate_p_hat(
        wide_p_hat_draws, 
        .x
      ) %>% as_tibble() %>% select(p_hat, treatment, dist, fit_version, obs_version)
    )
  )



long_p_hat_belief_data = p_hat_belief_data %>%
  unnest(p_hat)


long_p_hat_belief_data %>%
  filter(fit_version == 71) %>%
  filter(treatment == assigned_treatment) %>%
  filter(assigned_dist_group == "close") %>%
  group_by(
    fit_version,
    obs_version,
    assigned_dist_group, 
    assigned_treatment
  ) %>%
  summarise(
    p_hat = rvar_mean(p_hat)
  ) %>%
  mutate(
    p_hat = mean(p_hat)
  )


long_p_hat_belief_data %>%
  mutate(dist = dist*sd_of_dist) %>%
  mutate(p_hat = mean(p_hat)) %>%
  filter(assigned_treatment == treatment) %>%
  ggplot(aes(
    x = dist, 
    y = p_hat, 
    colour = treatment
  )) +
  geom_point() +
  facet_wrap(~fit_version, ncol = 1)  +
  geom_point(
    data = . %>%
      group_by(
        fit_version, 
        treatment, 
        assigned_dist_group
      ) %>%
      summarise(p_hat = mean(p_hat)) %>%
      mutate(dist = if_else(assigned_dist_group == "close", 1250/2, 1250 + (2500/2))), 
    colour = "black", 
    size = 2
  ) +
  geom_point(
    data = belief_data %>%
      select(obs_know_person, knows_other_dewormed, assigned_treatment, assigned_dist_group) %>%
      group_by(
        assigned_treatment, 
        assigned_dist_group
      ) %>%
      summarise(
        p_hat = mean(knows_other_dewormed/obs_know_person)
      ) %>%
      mutate(dist = if_else(assigned_dist_group == "close", 1250/2 + 200, 200 + 1250 + (2500/2))), 
      colour = "red", 
      size = 4
  )



belief_data %>%
  select(obs_know_person, knows_other_dewormed, assigned_treatment, assigned_dist_group) %>%
  group_by(
    assigned_treatment, 
    assigned_dist_group
  ) %>%
  summarise(
    pr = mean(knows_other_dewormed/obs_know_person),
    sd =  1/sqrt(n())*pr*(1-pr), 
    conf.low = pr - 1.96*sd,
    conf.high = pr + 1.96*sd
  ) %>%
  ggplot(aes(
    x = pr, 
    xmin = conf.low, 
    xmax = conf.high,
    y = assigned_treatment, 
    colour = assigned_dist_group
  )) +
  geom_pointrange(
    position = position_dodge(0.5)
  ) +
  # facet_wrap(~assigned_dist_group, ncol = 1)  +
  theme_bw() + 
  guides(colour = "none") +
  geom_vline(xintercept = 0.8) +
  geom_vline(xintercept = 0.85, linetype = "longdash") 

belief_data %>%
  select(obs_know_person, knows_other_dewormed, assigned_treatment, assigned_dist_group) %>%
  glm(
    data = .,
    knows_other_dewormed/obs_know_person ~ 0 + assigned_treatment:assigned_dist_group,
    weights = obs_know_person, 
    family = binomial(link = "logit")
  ) %>%
  tidy()


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
load(file.path("temp-data", str_interp("processed_dist_fit80_lite.RData")))

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
