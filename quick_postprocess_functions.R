library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)


## Load analysis data
load(file.path("data", "analysis.RData"))
standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
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

sd_of_dist = sd(analysis_data$cluster.dist.to.pot)

## Load Stan output
load_param_draws = function(fit_version, 
                            model, 
                            chain, 
                            prior_predictive = FALSE, 
                            input_path,
                            ...) {
  if (prior_predictive == TRUE) {
    fit_str = file.path(
        input_path,
        str_glue(
        "dist_prior{fit_version}_{model}-{chain}.csv"
        )
    )
  } else {
    fit_str = file.path(
        input_path,
        str_glue(
        "dist_fit{fit_version}_{model}-{chain}.csv"
        )
    )
  }

  fit_obj = as_cmdstan_fit(fit_str)
  draws = spread_rvars(
    fit_obj,
    ...
  ) %>%
    mutate(model = model, fit_version = fit_version, fit_type = if_else(prior_predictive, "prior-predict", "fit"))
  return(draws)
}