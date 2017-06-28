library(tidyverse)
library(rstan)

options(mc.cores = max(1, parallel::detectCores()))

# Data prep ---------------------------------------------------------------
 
test_data <- tibble(treated = rep(0:2, 2000),
                    re_intercept = rnorm(length(treated), 0.2, 1), # 0.2 intercept and epsilon ~ N(0, 1)
                    re_te1 = rnorm(length(treated), 0.5, 0.5),
                    re_te2 = rnorm(length(treated), 0.25, 0.1),
                    u = re_intercept + 
                      (treated == 1) * re_te1 + 
                      (treated == 2) * re_te2 + 
                      (treated == 2) * 0.25, 
                    y = 1*(u > 0)) %>% # observed outcome
  arrange(treated) %>% 
  mutate(obs_index = seq(n()))

# Analyze -----------------------------------------------------------------

test_model <- stan_model(file = "test.stan", model_name = "test_model")
test_model_fit <- sampling(test_model, iter = 100, 
                           data = 
                             lst(num_obs = nrow(test_data), 
                                 treated1 = test_data$treated == 1, 
                                 treated2 = test_data$treated == 2, 
                                 num_treated2 = sum(treated2),
                                 treated2_ids = test_data %>% filter(treated == 2) %>% pull(obs_index),
                                 mu_diff = 0.25 * 10,
                                 sigma_diff = sqrt(0.5^2 + 0.1^2) * 10,
                                 y = test_data$y),
                           sample_file = file.path("stanfit", "test_model_0.csv"),
                           control = lst(adapt_delta = 0.99, max_treedepth = 15))

print(test_model_fit, pars = c("latent_diff_raw", "latent_diff"), include = FALSE)
