library(tidyverse)
library(rstan)

options(mc.cores = max(1, parallel::detectCores()))

test_data <- tibble(treated = rep(0:1, 2000),
                    re_intercept = rnorm(length(treated), 0.2, 1),
                    re_te = rnorm(length(treated), 0.5, 0.5),
                    u = re_intercept + treated * re_te,
                    y = 1*(u > 0))

test_model <- stan_model(file = "test.stan", model_name = "test_model")
test_model_fit <- sampling(test_model, iter = 8000, data = lst(num_obs = nrow(test_data), treated = test_data$treated, y = test_data$y))

print(test_model_fit, pars = c("te_noise", "te_noise_raw", "epsilon", "y_ast"), include = FALSE)
