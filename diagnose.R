library(bayesplot)
color_scheme_set("darkgray")

post <- dist_fit$STRUCTURAL_LINEAR_U_SHOCKS$draws()
np <- nuts_params(dist_fit$STRUCTURAL_LINEAR_U_SHOCKS)

mcmc_parcoord(post, np = np)

# c("raw_u_sd[1]", "raw_u_sd[2]", "raw_u_sd[3]",  "wtp_sigma", "wtp_value_utility") %>% 
c("raw_u_sd[1]", "raw_u_sd[2]", "raw_u_sd[3]", "structural_cluster_obs_v[100]",  "wtp_sigma", "wtp_value_utility") %>% 
  mcmc_pairs(post[,, .], np = np, pars = ., off_diag_args = list(size = 0.75))
