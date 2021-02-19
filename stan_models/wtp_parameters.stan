real hyper_wtp_mu;

vector[use_strata_levels ? num_strata : 0] raw_strata_wtp_mu;
real<lower = 0> strata_wtp_mu_tau;

real<lower = 0> wtp_sigma;