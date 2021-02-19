vector[num_strata] strata_effect_wtp_mu = use_strata_levels ? raw_strata_wtp_mu * strata_wtp_mu_tau : rep_vector(0, num_strata);
vector[num_strata] strata_wtp_mu = hyper_wtp_mu + strata_effect_wtp_mu;