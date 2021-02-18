vector[num_strata] strata_effect_wtp_mu = raw_strata_wtp_mu * strata_wtp_mu_tau;
vector[num_strata] strata_wtp_mu = hyper_wtp_mu + strata_effect_wtp_mu;