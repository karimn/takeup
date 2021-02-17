vector<lower = 0, upper = 1>[num_strata] strata_prop = to_vector(wtp_strata_sizes) / sum(wtp_strata_sizes);
real strata_mu_student_df = 3;

vector<lower = 0>[num_wtp_obs] scaled_wtp_offer = wtp_offer / max(wtp_offer);
vector[num_preference_value_diff] scaled_preference_value_diff = preference_value_diff / max(wtp_offer);