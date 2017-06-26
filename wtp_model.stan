data {
  int<lower = 0> num_wtp_obs;
  int<lower = 0> num_strata;
  
  int<lower = 1, upper = num_wtp_obs> wtp_strata_sizes[num_strata]; 

  vector<lower = -1, upper = 1>[num_wtp_obs] gift_choice;
  vector<lower = -1, upper = 1>[num_wtp_obs] wtp_response;
  vector<lower = 0>[num_wtp_obs] wtp_offer;
  
  real<lower = 0> tau_mu_wtp_diff;
  real<lower = 0> mu_wtp_df_student_t;
  
  real<lower = 0> tau_sigma_wtp_diff;
  real<lower = 0> sigma_wtp_df_student_t;
  
  real<lower = 0, upper = 10> wtp_utility_df; // TODO put hyperprior on this parameter
}

parameters {
  real hyper_mu_wtp_diff_raw;
  vector[num_strata] mu_wtp_diff_raw;
  real<lower = 0> sigma_wtp_diff;
}

transformed parameters {
  real hyper_mu_wtp_diff = hyper_mu_wtp_diff_raw * tau_mu_wtp_diff;
  vector[num_strata] mu_wtp_diff = hyper_mu_wtp_diff + mu_wtp_diff_raw * tau_mu_wtp_diff;
}

model {
  hyper_mu_wtp_diff_raw ~ student_t(mu_wtp_df_student_t, 0, 1);
  mu_wtp_diff_raw ~ student_t(mu_wtp_df_student_t, 0, 1);
  sigma_wtp_diff ~ student_t(sigma_wtp_df_student_t, 0, tau_sigma_wtp_diff); 
 
  {
    int stratum_pos = 1;
    
    for (stratum_index in 1:num_strata) {
      int curr_stratum_size = wtp_strata_sizes[stratum_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      for (i in stratum_pos:stratum_end) {
        if (wtp_response[i] == -1) {
          if (gift_choice[i] == -1) {
            target += student_t_lcdf(- wtp_offer[i] | wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff);
          } else {
            target += student_t_lccdf(wtp_offer[i] | wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff);
          }
        } else {
          target += log(gift_choice[i] * (student_t_cdf(gift_choice[i] * wtp_offer[i], wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff) 
                    - student_t_cdf(0, wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff)));
        } 
      }
      
      stratum_pos = stratum_end + 1;
    }
  }
}

generated quantities {
  real hyper_prob_prefer_calendar = 1 - student_t_cdf(0, wtp_utility_df, hyper_mu_wtp_diff, sigma_wtp_diff);
  real prob_prefer_calendar[num_strata];
  
  for (stratum_index in 1:num_strata) {
    prob_prefer_calendar[stratum_index] = 1 - student_t_cdf(0, wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff);
  } 
}
