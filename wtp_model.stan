data {
  int<lower = 0> num_obs;
  int<lower = 0> num_strata;
  
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 

  vector<lower = -1, upper = 1>[num_obs] gift_choice;
  vector<lower = -1, upper = 1>[num_obs] response;
  vector<lower = 0>[num_obs] offer;
  
  real<lower = 0> tau_mu_diff;
  real<lower = 0> mu_df_student_t;
  
  real<lower = 0> tau_sigma_diff;
  real<lower = 0> sigma_df_student_t;
}

transformed data {
  real<lower = 0, upper = 10> utility_df = 3; // TODO put hyperprior on this parameter
}

parameters {
  real hyper_mu_diff_raw;
  vector[num_strata] mu_diff_raw;
  real<lower = 0> sigma_diff;
}

transformed parameters {
  real hyper_mu_diff = hyper_mu_diff_raw * tau_mu_diff;
  vector[num_strata] mu_diff = hyper_mu_diff + mu_diff_raw * tau_mu_diff;
}

model {
  hyper_mu_diff_raw ~ student_t(mu_df_student_t, 0, 1);
  mu_diff_raw ~ student_t(mu_df_student_t, 0, 1);
  sigma_diff ~ student_t(sigma_df_student_t, 0, tau_sigma_diff); 
 
  {
    int stratum_pos = 1;
    
    for (stratum_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[stratum_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      for (i in stratum_pos:stratum_end) {
        if (response[i] == -1) {
          if (gift_choice[i] == -1) {
            target += student_t_lcdf(- offer[i] | utility_df, mu_diff[stratum_index], sigma_diff);
          } else {
            target += student_t_lccdf(offer[i] | utility_df, mu_diff[stratum_index], sigma_diff);
          }
        } else {
          target += log(gift_choice[i] * (student_t_cdf(gift_choice[i] * offer[i], utility_df, mu_diff[stratum_index], sigma_diff) 
                    - student_t_cdf(0, utility_df, mu_diff[stratum_index], sigma_diff)));
        } 
      }
      
      stratum_pos = stratum_end + 1;
    }
  }
}

generated quantities {
  // real prob_prefer_calendar = 1 - normal_cdf(0, mu_diff, sigma_diff);
  real hyper_prob_prefer_calendar = 1 - student_t_cdf(0, utility_df, hyper_mu_diff, sigma_diff);
  real prob_prefer_calendar[num_strata];
  
  for (stratum_index in 1:num_strata) {
    prob_prefer_calendar[stratum_index] = 1 - student_t_cdf(0, utility_df, mu_diff[stratum_index], sigma_diff);
  } 
}
