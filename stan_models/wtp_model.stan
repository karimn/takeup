data {
  int<lower = 0> num_wtp_obs;
  int<lower = 0> num_strata;
  int<lower = 1> num_clusters; // In entire study
  int<lower = 1> num_wtp_clusters; // In WTP study only
  
  int<lower = 1, upper = num_clusters> wtp_cluster_id[num_wtp_obs];
  int<lower = 1, upper = num_wtp_obs> wtp_strata_sizes[num_strata]; 
  int<lower = 1, upper = num_wtp_obs> wtp_cluster_sizes[num_wtp_clusters];

  vector<lower = -1, upper = 1>[num_wtp_obs] gift_choice;
  vector<lower = -1, upper = 1>[num_wtp_obs] wtp_response;
  vector<lower = 0>[num_wtp_obs] wtp_offer;
  
  // real<lower = 0> tau_mu_wtp_diff;
  // real<lower = 0> mu_wtp_df_student_t;
  // 
  // real<lower = 0> tau_sigma_wtp_diff;
  // real<lower = 0> sigma_wtp_df_student_t;
  // 
  // real<lower = 0, upper = 10> wtp_utility_df; 
}

parameters {
  // real hyper_mu_wtp_diff;
  real hyper_mu;
  real<lower = 2, upper = 8> strata_mu_student_df;
  
  // vector[num_strata] mu_wtp_diff;
  real strata_mu[num_strata];
  real<lower = 0> strata_mu_tau;
  
  // real<lower = 0> sigma_wtp_diff;
  real<lower = 2, upper = 8> strata_student_df[num_strata];
  real<lower = 0> strata_sigma[num_strata];
}

model {
  // hyper_mu_wtp_diff ~ student_t(mu_wtp_df_student_t, 0, tau_mu_wtp_diff);
  hyper_mu ~ normal(0, 1000);
  // mu_wtp_diff ~ student_t(mu_wtp_df_student_t, 0, tau_mu_wtp_diff);
 
  strata_mu_tau ~ normal(0, 100); 
  strata_mu ~ student_t(strata_mu_student_df, 0, strata_mu_tau);
  
  // sigma_wtp_diff ~ student_t(sigma_wtp_df_student_t, 0, tau_sigma_wtp_diff); 
  strata_sigma ~ normal(0, 1000); 
 
  {
    int wtp_stratum_pos = 1;

    for (stratum_index in 1:num_strata) {
      int curr_wtp_stratum_size = wtp_strata_sizes[stratum_index];
      int wtp_stratum_end = wtp_stratum_pos + curr_wtp_stratum_size - 1;

      for (wtp_obs_index in wtp_stratum_pos:wtp_stratum_end) {
        if (wtp_response[wtp_obs_index] == -1) {
          if (gift_choice[wtp_obs_index] == -1) {
            // target += student_t_lcdf(- wtp_offer[wtp_obs_index] | wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff);
            target += student_t_lcdf(- wtp_offer[wtp_obs_index] | strata_student_df[stratum_index], strata_mu[stratum_index], strata_sigma[stratum_index]);
          } else {
            // target += student_t_lccdf(wtp_offer[wtp_obs_index] | wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff);
            target += student_t_lcdf(wtp_offer[wtp_obs_index] | strata_student_df[stratum_index], strata_mu[stratum_index], strata_sigma[stratum_index]);
          }
        } else {
          // target += log(gift_choice[wtp_obs_index] * (student_t_cdf(gift_choice[wtp_obs_index] * wtp_offer[wtp_obs_index], wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff) -
          //                                             student_t_cdf(0, wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff)));
          target += log(gift_choice[wtp_obs_index] * 
            (student_t_cdf(gift_choice[wtp_obs_index] * wtp_offer[wtp_obs_index], strata_student_df[stratum_index], strata_mu[stratum_index], strata_sigma[stratum_index]) - student_t_cdf(0, strata_student_df[stratum_index], strata_mu[stratum_index], strata_sigma[stratum_index])));
        }
      }

      wtp_stratum_pos = wtp_stratum_end + 1;
    }
  }
}

generated quantities {
  // real hyper_prob_prefer_calendar = 1 - student_t_cdf(0, wtp_utility_df, hyper_mu_wtp_diff, sigma_wtp_diff);
  // real hyper_prob_prefer_calendar = 1 - student_t_cdf(0, wtp_utility_df, hyper_mu_wtp_diff, sigma_wtp_diff);
  real prob_prefer_calendar[num_strata];
  
  vector<lower = -1, upper = 1>[num_wtp_obs] sim_gift_choice;
  vector<lower = -1, upper = 1>[num_wtp_obs] sim_wtp_response; 
  
  {
    int stratum_pos = 1;
    
    for (stratum_index in 1:num_strata) {
      int curr_stratum_size = wtp_strata_sizes[stratum_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      // prob_prefer_calendar[stratum_index] = 1 - student_t_cdf(0, wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff);
      prob_prefer_calendar[stratum_index] = 1 - student_t_cdf(0, strata_student_df[stratum_index], strata_mu[stratum_index], strata_sigma[stratum_index]);
     
      for (i in 1:curr_stratum_size) { 
        // real val_diff = student_t_rng(wtp_utility_df, mu_wtp_diff[stratum_index], sigma_wtp_diff);
        real val_diff = student_t_rng(strata_student_df[stratum_index], strata_mu[stratum_index], strata_sigma[stratum_index]);
        
        sim_gift_choice[stratum_pos + i - 1] = val_diff > 0 ? 1 : -1;
        sim_wtp_response[stratum_pos + i - 1] = wtp_offer[stratum_pos + i - 1] > fabs(val_diff) ? 1 : -1;
      }
      
      stratum_pos = stratum_end + 1;
    }
  }
}
