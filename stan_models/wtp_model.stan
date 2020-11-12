data {
  int<lower = 0> num_wtp_obs;
  int<lower = 0> num_strata;
  // int<lower = 1> num_clusters; // In entire study
  // int<lower = 1> num_wtp_clusters; // In WTP study only
  
  // int<lower = 1, upper = num_clusters> wtp_cluster_id[num_wtp_obs];
  int<lower = 1, upper = num_wtp_obs> wtp_strata_sizes[num_strata]; 
  // int<lower = 1, upper = num_wtp_obs> wtp_cluster_sizes[num_wtp_clusters];

  // Calendar = 1
  // Bracelet = -1
  vector<lower = -1, upper = 1>[num_wtp_obs] gift_choice;
  
  // Switch = 1
  // Keep = -1
  vector<lower = -1, upper = 1>[num_wtp_obs] wtp_response;
  
  vector<lower = 0>[num_wtp_obs] wtp_offer;
}

transformed data {
  row_vector<lower = 0, upper = 1>[num_strata] strata_prop = to_row_vector(wtp_strata_sizes) / sum(wtp_strata_sizes);
  real strata_mu_student_df = 3;
  
  vector<lower = 0>[num_wtp_obs] scaled_wtp_offer = wtp_offer / max(wtp_offer);
}

parameters {
  real hyper_mu;
  
  vector[num_strata] raw_strata_mu;
  real<lower = 0> strata_mu_tau;
  
  real<lower = 0> sigma;
}

transformed parameters {
  vector[num_strata] strata_mu = hyper_mu + raw_strata_mu * strata_mu_tau;
}

model {
  hyper_mu ~ normal(0, 2);
 
  raw_strata_mu ~ std_normal();
  strata_mu_tau ~ normal(0, 1);
  
  sigma ~ normal(0, 1); 
 
  { 
    int wtp_stratum_pos = 1;

    for (stratum_index in 1:num_strata) {
      int curr_wtp_stratum_size = wtp_strata_sizes[stratum_index];
      int wtp_stratum_end = wtp_stratum_pos + curr_wtp_stratum_size - 1;

      for (wtp_obs_index in wtp_stratum_pos:wtp_stratum_end) {
        if (wtp_response[wtp_obs_index] == -1) { // Decided to KEEP
          if (gift_choice[wtp_obs_index] == -1) { // Initial choice: BRACELET
            target += log(Phi_approx((- scaled_wtp_offer[wtp_obs_index] - strata_mu[stratum_index]) / sigma));
          } else { // Initial choice: CALENDAR
            target += log(1 - Phi_approx((scaled_wtp_offer[wtp_obs_index] - strata_mu[stratum_index]) / sigma));
          }
        } else { // Decided to SWITCH
          target += log(gift_choice[wtp_obs_index] *
            (Phi_approx((gift_choice[wtp_obs_index] * scaled_wtp_offer[wtp_obs_index] - strata_mu[stratum_index]) / sigma)
             - Phi_approx(- strata_mu[stratum_index] / sigma)));
        }
      }

      wtp_stratum_pos = wtp_stratum_end + 1;
    }
  }
}

generated quantities {
  real<lower = 0, upper = 1> prob_prefer_calendar;
  vector<lower = 0, upper = 1>[num_strata] strata_prob_prefer_calendar;
  
  // vector<lower = -1, upper = 1>[num_wtp_obs] sim_gift_choice;
  // vector<lower = -1, upper = 1>[num_wtp_obs] sim_wtp_response; 
  
  {
    int stratum_pos = 1;
    
    for (stratum_index in 1:num_strata) {
      int curr_stratum_size = wtp_strata_sizes[stratum_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      strata_prob_prefer_calendar[stratum_index] = 1 - normal_cdf(0, strata_mu[stratum_index], sigma);
     
      // for (i in 1:curr_stratum_size) { 
      //   real val_diff = student_t_rng(strata_student_df[stratum_index], strata_mu[stratum_index], strata_sigma[stratum_index]);
      //   
      //   sim_gift_choice[stratum_pos + i - 1] = val_diff > 0 ? 1 : -1;
      //   sim_wtp_response[stratum_pos + i - 1] = wtp_offer[stratum_pos + i - 1] > fabs(val_diff) ? 1 : -1;
      // }
      
      stratum_pos = stratum_end + 1;
    }
  }
  
  prob_prefer_calendar = strata_prop * strata_prob_prefer_calendar;
}
