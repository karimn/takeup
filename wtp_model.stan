data {
  int<lower = 0> num_obs;

  vector<lower = -1, upper = 1>[num_obs] gift_choice;
  vector<lower = -1, upper = 1>[num_obs] response;
  vector<lower = 0>[num_obs] offer;
  
  real<lower = 0> tau_mu_diff;
  real<lower = 0> mu_df_student_t;
  
  real<lower = 0> tau_sigma_diff;
  real<lower = 0> sigma_df_student_t;
}

parameters {
  real mu_diff_raw;
  real<lower = 0> sigma_diff;
}

transformed parameters {
  real mu_diff = mu_diff_raw * tau_mu_diff;
}

model {
  mu_diff_raw ~ student_t(mu_df_student_t, 0, 1);
  sigma_diff ~ student_t(sigma_df_student_t, 0, tau_sigma_diff); 
  
  for (i in 1:num_obs) {
    if (response[i] == -1) {
      if (gift_choice[i] == -1) {
        target += normal_lcdf(- offer[i] | mu_diff, sigma_diff);
      } else {
        target += normal_lccdf(offer[i] | mu_diff, sigma_diff);
      }
    } else {
      target += log(gift_choice[i] * (normal_cdf(gift_choice[i] * offer[i], mu_diff, sigma_diff) - normal_cdf(0, mu_diff, sigma_diff)));
    } 
  }
}

generated quantities {
  real prob_prefer_calendar = 1 - normal_cdf(0, mu_diff, sigma_diff);
}
