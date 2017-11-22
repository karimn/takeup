data {
  int<lower = 0> num_obs;
  int<lower = 0> num_treated2;
  
  int treated2_ids[num_treated2];
  
  int<lower = 0, upper = 1> y[num_obs];
  vector<lower = 0, upper = 1>[num_obs] treated1;
  vector<lower = 0, upper = 1>[num_obs] treated2;
 
  # Assume already know from other data set 
  real mu_diff;
  real sigma_diff;
}

transformed data {
  real<lower = 0, upper = 1> prob_prefer_te = 1 - normal_cdf(0, mu_diff, sigma_diff); 
}

parameters {
  real alpha_raw;
  real te_raw;
  // real extra_te_raw;
  real delta_raw;
  
  // vector[num_treated2] latent_diff_raw;
  
  // real<lower = 0> tau_alpha;
  // real<lower = 0> tau_te;
}

transformed parameters {
  real alpha = alpha_raw * 5;
  real te = te_raw * 5;
  // real extra_te = extra_te_raw * 5;
  real delta = delta_raw * 5;
  
  // vector[num_treated2] latent_diff = mu_diff + latent_diff_raw * sigma_diff;
}

model {
  alpha_raw ~ student_t(7, 0, 1);
  te_raw ~ student_t(7, 0, 1);
  // extra_te_raw ~ student_t(7, 0, 1);
  delta_raw ~ student_t(7, 0, 1);
  
  // tau_alpha ~ student_t(7, 0, 10);
  // tau_te ~ student_t(7, 0, 10);
  
  for (i in 1:num_obs) {
    target += log_mix(1 - prob_prefer_te, bernoulli_lpmf(y[i] | Phi(alpha + te * (treated1[i] + treated2[i]))), # + extra_te * treated2[i])),
                                          bernoulli_lpmf(y[i] | Phi(alpha + te * (treated1[i] + treated2[i]) + delta * treated2[i])));
  }
}

generated quantities {
  real t0 = Phi(alpha);
  real t1 = Phi(alpha + te); 
  real t2_pref_te = Phi(alpha + te);
  real t2_pref_delta = Phi(alpha + te + delta); 
  real t2 = prob_prefer_te * Phi(alpha + te) + (1 - prob_prefer_te) * Phi(alpha + te + delta);
}
