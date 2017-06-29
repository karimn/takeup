functions {
}

data {
  int<lower = 0> num_obs;
  int<lower = 0> num_treated2;
  
  int treated2_ids[num_treated2];
  
  int<lower = 0, upper = 1> y[num_obs];
  vector<lower = 0, upper = 1>[num_obs] treated1;
  vector<lower = 0, upper = 1>[num_obs] treated2;
  
  real mu_diff;
  real sigma_diff;
}

transformed data {
}

parameters {
  real alpha_raw;
  real te_raw;
  real<lower = 0> gamma_raw;
  real delta_raw;
  
  vector[num_treated2] latent_diff_raw;
  
  // real<lower = 0> tau_alpha;
  // real<lower = 0> tau_te;
}

transformed parameters {
  real alpha = alpha_raw * 5;
  real te = te_raw * 5;
  real<lower = 0> gamma = gamma_raw * 1;
  real delta = delta_raw * 5;
  
  vector[num_treated2] latent_diff = mu_diff + latent_diff_raw * sigma_diff;
}

model {
  alpha_raw ~ student_t(7, 0, 1);
  te_raw ~ student_t(7, 0, 1);
  gamma_raw ~ student_t(7, 0, 1);
  delta_raw ~ student_t(7, 0, 1);
  
  latent_diff_raw ~ normal(0, 1);
  
  // tau_alpha ~ student_t(7, 0, 10);
  // tau_te ~ student_t(7, 0, 10);
  
  { 
    vector[num_obs] full_latent_diff = rep_vector(0, num_obs);
    full_latent_diff[treated2_ids] = latent_diff;
    
    y ~ bernoulli(Phi(alpha + gamma * te * (treated1 + treated2) - gamma * (full_latent_diff .* treated2) + delta * treated2));
  }
}

generated quantities {
  real t0 = Phi(alpha);
  real t1 = Phi(alpha + gamma * te); 
  real t2; 
  real ape = t1 - t0;
  
  {
    vector[num_obs] simul_full_latent_diff = rep_vector(0, num_obs);
  
    for (i in 1:num_treated2) { 
      simul_full_latent_diff[treated2_ids[i]] = normal_rng(mu_diff, sigma_diff);
    }
    
    t2 = mean(Phi(alpha + gamma * te - gamma * simul_full_latent_diff + delta));
  }
}
