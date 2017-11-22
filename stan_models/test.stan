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
  real<lower = 0> te_raw;
  real<lower = 0> gamma_raw;
  real<lower = 0> delta_raw;
  
  vector[num_treated2] latent_diff_raw;
  
  // real<lower = 0> tau_alpha;
  // real<lower = 0> tau_te;
}

transformed parameters {
  real alpha = alpha_raw * 5;
  real te = te_raw * 5;
  real gamma = gamma_raw * 1;
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
    
    y ~ bernoulli(Phi(alpha + te * (treated1 + treated2) + gamma * full_latent_diff + delta * treated2));
  }
}

generated quantities {
  real t0 = Phi(alpha);
  // real util_te = gamma * te;
  real t1 = Phi(alpha + te); 
  real t2; // = Phi(alpha + te + delta);
  // real t2_mean_latent = Phi(alpha + util_te - gamma * mu_diff + delta); 
  // real ape1 = t1 - t0;
  
  {
    vector[num_treated2] simul_latent_diff = rep_vector(0, num_treated2);

    for (i in 1:num_treated2) {
      simul_latent_diff[i] = normal_rng(mu_diff, sigma_diff);
    }

    // t2 = mean(Phi(alpha + util_te - gamma * simul_latent_diff));
    // t2 = mean(Phi(alpha + util_te - gamma * simul_latent_diff + delta));
    t2 = mean(Phi(alpha + te + gamma * simul_latent_diff + delta));
  }
}
