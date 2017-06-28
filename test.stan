functions {
}

data {
  int<lower = 0> num_obs;
  
  int<lower = 0, upper = 1> y[num_obs];
  vector<lower = 0, upper = 1>[num_obs] treated;
}

transformed data {
}

parameters {
  real alpha_raw;
  real te_raw;
  // vector[num_obs] epsilon;
  vector[num_obs] te_noise_raw;
  
  real<lower = 0> sigma_te;
  
  // real<lower = 0> tau_alpha;
  // real<lower = 0> tau_te;
}

transformed parameters {
  real alpha = alpha_raw * 5;
  real te = te_raw * 5;
  vector[num_obs] te_noise = te_noise_raw * sigma_te;
  
  // vector[num_obs] y_ast = alpha + (te + te_noise) .* treated + epsilon;
}

model {
  alpha_raw ~ student_t(7, 0, 2.5);
  te_raw ~ student_t(7, 0, 2.5);
  
  // tau_alpha ~ student_t(7, 0, 10);
  // tau_te ~ student_t(7, 0, 10);
  
  sigma_te ~ student_t(7, 0, 2.5);
    
  te_noise_raw ~ normal(0, 1);
  // epsilon ~ normal(0, 1);
    
  // y ~ bernoulli(Phi((alpha + te * treated)./sqrt(1 + square(sigma_te) * treated)));
  // y ~ bernoulli(Phi(alpha + (te + te_noise) .* treated));
  
  for (i in 1:num_obs) {
    if (y[i] == 1) {
      target += normal_lcdf(alpha + (te + te_noise[i]) * treated[i] | 0, 1);
    } else {
      target += normal_lccdf(alpha + (te + te_noise[i]) * treated[i] | 0, 1);
    }
  }
    
  // y ~ bernoulli_logit(alpha + (te + te_noise * sigma_te) .* treated);
  // y ~ bernoulli_logit((alpha + te * treated)./(sqrt(1 + square(sigma_te) * treated)));
}

generated quantities {
}
