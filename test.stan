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
  real<lower = 0> te_raw;
  // real<lower = 0> sigma_te;
  
  // real<lower = 0> tau_alpha;
  // real<lower = 0> tau_te;
}

transformed parameters {
  real alpha = alpha_raw * 5;
  real te = te_raw * 5;
}

model {
  alpha_raw ~ student_t(7, 0, 2.5);
  te_raw ~ student_t(7, 0, 2.5);
  
  // tau_alpha ~ student_t(7, 0, 10);
  // tau_te ~ student_t(7, 0, 10);
  
  // sigma_te ~ student_t(7, 0, 2.5);
    
  // te_noise ~ normal(0, 1);
    
  // y ~ bernoulli(Phi((alpha + te * treated)./sqrt(1 + square(sigma_te) * treated)));
  y ~ bernoulli(Phi(alpha + te * treated));
    
  // y ~ bernoulli_logit(alpha + (te + te_noise * sigma_te) .* treated);
  // y ~ bernoulli_logit((alpha + te * treated)./(sqrt(1 + square(sigma_te) * treated)));
}

generated quantities {
}
