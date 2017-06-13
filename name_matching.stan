functions {
}

data {
  int<lower = 0> num_obs;
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  // Constants for priors
  
  real<lower = 0> name_match_false_pos_alpha;
  real<lower = 0> name_match_false_pos_beta;
  
  int <lower = 0, upper = num_obs> num_name_matched;
  int <lower = 0, upper = num_obs> num_name_matching_errors_ids;
  vector<lower = 0, upper = 1>[num_obs] name_matched;
  int<lower = 1, upper = num_obs> name_matched_id[num_name_matched];
  int<lower = 1, upper = num_obs> name_matching_error_ids[num_name_matching_errors_ids];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  
  // Binary deworming outcome 
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
}

transformed data {
  vector<lower = 0, upper = 1>[num_obs] monitored = 1 - name_matched;
}

parameters {
  real<lower = 0, upper = 1> name_match_false_pos; 
  real<lower = 0, upper = 1> name_match_false_neg; 
  
  // real name_match_error_intercept;
  // real name_match_error_diff;
  real<lower = 0, upper = 1> name_match_error_intercept;
}

transformed parameters {
  real<lower = 0, upper = 1> name_match_error_diff =
    (name_match_error_intercept + name_match_error_intercept) * (name_match_false_neg - name_match_false_pos) / (1 - name_match_false_neg + name_match_false_pos);
  // real<lower = 0, upper = 1> name_match_false_neg = 1 - (inv_logit(name_match_error_intercept)/inv_logit(name_match_error_intercept + name_match_error_diff)) + name_match_false_pos;
  // real<lower = 0, upper = 1> name_match_false_neg = 1 - (inv_Phi(name_match_error_intercept)/inv_Phi(name_match_error_intercept + name_match_error_diff)) + name_match_false_pos; 
}

model {
  name_match_false_pos ~ beta(name_match_false_pos_alpha, name_match_false_pos_beta);
  
  dewormed_any[name_matching_error_ids] ~ bernoulli(name_match_error_intercept + name_match_error_diff * (monitored[name_matching_error_ids]));
  // dewormed_any[name_matching_error_ids] ~ bernoulli_logit(name_match_error_intercept + name_match_error_diff * (monitored[name_matching_error_ids]));
  // dewormed_any[name_matching_error_ids] ~ bernoulli(Phi(name_match_error_intercept + name_match_error_diff * (monitored[name_matching_error_ids])));
}

generated quantities {
}
