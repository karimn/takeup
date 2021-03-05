functions {
}

data {
  int<lower = 0, upper = 1> use_obs_level;
  int<lower = 0, upper = 1> use_cluster_level;
  int<lower = 0, upper = 1> use_stratum_level;
  
  int know_table_A_sample_size;
  
  int<lower = 0> num_obs;
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  int<lower = 1> num_treatments; // # of treatment cells
  
  int<lower = 1, upper = num_clusters> cluster_index[num_obs];
  int<lower = 1, upper = num_strata> stratum_index[num_obs];
  
  int<lower = 0, upper = know_table_A_sample_size> num_recognized[num_obs];
  int<lower = 0, upper = know_table_A_sample_size> num_knows_1ord[num_obs];
  int<lower = 0, upper = know_table_A_sample_size> num_knows_2ord[num_obs];
  
  int<lower = 1, upper = num_treatments> treatment_id[num_obs];
  
  matrix[num_treatments, num_treatments] treatment_map_design_matrix;
 
  int<lower = 0> num_ate_pairs; 
  int<lower = 1, upper = num_treatments> ate_pairs[num_ate_pairs, 2];
}

transformed data {
  matrix[num_obs, num_treatments] treatment_design_matrix = treatment_map_design_matrix[treatment_id];
}

parameters {
  row_vector[num_treatments] hyper_beta_1ord;
  
  matrix[use_stratum_level ? num_strata : 0, num_treatments] stratum_beta_1ord_raw;
  matrix[use_cluster_level ? num_clusters : 0, num_treatments] cluster_beta_1ord_raw;
  matrix[use_obs_level ? num_obs : 0, num_treatments] obs_beta_1ord_raw;
  
  row_vector<lower = 0>[use_stratum_level ? num_treatments : 0] stratum_beta_1ord_sd;
  row_vector<lower = 0>[use_cluster_level ? num_treatments : 0] cluster_beta_1ord_sd;
  row_vector<lower = 0>[use_obs_level ? num_treatments : 0] obs_beta_1ord_sd;
  
  row_vector[num_treatments] hyper_beta_2ord;
  
  matrix[use_stratum_level? num_strata : 0, num_treatments] stratum_beta_2ord_raw;
  matrix[use_cluster_level ? num_clusters : 0, num_treatments] cluster_beta_2ord_raw;
  matrix[use_obs_level ? num_obs : 0, num_treatments] obs_beta_2ord_raw;
  
  row_vector<lower = 0>[use_stratum_level ? num_treatments : 0] stratum_beta_2ord_sd;
  row_vector<lower = 0>[use_cluster_level ? num_treatments: 0] cluster_beta_2ord_sd;
  row_vector<lower = 0>[use_obs_level ? num_treatments : 0] obs_beta_2ord_sd;
  
  vector[num_obs] obs_beta_common_raw;
  real<lower = 0> obs_beta_common_sd;
}

transformed parameters {
  matrix[num_obs, num_treatments] centered_obs_beta_1ord; 
  matrix[num_obs, num_treatments] centered_obs_beta_2ord;  
  
  {
    matrix[num_strata, num_treatments] stratum_beta_1ord;
    matrix[num_clusters, num_treatments] cluster_beta_1ord;
    matrix[num_obs, num_treatments] obs_beta_1ord;
    
    matrix[num_strata, num_treatments] stratum_beta_2ord;
    matrix[num_clusters, num_treatments] cluster_beta_2ord;
    matrix[num_obs, num_treatments] obs_beta_2ord;
    
    vector[num_obs] obs_beta_common = obs_beta_common_raw * obs_beta_common_sd;
    
    stratum_beta_1ord = use_stratum_level ? stratum_beta_1ord_raw .* rep_matrix(stratum_beta_1ord_sd, num_strata) : rep_matrix(0, num_strata, num_treatments);
    cluster_beta_1ord = use_cluster_level ? cluster_beta_1ord_raw .* rep_matrix(cluster_beta_1ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_treatments);
    obs_beta_1ord = use_obs_level ? obs_beta_1ord_raw .* rep_matrix(obs_beta_1ord_sd, num_obs) : rep_matrix(0, num_obs, num_treatments);
    
    stratum_beta_2ord = use_stratum_level ? stratum_beta_2ord_raw .* rep_matrix(stratum_beta_2ord_sd, num_strata) : rep_matrix(0, num_strata, num_treatments);
    cluster_beta_2ord = use_cluster_level ? cluster_beta_2ord_raw .* rep_matrix(cluster_beta_2ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_treatments);
    obs_beta_2ord = use_obs_level ? obs_beta_2ord_raw .* rep_matrix(obs_beta_2ord_sd, num_obs) : rep_matrix(0, num_obs, num_treatments);
    
    centered_obs_beta_1ord = 
      rep_matrix(hyper_beta_1ord, num_obs) + stratum_beta_1ord[stratum_index] + cluster_beta_1ord[cluster_index] + obs_beta_1ord + rep_matrix(obs_beta_common, num_treatments);
    
    centered_obs_beta_2ord = 
      rep_matrix(hyper_beta_2ord, num_obs) + stratum_beta_2ord[stratum_index] + cluster_beta_2ord[cluster_index] + obs_beta_2ord + rep_matrix(obs_beta_common, num_treatments);
  }
}

model {
  obs_beta_common_raw ~ std_normal();
  obs_beta_common_sd ~ normal(0, 0.125);
  
  hyper_beta_1ord[1] ~ normal(0, 2);
  hyper_beta_1ord[2:] ~ normal(0, 1);
 
  if (use_obs_level) { 
    to_vector(obs_beta_1ord_raw) ~ std_normal();
    obs_beta_1ord_sd ~ normal(0, 0.125);
  }
  
  if (use_cluster_level) {
    to_vector(cluster_beta_1ord_raw) ~ std_normal();
    cluster_beta_1ord_sd ~ normal(0, 0.25);
  }
  
  if (use_stratum_level) {
    to_vector(stratum_beta_1ord_raw) ~ std_normal();
    stratum_beta_1ord_sd ~ normal(0, 0.5);
  }
  
  num_knows_1ord ~ binomial_logit(num_recognized, rows_dot_product(treatment_design_matrix, centered_obs_beta_1ord));
  
  hyper_beta_2ord[1] ~ normal(0, 2);
  hyper_beta_2ord[2:] ~ normal(0, 1);
  
  if (use_obs_level) { 
    to_vector(obs_beta_2ord_raw) ~ std_normal();
    obs_beta_2ord_sd ~ normal(0, 0.125);
  }
  
  if (use_cluster_level) {
    to_vector(cluster_beta_2ord_raw) ~ std_normal();
    cluster_beta_2ord_sd ~ normal(0, 0.25);
  }
  
  if (use_stratum_level) {
    to_vector(stratum_beta_2ord_raw) ~ std_normal();
    stratum_beta_2ord_sd ~ normal(0, 0.5);
  }
  
  num_knows_2ord ~ binomial_logit(num_recognized, rows_dot_product(treatment_design_matrix, centered_obs_beta_2ord));
}

generated quantities {
  matrix<lower = 0, upper = 1>[num_obs, num_treatments] obs_prob_1ord = inv_logit(centered_obs_beta_1ord * treatment_map_design_matrix');
  vector<lower = -1, upper = 1>[num_treatments] prob_1ord;
  vector<lower = -1, upper = 1>[num_ate_pairs] ate_1ord;
  
  matrix<lower = 0, upper = 1>[num_obs, num_treatments] obs_prob_2ord = inv_logit(centered_obs_beta_2ord * treatment_map_design_matrix');
  vector<lower = -1, upper = 1>[num_treatments] prob_2ord;
  vector<lower = -1, upper = 1>[num_ate_pairs] ate_2ord;
  
  for (treatment_index in 1:num_treatments) {
    prob_1ord[treatment_index] = mean(obs_prob_1ord[, treatment_index]);
    prob_2ord[treatment_index] = mean(obs_prob_2ord[, treatment_index]);
  }
  
  ate_1ord = prob_1ord[ate_pairs[, 1]] - prob_1ord[ate_pairs[, 2]];
  ate_2ord = prob_2ord[ate_pairs[, 1]] - prob_2ord[ate_pairs[, 2]];
}
