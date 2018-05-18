functions {
}

data {
  int know_table_A_sample_size;
  
  int<lower = 0> num_obs;
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  int<lower = 1> num_all_treatments; // # of treatment cells
  int<lower = 1> num_all_treatment_coef; // # of linear model coefficients minus intercept
  
  int<lower = 1, upper = num_clusters> num_know_table_A_clusters;
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  
  int<lower = 1, upper = num_obs> know_table_A_cluster_sizes[num_know_table_A_clusters];
  int<lower = 1, upper = num_obs> cluster_pop_sizes[num_clusters];
  
  int<lower = 1, upper = num_obs> num_know_table_A_obs;
  int<lower = 1, upper = num_obs> know_table_A_obs_ids[num_know_table_A_obs];
  int<lower = 0, upper = know_table_A_sample_size> num_know_table_A_recognized[num_know_table_A_obs];
  
  int<lower = 0, upper = know_table_A_sample_size> num_know_table_A_2ord_knows[num_know_table_A_obs];
  int<lower = 0, upper = know_table_A_sample_size> num_know_table_A_2ord_knows_yes[num_know_table_A_obs];
  
  int<lower = 0, upper = know_table_A_sample_size> num_know_table_A_1ord_knows[num_know_table_A_obs];
  int<lower = 0, upper = know_table_A_sample_size> num_know_table_A_1ord_knows_yes[num_know_table_A_obs];
  
  // int<lower = 1, upper = num_obs> num_know_table_B_obs;
  // int<lower = 1, upper = num_obs> know_table_B_obs_ids[num_know_table_B_obs];
  // int<lower = 0, upper = 20> num_know_table_B_recognized[num_know_table_B_obs];
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix; // Design matrix generated from treatment map
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix;
  
  // ATE
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_treatments;
  int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments, 1];  
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_pairs;
  int<lower = 1, upper = num_all_treatments * num_all_treatments> ate_pairs[num_ate_pairs, 2];
  
  int<lower = 0> num_know_table_A_missing_obs_ids;
  int<lower = 0, upper = num_obs> know_table_A_missing_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> know_table_A_missing_obs_ids[num_know_table_A_missing_obs_ids];
  
  int<lower = 0> num_know_table_A_observed_obs_ids;
  int<lower = 0, upper = num_obs> know_table_A_observed_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> know_table_A_observed_obs_ids[num_know_table_A_observed_obs_ids];
}

transformed data {
  int<lower = 1, upper = num_clusters> know_table_A_cluster_id[num_know_table_A_obs] = cluster_id[know_table_A_obs_ids];
  int<lower = 0, upper = num_know_table_A_clusters> know_table_A_cluster_id_dict[num_clusters] = rep_array(0, num_clusters);
  int<lower = 1, upper = num_know_table_A_clusters> know_table_A_cluster_index[num_know_table_A_obs];
  
  vector<lower = 1, upper = num_obs>[num_know_table_A_obs] obs_cluster_pop_size = to_vector(cluster_pop_sizes[know_table_A_cluster_id]);
  real<lower = 0> obs_cluster_pop_size_mean = mean(obs_cluster_pop_size);
  real<lower = 0> obs_cluster_pop_size_sd = sd(obs_cluster_pop_size);
  vector[num_know_table_A_obs] standard_obs_cluster_pop_size = (obs_cluster_pop_size - obs_cluster_pop_size_mean) / obs_cluster_pop_size_sd;
  
  int num_1ord_treatment_coef = num_all_treatment_coef + 1; 
  int num_2ord_treatment_coef = num_all_treatment_coef + 1; 
  
  matrix[num_know_table_A_obs, num_all_treatment_coef] know_table_A_design_matrix = treatment_design_matrix[know_table_A_obs_ids];
  
  matrix[num_know_table_A_obs, num_2ord_treatment_coef] beliefs_design_matrix;  
    
  beliefs_design_matrix[, 1:num_all_treatment_coef] = know_table_A_design_matrix;
  beliefs_design_matrix[, num_all_treatment_coef + 1] = standard_obs_cluster_pop_size;
  
  {
    int cluster_index = 1;
    int curr_cluster_id = 0; 
    
    for (obs_index in 1:num_know_table_A_obs) {
      if (know_table_A_cluster_id[obs_index] != curr_cluster_id) {
        curr_cluster_id = know_table_A_cluster_id[obs_index];
        know_table_A_cluster_id_dict[curr_cluster_id] = cluster_index;
        
        cluster_index += 1;
      }
    }
    
    know_table_A_cluster_index = know_table_A_cluster_id_dict[know_table_A_cluster_id];
  }
}

parameters {
  // Recognizes others
  
  real hyper_recognized_intercept;
  real hyper_recognized_beta_cluster_size;
  
  vector[num_know_table_A_clusters] cluster_recognized_intercept_raw;
  real<lower = 0> cluster_recognized_intercept_sd;
  // real<lower = 3> cluster_recognized_intercept_df;
  
  vector[num_know_table_A_clusters] cluster_recognized_beta_cluster_size_raw;
  real<lower = 0> cluster_recognized_beta_cluster_size_sd;
  // real<lower = 3> cluster_recognized_beta_cluster_size_df;
  
  vector[num_know_table_A_obs] obs_recognized_intercept_raw;
  real<lower = 0> obs_recognized_intercept_sd;
  // real<lower = 3> obs_recognized_intercept_df;
  
  // First order beliefs
  
  vector[num_1ord_treatment_coef] hyper_1ord_beta;
  
  matrix[num_1ord_treatment_coef, num_know_table_A_clusters] cluster_1ord_beta_raw;
  vector<lower = 0>[num_1ord_treatment_coef] cluster_1ord_beta_sd;
  
  matrix[num_1ord_treatment_coef, num_know_table_A_obs] obs_1ord_beta_raw;
  vector<lower = 0>[num_1ord_treatment_coef] obs_1ord_beta_sd;
  
  // Second order beliefs
  
  vector[num_2ord_treatment_coef] hyper_2ord_beta;
  
  // real hyper_2ord_intercept;
  // vector[1] hyper_2ord_beta_treat;
  
  matrix[num_2ord_treatment_coef, num_know_table_A_clusters] cluster_2ord_beta_raw;
  vector<lower = 0>[num_2ord_treatment_coef] cluster_2ord_beta_sd;
  
  matrix[num_2ord_treatment_coef, num_know_table_A_obs] obs_2ord_beta_raw;
  vector<lower = 0>[num_2ord_treatment_coef] obs_2ord_beta_sd;
}

transformed parameters {
  vector[num_know_table_A_clusters] cluster_recognized_intercept = hyper_recognized_intercept + cluster_recognized_intercept_raw * cluster_recognized_intercept_sd;
  vector[num_know_table_A_clusters] cluster_recognized_beta_cluster_size = 
    hyper_recognized_beta_cluster_size + cluster_recognized_beta_cluster_size_raw * cluster_recognized_beta_cluster_size_sd;
    
  vector[num_know_table_A_obs] obs_recognized_intercept = cluster_recognized_intercept[know_table_A_cluster_index] + obs_recognized_intercept_raw * obs_recognized_intercept_sd;
  
  matrix[num_1ord_treatment_coef, num_know_table_A_clusters] cluster_1ord_beta = 
    rep_matrix(hyper_1ord_beta, num_know_table_A_clusters) + diag_matrix(cluster_1ord_beta_sd) * cluster_1ord_beta_raw;
  matrix[num_1ord_treatment_coef, num_know_table_A_obs] obs_1ord_beta = 
    cluster_1ord_beta[, know_table_A_cluster_index] + diag_matrix(obs_1ord_beta_sd) * obs_1ord_beta_raw;
  
  matrix[num_2ord_treatment_coef, num_know_table_A_clusters] cluster_2ord_beta = 
    rep_matrix(hyper_2ord_beta, num_know_table_A_clusters) + diag_matrix(cluster_2ord_beta_sd) * cluster_2ord_beta_raw;
  matrix[num_2ord_treatment_coef, num_know_table_A_obs] obs_2ord_beta = 
    cluster_2ord_beta[, know_table_A_cluster_index] + diag_matrix(obs_2ord_beta_sd) * obs_2ord_beta_raw;
    
  vector[num_know_table_A_obs] recognized_latent_var = 
    obs_recognized_intercept + standard_obs_cluster_pop_size .* cluster_recognized_beta_cluster_size[know_table_A_cluster_index];
    
  vector<lower = 0, upper = 1>[num_know_table_A_obs] rep_know_table_A_prop_recognized = inv_logit(recognized_latent_var);
  vector<lower = 0>[num_know_table_A_obs] degree = 
    to_vector(num_know_table_A_recognized) + (obs_cluster_pop_size - know_table_A_sample_size) .* rep_know_table_A_prop_recognized;
    
  vector[num_know_table_A_obs] beliefs_1ord_latent_var;
  vector<lower = 0, upper = 1>[num_know_table_A_obs] rep_know_table_A_1ord_prop_know;
  vector<lower = 0, upper = 1>[num_know_table_A_obs] beliefs_1ord_prop_know;
  
  vector[num_know_table_A_obs] beliefs_2ord_latent_var;
  vector<lower = 0, upper = 1>[num_know_table_A_obs] rep_know_table_A_2ord_prop_know;
  vector<lower = 0, upper = 1>[num_know_table_A_obs] beliefs_2ord_prop_know;
  
  {
    vector[num_know_table_A_obs] standard_rep_know_table_A_prop_recognized = 
      (rep_know_table_A_prop_recognized - mean(rep_know_table_A_prop_recognized)) / sd(rep_know_table_A_prop_recognized);
      
    beliefs_1ord_latent_var = rows_dot_product(beliefs_design_matrix, obs_1ord_beta');
    rep_know_table_A_1ord_prop_know = inv_logit(beliefs_1ord_latent_var);
    beliefs_1ord_prop_know = 
      (to_vector(num_know_table_A_1ord_knows) + rep_know_table_A_1ord_prop_know .* rep_know_table_A_prop_recognized .* (obs_cluster_pop_size - know_table_A_sample_size)) 
      ./ degree;
      
    beliefs_2ord_latent_var = rows_dot_product(beliefs_design_matrix, obs_2ord_beta');
    rep_know_table_A_2ord_prop_know = inv_logit(beliefs_2ord_latent_var);
    beliefs_2ord_prop_know = 
      (to_vector(num_know_table_A_2ord_knows) + rep_know_table_A_2ord_prop_know .* rep_know_table_A_prop_recognized .* (obs_cluster_pop_size - know_table_A_sample_size)) 
      ./ degree;
  }
}

model {
  // Priors
  
  hyper_recognized_intercept               ~ normal(0, 5);
  hyper_recognized_beta_cluster_size       ~ normal(0, 0.5);
  
  // cluster_recognized_intercept_df          ~ gamma(2, 0.1);
  // cluster_recognized_intercept_raw         ~ student_t(cluster_recognized_intercept_df, 0, 1);
  cluster_recognized_intercept_raw         ~ normal(0, 1);
  cluster_recognized_intercept_sd          ~ normal(0, 0.25);
  
  // cluster_recognized_beta_cluster_size_df  ~ gamma(2, 0.1);
  // cluster_recognized_beta_cluster_size_raw ~ student_t(cluster_recognized_beta_cluster_size_df, 0, 1);
  cluster_recognized_beta_cluster_size_raw ~ normal(0, 1);
  cluster_recognized_beta_cluster_size_sd  ~ normal(0, 0.25);
  
  // obs_recognized_intercept_df              ~ gamma(2, 0.1);
  // obs_recognized_intercept_raw             ~ student_t(obs_recognized_intercept_df, 0, 1);
  obs_recognized_intercept_raw             ~ normal(0, 1);
  obs_recognized_intercept_sd              ~ normal(0, 0.25);
  
  hyper_1ord_beta[1]                         ~ normal(0, 5);
  hyper_1ord_beta[2:num_1ord_treatment_coef] ~ normal(0, 0.5);
  
  to_vector(cluster_1ord_beta_raw)         ~ normal(0, 1);
  cluster_1ord_beta_sd                     ~ normal(0, 0.25);
  
  to_vector(obs_1ord_beta_raw)         ~ normal(0, 1);
  obs_1ord_beta_sd                     ~ normal(0, 0.25);
  
  hyper_2ord_beta[1]                         ~ normal(0, 5);
  hyper_2ord_beta[2:num_2ord_treatment_coef] ~ normal(0, 0.5);
  
  to_vector(cluster_2ord_beta_raw)         ~ normal(0, 1);
  cluster_2ord_beta_sd                     ~ normal(0, 0.25);
  
  to_vector(obs_2ord_beta_raw)         ~ normal(0, 1);
  obs_2ord_beta_sd                     ~ normal(0, 0.25);
  
  // Likelihood 
  
  num_know_table_A_recognized ~ binomial_logit(know_table_A_sample_size, recognized_latent_var); 
  num_know_table_A_1ord_knows ~ binomial_logit(num_know_table_A_recognized, beliefs_1ord_latent_var);
  num_know_table_A_2ord_knows ~ binomial_logit(num_know_table_A_recognized, beliefs_2ord_latent_var);
}

generated quantities {
  vector<lower = 0, upper = 1>[num_ate_treatments] fp_1ord_prop_knows_mean = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] fp_1ord_prop_knows_ate = rep_vector(0, num_ate_pairs);
  
  // vector<lower = 0, upper = 1>[num_ate_treatments] fp_2ord_prop_knows_mean = rep_vector(0, num_ate_treatments);
  vector[num_ate_treatments] fp_2ord_prop_knows_mean = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] fp_2ord_prop_knows_ate = rep_vector(0, num_ate_pairs);
  
  {
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    
    vector[num_know_table_A_obs] rep_num_know_table_A_recognized = rep_know_table_A_prop_recognized .* obs_cluster_pop_size ./ degree;
  
    for (treatment_ids_index in 1:num_ate_treatments) {
      int curr_missing_treatment_size = know_table_A_missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = know_table_A_observed_treatment_sizes[treatment_ids_index];
      int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
      
      int curr_missing_obs_ids[curr_missing_treatment_size] = know_table_A_missing_obs_ids[missing_treatment_pos:missing_treatment_end];
      int curr_observed_obs_ids[curr_observed_treatment_size] = know_table_A_observed_obs_ids[observed_treatment_pos:observed_treatment_end];
      
      vector[curr_missing_treatment_size] missing_1ord_prop_knows = 
        inv_logit(rows_dot_product(append_col(rep_matrix(treatment_map_design_matrix[ate_treatments[treatment_ids_index, 1]], curr_missing_treatment_size),
                                              standard_obs_cluster_pop_size[curr_missing_obs_ids]),
                                   obs_1ord_beta[, curr_missing_obs_ids]')) .* rep_num_know_table_A_recognized[curr_missing_obs_ids];
                                   
      vector[curr_missing_treatment_size] missing_2ord_prop_knows = 
        inv_logit(rows_dot_product(append_col(rep_matrix(treatment_map_design_matrix[ate_treatments[treatment_ids_index, 1]], curr_missing_treatment_size),
                                              standard_obs_cluster_pop_size[curr_missing_obs_ids]),
                                   obs_2ord_beta[, curr_missing_obs_ids]')) .* rep_num_know_table_A_recognized[curr_missing_obs_ids];
                         
      fp_1ord_prop_knows_mean[treatment_ids_index] = (sum(missing_1ord_prop_knows) + sum(beliefs_1ord_prop_know[curr_observed_obs_ids])) 
                                                     / (curr_missing_treatment_size + curr_observed_treatment_size);
                                                     
      fp_2ord_prop_knows_mean[treatment_ids_index] = (sum(missing_2ord_prop_knows) + sum(beliefs_2ord_prop_know[curr_observed_obs_ids])) 
                                                     / (curr_missing_treatment_size + curr_observed_treatment_size);
      
      missing_treatment_pos = missing_treatment_end + 1;
      observed_treatment_pos = observed_treatment_end + 1;
    }
    
    fp_1ord_prop_knows_ate = fp_1ord_prop_knows_mean[ate_pairs[, 1]] - fp_1ord_prop_knows_mean[ate_pairs[, 2]];
    fp_2ord_prop_knows_ate = fp_2ord_prop_knows_mean[ate_pairs[, 1]] - fp_2ord_prop_knows_mean[ate_pairs[, 2]];
  }
}
