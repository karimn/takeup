data {
  int<lower = 0> num_obs;
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  int<lower = 1, upper = num_clusters> num_know_table_A_clusters;
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  
  int<lower = 1, upper = num_obs> know_table_A_cluster_sizes[num_know_table_A_clusters];
  int<lower = 1, upper = num_obs> cluster_pop_sizes[num_clusters];
  
  int<lower = 1, upper = num_obs> num_know_table_A_obs;
  int<lower = 1, upper = num_obs> know_table_A_obs_ids[num_know_table_A_obs];
  int<lower = 0, upper = 10> num_know_table_A_recognized[num_know_table_A_obs];
  
  int<lower = 0, upper = 10> num_know_table_A_2ord_knows[num_know_table_A_obs];
  
  // int<lower = 1, upper = num_obs> num_know_table_B_obs;
  // int<lower = 1, upper = num_obs> know_table_B_obs_ids[num_know_table_B_obs];
  // int<lower = 0, upper = 20> num_know_table_B_recognized[num_know_table_B_obs];
}

transformed data {
  int<lower = 1, upper = num_clusters> know_table_A_cluster_id[num_know_table_A_obs] = cluster_id[know_table_A_obs_ids];
  int<lower = 0, upper = num_know_table_A_clusters> know_table_A_cluster_id_dict[num_clusters] = rep_array(0, num_clusters);
  int<lower = 1, upper = num_know_table_A_clusters> know_table_A_cluster_index[num_know_table_A_obs];
  
  vector<lower = 1, upper = num_obs>[num_know_table_A_obs] obs_cluster_pop_size = to_vector(cluster_pop_sizes[know_table_A_cluster_id]);
  real<lower = 0> obs_cluster_pop_size_mean = mean(obs_cluster_pop_size);
  real<lower = 0> obs_cluster_pop_size_sd = sd(obs_cluster_pop_size);
  vector[num_know_table_A_obs] standard_obs_cluster_pop_size = (obs_cluster_pop_size - obs_cluster_pop_size_mean) / obs_cluster_pop_size_sd;
  
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
  
  // Second order beliefs
  
  real hyper_2ord_intercept;
  vector[1] hyper_2ord_beta_treat;
  
  matrix[2, num_know_table_A_clusters] cluster_2ord_beta_raw;
  vector<lower = 0>[2] cluster_2ord_beta_sd;
  
  matrix[2, num_know_table_A_obs] obs_2ord_beta_raw;
  vector<lower = 0>[2] obs_2ord_beta_sd;
}

transformed parameters {
  vector[num_know_table_A_clusters] cluster_recognized_intercept = hyper_recognized_intercept + cluster_recognized_intercept_raw * cluster_recognized_intercept_sd;
  vector[num_know_table_A_clusters] cluster_recognized_beta_cluster_size = 
    hyper_recognized_beta_cluster_size + cluster_recognized_beta_cluster_size_raw * cluster_recognized_beta_cluster_size_sd;
    
  vector[num_know_table_A_obs] obs_recognized_intercept = cluster_recognized_intercept[know_table_A_cluster_index] + obs_recognized_intercept_raw * obs_recognized_intercept_sd;
    
  vector[num_know_table_A_obs] recognized_latent_var = 
    obs_recognized_intercept + standard_obs_cluster_pop_size .* cluster_recognized_beta_cluster_size[know_table_A_cluster_index];
    
  vector<lower = 0, upper = 1>[num_know_table_A_obs] rep_know_table_A_prop_recognized = inv_logit(recognized_latent_var);
  vector<lower = 0>[num_know_table_A_obs]            degree = 
    to_vector(num_know_table_A_recognized) + (obs_cluster_pop_size - 10) .* rep_know_table_A_prop_recognized;
    
  vector[num_know_table_A_obs] beliefs_2ord_latent_var;
  vector<lower = 0, upper = 1>[num_know_table_A_obs] rep_know_table_A_2ord_prop_know;
  
  vector[2] hyper_2ord_beta = append_row(hyper_2ord_intercept, hyper_2ord_beta_treat);
  matrix[2, num_know_table_A_clusters] cluster_2ord_beta = 
    rep_matrix(hyper_2ord_beta, num_know_table_A_clusters) + diag_matrix(cluster_2ord_beta_sd) * cluster_2ord_beta_raw;
  matrix[2, num_know_table_A_obs] obs_2ord_beta = 
    cluster_2ord_beta[, know_table_A_cluster_index] + diag_matrix(obs_2ord_beta_sd) * obs_2ord_beta_raw;
  
  vector<lower = 0, upper = 1>[num_know_table_A_obs] beliefs_2ord_prop_know;
  
  {
    vector[num_know_table_A_obs] standard_rep_know_table_A_prop_recognized = 
      (rep_know_table_A_prop_recognized - mean(rep_know_table_A_prop_recognized)) / sd(rep_know_table_A_prop_recognized);
      
    matrix[num_know_table_A_obs, 2] beliefs_2ord_dm;  
      
    beliefs_2ord_dm[, 1] = rep_vector(1, num_know_table_A_obs);
    beliefs_2ord_dm[, 2] = standard_obs_cluster_pop_size;
      
    // beliefs_2ord_latent_var = rows_dot_product(beliefs_2ord_dm, cluster_2ord_beta[, know_table_A_cluster_index]');
    beliefs_2ord_latent_var = rows_dot_product(beliefs_2ord_dm, obs_2ord_beta');
    rep_know_table_A_2ord_prop_know = inv_logit(beliefs_2ord_latent_var);
    beliefs_2ord_prop_know = (to_vector(num_know_table_A_2ord_knows) + rep_know_table_A_2ord_prop_know .* rep_know_table_A_prop_recognized .* (obs_cluster_pop_size - 10)) 
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
  
  hyper_2ord_intercept                     ~ normal(0, 5);
  hyper_2ord_beta_treat                    ~ normal(0, 0.5);
  
  to_vector(cluster_2ord_beta_raw)         ~ normal(0, 1);
  cluster_2ord_beta_sd                     ~ normal(0, 0.25);
  
  to_vector(obs_2ord_beta_raw)         ~ normal(0, 1);
  obs_2ord_beta_sd                     ~ normal(0, 0.25);
  
  // Likelihood 
  
  num_know_table_A_recognized              ~ binomial_logit(10, recognized_latent_var); 
  
  num_know_table_A_2ord_knows              ~ binomial_logit(num_know_table_A_recognized, beliefs_2ord_latent_var);
}

generated quantities {
}
