functions {
  matrix treatment_deworming_day_rng(matrix missing_latent_var, int use_logit) {
    int curr_treatment_size = rows(missing_latent_var); 
    int num_deworming_days = cols(missing_latent_var);
    
    matrix[curr_treatment_size, num_deworming_days + 1] deworming_days_mask = rep_matrix(0, curr_treatment_size, num_deworming_days + 1); 
    
    for (obs_ids_index in 1:curr_treatment_size) {
      vector[num_deworming_days + 1] deworming_day_prob = rep_vector(1, num_deworming_days + 1);
      vector[num_deworming_days + 1] cumul_prod_alphas = rep_vector(1, num_deworming_days + 1);
    
      real prev_not_deworm_prob = 1; 
      
      for (day_index in 1:num_deworming_days) {
        cumul_prod_alphas[day_index] = prev_not_deworm_prob;
        
        if (use_logit) {
          deworming_day_prob[day_index] = inv_logit(missing_latent_var[obs_ids_index, day_index]);
        } else {
          deworming_day_prob[day_index] = 1 - gumbel_cdf(- missing_latent_var[obs_ids_index, day_index], 0, 1); 
        }
        
        prev_not_deworm_prob = prev_not_deworm_prob * (1 - deworming_day_prob[day_index]);
      }
      
      cumul_prod_alphas[num_deworming_days + 1] = prev_not_deworm_prob;
      
      deworming_day_prob = deworming_day_prob .* cumul_prod_alphas; 
      
      deworming_days_mask[obs_ids_index, categorical_rng(deworming_day_prob)] = 1;
    }
  
    return deworming_days_mask; 
  }
  
  matrix treatment_cell_deworming_day_prop_rng(int[,] treatment_ids, 
                                               int[] missing_obs_ids, 
                                               // int[] missing_cluster_id,
                                               int[] missing_level_id, // Either stratum or cluster
                                               int[] missing_treatment_sizes,
                                               int[] observed_treatment_sizes,
                                               // matrix census_covar_latent_var,
                                               matrix[] latent_var_map, // matrix[num_clusters, num_deworming_days][num_ate_treatments]
                                               int[] observed_dewormed_any,
                                               int use_logit) {
    int num_treatment_ids = size(treatment_ids);
    int num_deworming_days = cols(latent_var_map[1]);
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    matrix[num_treatment_ids, num_deworming_days + 1] cell_deworming_day = rep_matrix(0, num_treatment_ids, num_deworming_days + 1);
    
    matrix[num_deworming_days + 1, num_deworming_days + 1] days_diag = diag_matrix(rep_vector(1, num_deworming_days + 1));
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      
      matrix[curr_missing_treatment_size, num_deworming_days + 1] missing_deworming_days_mask = 
        treatment_deworming_day_rng(latent_var_map[treatment_ids_index, missing_level_id[missing_treatment_pos:missing_treatment_end]], use_logit);
                                    
      if (curr_observed_treatment_size > 0) {
        int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
        
        matrix[curr_observed_treatment_size, num_deworming_days + 1] observed_deworming_days_mask = 
          days_diag[observed_dewormed_any[observed_treatment_pos:observed_treatment_end]];
      
        cell_deworming_day[treatment_ids_index, 1:(num_deworming_days + 1)] = (diagonal(crossprod(missing_deworming_days_mask)) 
          + diagonal(crossprod(observed_deworming_days_mask)))' / (curr_missing_treatment_size + curr_observed_treatment_size);
          
        observed_treatment_pos = observed_treatment_end + 1;
      } else {
        cell_deworming_day[treatment_ids_index, 1:(num_deworming_days + 1)] = 
          diagonal(crossprod(missing_deworming_days_mask))' / curr_missing_treatment_size;
      }
        
      missing_treatment_pos = missing_treatment_end + 1;
    }
    
    return(cell_deworming_day);
  }
}

data {
  int<lower = 0> num_obs;
  int<lower = 1> num_all_treatments; // # of treatment cells
  int<lower = 1> num_all_treatment_coef; // # of linear model coefficients minus intercept
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  // Below references to "maps" refer to a finite set of treatment cells that has a one-to-many relationship to actual observations
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix; // Design matrix generated from treatment map
  int<lower = 1, upper = num_all_treatments> obs_treatment[num_obs]; // ID of observed treatment (from treatment map)
  int<lower = 1> num_subgroups;
  int<lower = 1, upper = num_all_treatment_coef> subgroup_treatment_col_sizes[num_subgroups];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs]; // Observations' stratum IDs
  int<lower = 1, upper = num_strata> cluster_stratum_ids[num_clusters]; // Clusters' stratum IDs
  int<lower = 1, upper = num_clusters> strata_cluster_ids[num_clusters]; // Cluster IDs ordered by stratum
  int<lower = 1, upper = num_obs> cluster_obs_ids[num_obs]; // Observation IDs ordered by cluster
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; // Observations in strata 
  int<lower = 1, upper = num_clusters> strata_num_clusters[num_strata]; // Clusters in strata
  int<lower = 1, upper = num_obs> cluster_sizes[num_clusters]; // Observations in clusters
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  
  // Deworming outcomes
  
  int num_deworming_days;
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
  int<lower = 1, upper = num_deworming_days + 1> dewormed_day_any[num_obs];
 
  int num_dewormed; 
  int strata_dewormed_sizes[num_strata];
  int dewormed_ids[num_dewormed];
  int<lower = 1> stratum_dewormed_index[num_dewormed];
  
  // Covariates
  
  int<lower = 1> num_census_covar_coef;
  matrix[num_obs, num_census_covar_coef] census_covar_dm;
  
  // Dynamics
  
  int relevant_latent_var_map[num_deworming_days + 1, num_deworming_days];
 
  // Counterfactuals and ATE
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_treatments;
  int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments, 2]; 
  
  int<lower = 0> num_missing_obs_ids;
  int<lower = 0, upper = num_obs> missing_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> missing_obs_ids[num_missing_obs_ids];
  
  int<lower = 0> num_observed_obs_ids;
  int<lower = 0, upper = num_obs> observed_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> observed_obs_ids[num_observed_obs_ids];
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_pairs;
  
  int<lower = 1, upper = num_all_treatments * num_all_treatments> ate_pairs[num_ate_pairs, 2];
  
  int missing_treatment_stratum_id[num_missing_obs_ids]; 
  int missing_treatment_cluster_id[num_missing_obs_ids]; 
  
  // matrix[num_missing_obs_ids, num_census_covar_coef] missing_census_covar_dm;
  
  int observed_dewormed_day[num_observed_obs_ids];
  
  // Constants for hyperpriors 
  
  real<lower = 0> scale_sigma;
  
  real<lower = 0> hyper_coef_sigma;
  real<lower = 0> hyper_intercept_sigma;
  
  real<lower = 0> lkj_df;
  
  // Configuration
  
  int<lower = 0, upper = 1> estimate_ate;
  int<lower = 0, upper = 1> use_logit; 
  
  int<lower = 0> param_poly_order; 
}

transformed data {
  real delta = 1e-9; 
  
  int num_non_phone_treat_col = subgroup_treatment_col_sizes[1];
  int num_phone_treat_col = subgroup_treatment_col_sizes[2];
  
  int non_phone_treat_col[num_non_phone_treat_col];
  int phone_treat_col[num_phone_treat_col];
  
  vector[num_deworming_days - 1] mu_dyn_zero = rep_vector(0, num_deworming_days - 1);
  vector[num_all_treatment_coef] mu_treat_zero = rep_vector(0, num_all_treatment_coef);
  row_vector[num_all_treatment_coef] mu_treat_zero_row = rep_row_vector(0, num_all_treatment_coef);
  matrix[num_deworming_days - 1, num_deworming_days - 1] dyn_identity_mat = diag_matrix(rep_vector(1, num_deworming_days - 1));
  matrix[num_deworming_days, num_deworming_days] treat_identity_mat = diag_matrix(rep_vector(1, num_deworming_days));
  
  int<lower = 1, upper = num_obs * num_deworming_days> num_relevant_obs_days = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any]));
  
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  
  matrix[num_relevant_obs_days, num_census_covar_coef] census_covar_dm_long;
  matrix[num_relevant_obs_days, num_all_treatment_coef] treatment_design_matrix_long;
  int cluster_id_long[num_relevant_obs_days];
  
  int<lower = 0, upper = 1> relevant_dewormed_any_daily[num_relevant_obs_days] = rep_array(0, num_relevant_obs_days);
  int<lower = 1, upper = num_deworming_days> dewormed_day_long[num_relevant_obs_days];
  int<lower = 0> relevant_daily_strata_sizes[num_strata];
  int<lower = 0> relevant_daily_cluster_sizes[num_clusters]; // By cluster ID
  
  int<lower = 1> actual_poly_order = max(1, param_poly_order);
  
  matrix[num_deworming_days, actual_poly_order] days_poly_trend;
  matrix[num_deworming_days, actual_poly_order] Q_days_poly_trend;
  matrix[actual_poly_order, actual_poly_order] R_inv_days_poly_trend;
  int<lower = 1> num_dyn_rows = param_poly_order > 0 ? param_poly_order : num_deworming_days - 1;
  
  for (col_index in 1:num_non_phone_treat_col) non_phone_treat_col[col_index] = col_index;
  for (col_index in 1:num_phone_treat_col) phone_treat_col[col_index] = num_non_phone_treat_col + col_index;
  
  if (param_poly_order > 0) {
    matrix[param_poly_order, param_poly_order] R_days_poly_trend;
    vector[num_deworming_days] poly_1_trend = cumulative_sum(rep_vector(1, num_deworming_days)) - 1;
    
    days_poly_trend = rep_matrix(poly_1_trend / sd(poly_1_trend), actual_poly_order);  // scale trend
    // days_poly_trend = rep_matrix((poly_1_trend - mean(poly_1_trend)) / sd(poly_1_trend), actual_poly_order);  // Standardize trend
    
    for (poly_index in 2:param_poly_order) {
      days_poly_trend[, poly_index] = days_poly_trend[, poly_index] .* days_poly_trend[, poly_index - 1];
    }
    
    Q_days_poly_trend = qr_Q(days_poly_trend)[, 1:param_poly_order] * num_deworming_days; 
    R_days_poly_trend = qr_R(days_poly_trend)[1:param_poly_order, ] / num_deworming_days; 
    R_inv_days_poly_trend = inv(R_days_poly_trend);
  }
  
  {
    int stratum_pos = 1;
    int cluster_pos = 1;
    int cluster_obs_pos = 1;
    int relevant_day_obs_pos = 1;
   
    for (stratum_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[stratum_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;

      for (cluster_pos_index in cluster_pos:cluster_end) {
        int curr_cluster_id = strata_cluster_ids[cluster_pos_index];
        int curr_cluster_num_obs = cluster_sizes[curr_cluster_id];
        
        int cluster_obs_end = cluster_obs_pos + curr_cluster_num_obs - 1;
        
        relevant_daily_cluster_sizes[curr_cluster_id] = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any[cluster_obs_pos:cluster_obs_end]]));
        
        cluster_obs_pos = cluster_obs_end + 1;
      }

      relevant_daily_strata_sizes[stratum_index] = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any[stratum_pos:stratum_end]]));

      stratum_pos = stratum_end + 1;
      cluster_pos = cluster_end + 1;
    }

    if (sum(relevant_daily_strata_sizes) != num_relevant_obs_days) {
      reject("Unexpected relevant strata sizes ", sum(relevant_daily_strata_sizes), " ", num_relevant_obs_days);
    }
    
    for (obs_index in 1:num_obs) {
      int relevant_latent_var_array[num_deworming_days] = relevant_latent_var_map[dewormed_day_any[obs_index]];
      int curr_relevant_days = sum(relevant_latent_var_array);
      int relevant_day_obs_end = relevant_day_obs_pos + curr_relevant_days - 1;
      
      census_covar_dm_long[relevant_day_obs_pos:relevant_day_obs_end] = rep_matrix(census_covar_dm[obs_index], curr_relevant_days);
      
      treatment_design_matrix_long[relevant_day_obs_pos:relevant_day_obs_end] = rep_matrix(treatment_map_design_matrix[obs_treatment[obs_index]], curr_relevant_days);
        
      cluster_id_long[relevant_day_obs_pos:relevant_day_obs_end] = rep_array(cluster_id[obs_index], curr_relevant_days);
      
      for (obs_rel_day_index in 1:curr_relevant_days) {
        int curr_obs_day_pos = relevant_day_obs_pos + obs_rel_day_index - 1;
        
        dewormed_day_long[curr_obs_day_pos] = obs_rel_day_index;
      }
     
      if (dewormed_any[obs_index]) {
        relevant_dewormed_any_daily[relevant_day_obs_pos + dewormed_day_any[obs_index] - 1] = 1;
      }
      
      relevant_day_obs_pos = relevant_day_obs_end + 1;
    }
  }
}

parameters {
  row_vector[num_all_treatment_coef] hyper_beta_day1; 
  vector[num_deworming_days - 1] hyper_baseline_nonparam_dyn_effect;
  matrix[num_dyn_rows, num_all_treatment_coef - 1] hyper_treat_beta_dyn_effect;
  
  // matrix[num_deworming_days - 1, num_nonparam_col] hyper_beta_dyn_effect; 
  
  vector<lower = 0>[num_all_treatment_coef] strata_beta_day1_tau;
  cholesky_factor_corr[num_non_phone_treat_col] strata_beta_day1_L_corr_mat_non_phone;
  cholesky_factor_corr[num_phone_treat_col] strata_beta_day1_L_corr_mat_phone;
  
  // matrix[num_strata, num_all_treatment_coef] strata_beta_day1;
  matrix[num_all_treatment_coef, num_strata] strata_beta_day1_raw;
  
  // vector<lower = 0>[num_deworming_days - 1] strata_beta_dyn_tau_raw[num_all_treatment_coef];
  // cholesky_factor_corr[num_deworming_days - 1] strata_beta_dyn_L_corr_mat[num_all_treatment_coef];
  // matrix[num_deworming_days - 1, num_all_treatment_coef] strata_beta_dyn_effect[num_strata]; 
  
  // vector<lower = 0>[num_all_treatment_coef] cluster_beta_day1_tau;
  vector<lower = 0>[num_all_treatment_coef] cluster_beta_day1_tau;
  // cholesky_factor_corr[num_non_phone_treat_col] cluster_beta_day1_L_corr_mat_non_phone;
  // cholesky_factor_corr[num_phone_treat_col] cluster_beta_day1_L_corr_mat_phone;
  cholesky_factor_corr[num_phone_treat_col] cluster_beta_day1_L_corr_mat;
  
  matrix[num_clusters, num_all_treatment_coef] cluster_beta_day1;
  // matrix[num_clusters, num_phone_treat_col] cluster_beta_day1;
}

transformed parameters {
  matrix[num_dyn_rows, num_all_treatment_coef - 1] QR_hyper_treat_beta_dyn_effect;
  
  matrix[num_deworming_days, num_all_treatment_coef] hyper_full_beta_dyn_effect = rep_matrix(0, num_deworming_days, num_all_treatment_coef); 
  // matrix[num_deworming_days, num_all_treatment_coef] QR_hyper_full_beta_dyn_effect; 
  
  matrix[num_deworming_days, num_all_treatment_coef] hyper_beta = rep_matrix(0, num_deworming_days, num_all_treatment_coef); 
  // matrix[num_deworming_days, num_all_treatment_coef] hyper_full_beta_dyn_effect = append_row(mu_treat_zero_row, hyper_beta_dyn_effect); 
  // matrix[num_deworming_days, num_all_treatment_coef] hyper_beta = rep_matrix(hyper_beta_day1, num_deworming_days) + hyper_full_beta_dyn_effect; 
  
  matrix[num_strata, num_all_treatment_coef] strata_beta_day1;
  
  matrix[num_deworming_days, num_all_treatment_coef] strata_beta[num_strata];
  
  matrix[num_deworming_days, num_all_treatment_coef] cluster_beta[num_clusters];
  // matrix[num_deworming_days, num_all_treatment_coef] cluster_beta[num_clusters];
   
  // for (beta_index in 1:num_all_treatment_coef) {
  //   strata_beta_dyn_tau[beta_index] = strata_beta_dyn_tau_raw[beta_index] * scale_sigma;
  // }
 
  hyper_full_beta_dyn_effect[2:num_deworming_days, 1] = hyper_baseline_nonparam_dyn_effect; 
  
  if (param_poly_order > 0) {
    // hyper_full_beta_dyn_effect[2:num_deworming_days, 2:num_all_treatment_coef] = days_poly_trend * hyper_treat_beta_dyn_effect;
    hyper_full_beta_dyn_effect[, 2:num_all_treatment_coef] = Q_days_poly_trend * hyper_treat_beta_dyn_effect;
    QR_hyper_treat_beta_dyn_effect = R_inv_days_poly_trend * hyper_treat_beta_dyn_effect;
    
    // QR_hyper_full_beta_dyn_effect[2:num_deworming_days, 1] = hyper_baseline_nonparam_dyn_effect; 
    // QR_hyper_full_beta_dyn_effect[, 2:num_all_treatment_coef] = days_poly_trend * QR_hyper_treat_beta_dyn_effect;
    
    // hyper_beta = rep_matrix(hyper_beta_day1, num_deworming_days) + QR_hyper_full_beta_dyn_effect;
  } else {
    hyper_full_beta_dyn_effect[2:num_deworming_days, 2:num_all_treatment_coef] = hyper_treat_beta_dyn_effect;
  }
  
  hyper_beta = rep_matrix(hyper_beta_day1, num_deworming_days) + hyper_full_beta_dyn_effect;

  {
    int cluster_pos = 1;
    
    matrix[num_all_treatment_coef, num_all_treatment_coef] strata_beta_day1_L_vcov = rep_matrix(0, num_all_treatment_coef, num_all_treatment_coef);
    
    strata_beta_day1_L_vcov[non_phone_treat_col, non_phone_treat_col] = diag_pre_multiply(strata_beta_day1_tau[non_phone_treat_col], 
                                                                                          strata_beta_day1_L_corr_mat_non_phone);
    strata_beta_day1_L_vcov[phone_treat_col, phone_treat_col] = diag_pre_multiply(strata_beta_day1_tau[phone_treat_col], strata_beta_day1_L_corr_mat_phone);
    // matrix[7, 7] strata_beta_day1_L_vcov_non_phone = diag_pre_multiply(strata_beta_day1_tau[1:7], strata_beta_day1_L_corr_mat_non_phone);
    // matrix[19, 19] strata_beta_day1_L_vcov_phone = diag_pre_multiply(strata_beta_day1_tau[8:26], strata_beta_day1_L_corr_mat_phone);
    
    strata_beta_day1 = rep_matrix(hyper_beta_day1, num_strata) + (strata_beta_day1_L_vcov * strata_beta_day1_raw)';
    
    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
     
      strata_beta[stratum_index] = rep_matrix(strata_beta_day1[stratum_index], num_deworming_days) + hyper_full_beta_dyn_effect;
      // strata_beta[stratum_index, 2:num_deworming_days] = strata_beta[stratum_index, 2:num_deworming_days] + hyper_beta_dyn_effect;
    
      for (cluster_index in cluster_pos:cluster_end) {
        int curr_cluster_id = strata_cluster_ids[cluster_index];
        
        // cluster_beta[curr_cluster_id, , non_phone_treat_col] = strata_beta[stratum_index, , non_phone_treat_col];
        // cluster_beta[curr_cluster_id, , phone_treat_col] = rep_matrix(cluster_beta_day1[curr_cluster_id], num_deworming_days);
        // cluster_beta[curr_cluster_id, 2:num_deworming_days, phone_treat_col] = 
        //   cluster_beta[curr_cluster_id, 2:num_deworming_days, phone_treat_col] + hyper_beta_dyn_effect[, phone_treat_col];
        
        cluster_beta[curr_cluster_id] = rep_matrix(cluster_beta_day1[curr_cluster_id], num_deworming_days) + hyper_full_beta_dyn_effect;
      }
      
      cluster_pos = cluster_end + 1;
    }
  }
}

model {
  hyper_beta_day1[1] ~ normal(0, hyper_intercept_sigma);
  
  hyper_beta_day1[2:num_all_treatment_coef] ~ normal(0, hyper_coef_sigma);
  
  hyper_baseline_nonparam_dyn_effect ~ normal(0, hyper_coef_sigma);
  
  if (param_poly_order > 0) {
    to_vector(QR_hyper_treat_beta_dyn_effect) ~ normal(0, hyper_coef_sigma);
  } else {
    to_vector(hyper_treat_beta_dyn_effect) ~ normal(0, hyper_coef_sigma);
  }
  
  strata_beta_day1_tau ~ normal(0, scale_sigma);
  
  strata_beta_day1_L_corr_mat_non_phone ~ lkj_corr_cholesky(lkj_df);
  strata_beta_day1_L_corr_mat_phone ~ lkj_corr_cholesky(lkj_df);
  
  cluster_beta_day1_tau ~ normal(0, scale_sigma);
  // cluster_beta_day1_L_corr_mat_non_phone ~ lkj_corr_cholesky(lkj_df);
  // cluster_beta_day1_L_corr_mat_phone ~ lkj_corr_cholesky(lkj_df);
  cluster_beta_day1_L_corr_mat ~ lkj_corr_cholesky(lkj_df);
  
  to_vector(strata_beta_day1_raw) ~ normal(0, 1);
  
  {
    vector[num_relevant_obs_days] latent_var = rep_vector(0, num_relevant_obs_days);
    
    int cluster_pos = 1;
    int relevant_daily_cluster_pos = 1;
    
    matrix[num_all_treatment_coef, num_all_treatment_coef] cluster_beta_day1_L_vcov = rep_matrix(0, num_all_treatment_coef, num_all_treatment_coef);
    // matrix[num_phone_treat_col, num_phone_treat_col] cluster_beta_day1_L_vcov = rep_matrix(0, num_phone_treat_col, num_phone_treat_col); 
    
    cluster_beta_day1_L_vcov[non_phone_treat_col, non_phone_treat_col] = diag_matrix(cluster_beta_day1_tau[non_phone_treat_col]);
    cluster_beta_day1_L_vcov[phone_treat_col, phone_treat_col] = diag_pre_multiply(cluster_beta_day1_tau[phone_treat_col], cluster_beta_day1_L_corr_mat);
    // cluster_beta_day1_L_vcov[1:7, 1:7] = diag_pre_multiply(cluster_beta_day1_tau[1:7], cluster_beta_day1_L_corr_mat_non_phone);
    // cluster_beta_day1_L_vcov[8:26, 8:26] = diag_pre_multiply(cluster_beta_day1_tau[8:26], cluster_beta_day1_L_corr_mat_phone); 
    // matrix[7, 7] cluster_beta_day1_L_vcov_non_phone = diag_pre_multiply(cluster_beta_day1_tau[1:7], cluster_beta_day1_L_corr_mat_non_phone);
    // matrix[19, 19] cluster_beta_day1_L_vcov_phone = diag_pre_multiply(cluster_beta_day1_tau[8:26], cluster_beta_day1_L_corr_mat_phone); 

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;

      for (cluster_pos_index in cluster_pos:cluster_end) {
        int curr_cluster_id = strata_cluster_ids[cluster_pos_index];
        int curr_cluster_num_obs = cluster_sizes[curr_cluster_id];
        
        int curr_relevant_daily_cluster_size = relevant_daily_cluster_sizes[curr_cluster_id];
        int relevant_daily_cluster_end = relevant_daily_cluster_pos + curr_relevant_daily_cluster_size - 1;

        // cluster_beta_day1[curr_cluster_id] ~ multi_normal_cholesky(strata_beta_day1[stratum_index, phone_treat_col], cluster_beta_day1_L_vcov);
        cluster_beta_day1[curr_cluster_id] ~ multi_normal_cholesky(strata_beta_day1[stratum_index], cluster_beta_day1_L_vcov);
        // cluster_beta_day1[curr_cluster_id, 1:7] ~ multi_normal_cholesky(strata_beta_day1[stratum_index, 1:7], cluster_beta_day1_L_vcov_non_phone);
        // cluster_beta_day1[curr_cluster_id, 8:26] ~ multi_normal_cholesky(strata_beta_day1[stratum_index, 8:26], cluster_beta_day1_L_vcov_phone);

        latent_var[relevant_daily_cluster_pos:relevant_daily_cluster_end] =
          rows_dot_product(treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end],
                           cluster_beta[curr_cluster_id, dewormed_day_long[relevant_daily_cluster_pos:relevant_daily_cluster_end]]);
                           
        relevant_daily_cluster_pos = relevant_daily_cluster_end + 1;
      }
    
      // latent_var[relevant_daily_stratum_pos:relevant_daily_stratum_end] =
      //   rows_dot_product(treatment_design_matrix_long[relevant_daily_stratum_pos:relevant_daily_stratum_end],
      //                    strata_beta[stratum_index, dewormed_day_long[relevant_daily_stratum_pos:relevant_daily_stratum_end]]);

      cluster_pos = cluster_end + 1;
    }
   
    //  Non-multilevel version 
    // vector[num_relevant_obs_days] latent_var = rows_dot_product(treatment_design_matrix_long, hyper_beta[dewormed_day_long]);
    
    if (use_logit) { 
      relevant_dewormed_any_daily ~ bernoulli_logit(latent_var);
    } else {
      relevant_dewormed_any_daily ~ bernoulli(inv_cloglog(latent_var));
    }
  }
}

generated quantities {
  matrix<lower = 0>[num_ate_treatments, num_deworming_days + 1] est_deworming_days = rep_matrix(0, num_ate_treatments, num_deworming_days + 1);
  vector<lower = 0, upper = 1>[num_ate_treatments] est_takeup = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] est_takeup_ate = rep_vector(0, num_ate_pairs);
  
  if (estimate_ate) {
    int stratum_pos = 1;

    matrix[num_clusters, num_ate_treatments] cluster_latent_var_map_day1 = cluster_beta_day1 * treatment_map_design_matrix[ate_treatments[, 1]]';
      // append_col(strata_beta_day1[cluster_stratum_ids], cluster_beta_day1) * treatment_map_design_matrix[ate_treatments[, 1]]';
      
    matrix[num_clusters, num_deworming_days] cluster_latent_var_map[num_ate_treatments];

    for (ate_treat_index in 1:num_ate_treatments) {
      row_vector[num_deworming_days] ate_hyper_dyn_effect;
      
      ate_hyper_dyn_effect = treatment_map_design_matrix[ate_treatments[ate_treat_index, 2]] * (hyper_full_beta_dyn_effect');
      
      cluster_latent_var_map[ate_treat_index] = 
        rep_matrix(cluster_latent_var_map_day1[, ate_treat_index], num_deworming_days) + rep_matrix(ate_hyper_dyn_effect, num_clusters);
      
      // cluster_latent_var_map[ate_treat_index, , 2:num_deworming_days] = 
      //   cluster_latent_var_map[ate_treat_index, , 2:num_deworming_days] + rep_matrix(ate_hyper_dyn_effect, num_clusters);
    }

    est_deworming_days =
      treatment_cell_deworming_day_prop_rng(ate_treatments,
                                            missing_obs_ids,
                                            missing_treatment_cluster_id,
                                            missing_treatment_sizes,
                                            observed_treatment_sizes,
                                            cluster_latent_var_map,
                                            observed_dewormed_day,
                                            use_logit);

    est_takeup = 1 - est_deworming_days[, 13];
    est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }
}