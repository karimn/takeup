functions {
  matrix treatment_deworming_day_rng(//matrix census_covar_latent_var, 
                                     //matrix day_constant_treat_latent_var,
                                     //row_vector dynamic_latent_var, 
                                     //vector cluster_effects) {
                                     matrix missing_latent_var,
                                     int use_logit) {
    int curr_treatment_size = rows(missing_latent_var); 
    // int curr_treatment_size = rows(census_covar_latent_var); 
    int num_deworming_days = cols(missing_latent_var);
    // int num_deworming_days = cols(census_covar_latent_var);
    
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
  
  /**
   * @return A (number of treatments) X (number of deworming days + 1) matrix. 
   */
  matrix treatment_cell_deworming_day_prop_rng(int[,] treatment_ids, 
                                               int[] missing_obs_ids, 
                                               int[] missing_stratum_id,
                                               int[] missing_cluster_id,
                                               int[] missing_treatment_sizes,
                                               int[] observed_treatment_sizes,
                                               // vector cluster_effects,
                                               // matrix census_covar_latent_var,
                                               matrix[] latent_var_map, // matrix[num_strata, num_deworming_days][num_ate_treatments]
                                               // matrix dynamic_latent_var,
                                               int[] observed_dewormed_any,
                                               int use_logit) {
    int num_treatment_ids = size(treatment_ids);
    int num_deworming_days = cols(latent_var_map[1]);
    // int num_deworming_days = cols(census_covar_latent_var);
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    matrix[num_treatment_ids, num_deworming_days + 1] cell_deworming_day = rep_matrix(0, num_treatment_ids, num_deworming_days + 1);
    
    matrix[num_deworming_days + 1, num_deworming_days + 1] days_diag = diag_matrix(rep_vector(1, num_deworming_days + 1));
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      // int curr_all_treatment_id = treatment_ids[treatment_ids_index];
      // int curr_dyn_treatment_id = all_treat_dyn_treat_id[treatment_ids_index];
      
      matrix[curr_missing_treatment_size, num_deworming_days + 1] missing_deworming_days_mask = 
        treatment_deworming_day_rng(latent_var_map[treatment_ids_index, missing_stratum_id[missing_treatment_pos:missing_treatment_end]], use_logit);
          
          // census_covar_latent_var[missing_obs_ids[missing_treatment_pos:missing_treatment_end]],
          //                           rep_matrix(day_constant_treat_latent_var[treatment_ids_index, missing_stratum_id[missing_treatment_pos:missing_treatment_end]]', num_deworming_days),
          //                           dynamic_latent_var[curr_dyn_treatment_id],
          //                           cluster_effects[missing_cluster_id[missing_treatment_pos:missing_treatment_end]]);
                                    
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
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 
  
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
  
  // int<lower = 0, upper = num_all_treatments> num_static_ate_treatments;
  // int<lower = 0, upper = num_all_treatments> num_dyn_ate_treatments;
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_treatments;
  // int<lower = 1, upper = num_all_treatments> static_ate_treatments[num_static_ate_treatments]; 
  // int<lower = 1, upper = num_all_treatments> dyn_ate_treatments[num_dyn_ate_treatments]; 
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
  
  // real<lower = 0> scale_df;
  real<lower = 0> scale_sigma;
  
  // real<lower = 0> coef_df;
  real<lower = 0> hyper_coef_sigma;
  real<lower = 0> hyper_intercept_sigma;
  
  real<lower = 0> lkj_df;
  
  // Configuration
  
  int<lower = 0, upper = 1> estimate_ate;
  int<lower = 0, upper = 1> use_logit; 
  
  int<lower = 1, upper = num_deworming_days> first_dynamics_day;
}

transformed data {
  real delta = 1e-9; 
  // matrix<lower = 0>[num_deworming_days - 1, num_all_treatment_coef] scale_sigma_mat = rep_matrix(scale_sigma, num_deworming_days - 1, num_all_treatment_coef);
  
  vector[num_deworming_days - 1] mu_dyn_zero = rep_vector(0, num_deworming_days - 1);
  matrix[num_deworming_days - 1, num_deworming_days - 1] dyn_identity_mat = diag_matrix(rep_vector(1, num_deworming_days - 1));
  
  int<lower = 1, upper = num_obs * num_deworming_days> num_relevant_obs_days = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any]));
  
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  
  matrix[num_relevant_obs_days, num_census_covar_coef] census_covar_dm_long;
  matrix[num_relevant_obs_days, num_all_treatment_coef] treatment_design_matrix_long;
  int cluster_id_long[num_relevant_obs_days];
  
  int<lower = 0, upper = 1> relevant_dewormed_any_daily[num_relevant_obs_days] = rep_array(0, num_relevant_obs_days);
  int<lower = 1, upper = num_deworming_days> dewormed_day_long[num_relevant_obs_days];
  int<lower = 0> relevant_daily_strata_sizes[num_strata];
  
  {
    int stratum_pos = 1;
    int relevant_day_obs_pos = 1;
   
    for (stratum_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[stratum_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;

      relevant_daily_strata_sizes[stratum_index] = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any[stratum_pos:stratum_end]]));

      stratum_pos = stratum_end + 1;
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
  real hyper_intercept_raw;
  row_vector[num_all_treatment_coef - 1] hyper_treat_beta_day1; 
  matrix[num_deworming_days - 1, num_all_treatment_coef] hyper_beta_dyn_effect; 
  
  vector<lower = 0>[num_all_treatment_coef] strata_beta_day1_tau_raw;
  vector<lower = 0>[num_deworming_days - 1] strata_beta_dyn_tau_raw[num_all_treatment_coef];
 
  cholesky_factor_corr[num_all_treatment_coef] strata_beta_day1_L_corr_mat;
  cholesky_factor_corr[num_deworming_days - 1] strata_beta_dyn_L_corr_mat[num_all_treatment_coef];
  
  matrix[num_strata, num_all_treatment_coef] strata_beta_day1; 
  matrix[num_deworming_days - 1, num_all_treatment_coef] strata_beta_dyn_effect[num_strata]; 
}

transformed parameters {
  real hyper_intercept = hyper_intercept_raw * hyper_intercept_sigma;
  row_vector[num_all_treatment_coef] hyper_beta_day1 = append_col(hyper_intercept, hyper_treat_beta_day1); 
  matrix[num_deworming_days, num_all_treatment_coef] hyper_beta = rep_matrix(hyper_beta_day1, num_deworming_days); 
  
  vector<lower = 0>[num_all_treatment_coef] strata_beta_day1_tau = strata_beta_day1_tau_raw * scale_sigma;
  vector<lower = 0>[num_deworming_days - 1] strata_beta_dyn_tau[num_all_treatment_coef];
  
  matrix[num_deworming_days, num_all_treatment_coef] strata_beta[num_strata];
  
  hyper_beta[2:num_deworming_days] = hyper_beta[2:num_deworming_days] + hyper_beta_dyn_effect;
   
  for (beta_index in 1:num_all_treatment_coef) {
    strata_beta_dyn_tau[beta_index] = strata_beta_dyn_tau_raw[beta_index] * scale_sigma;
  }
 
  for (stratum_index in 1:num_strata) {
    strata_beta[stratum_index] = rep_matrix(strata_beta_day1[stratum_index], num_deworming_days);
    strata_beta[stratum_index, 2:num_deworming_days] = strata_beta[stratum_index, 2:num_deworming_days] + strata_beta_dyn_effect[stratum_index];
  } 
}

model {
  hyper_intercept_raw ~ normal(0, 1);
  
  hyper_treat_beta_day1 ~ normal(0, hyper_coef_sigma);
  
  to_vector(hyper_beta_dyn_effect) ~ normal(0, hyper_coef_sigma);
  
  strata_beta_day1_tau_raw ~ normal(0, 1);
  strata_beta_dyn_tau_raw ~ multi_normal(mu_dyn_zero, dyn_identity_mat);

  strata_beta_day1_L_corr_mat ~ lkj_corr_cholesky(lkj_df);
  
  {
    matrix[num_all_treatment_coef, num_all_treatment_coef] strata_beta_day1_L_vcov = diag_pre_multiply(strata_beta_day1_tau, strata_beta_day1_L_corr_mat);
    matrix[num_deworming_days - 1, num_deworming_days - 1] strata_beta_dyn_L_vcov[num_all_treatment_coef];
    
    vector[num_relevant_obs_days] latent_var = rep_vector(0, num_relevant_obs_days);
    
    int relevant_daily_stratum_pos = 1;
    
    for (beta_index in 1:num_all_treatment_coef) {
      strata_beta_dyn_L_corr_mat[beta_index] ~ lkj_corr_cholesky(lkj_df);

      strata_beta_dyn_L_vcov[beta_index] = diag_pre_multiply(strata_beta_dyn_tau[beta_index], strata_beta_dyn_L_corr_mat[beta_index]);
    }
    
    for (stratum_index in 1:num_strata) {
      int curr_relevant_daily_stratum_size = relevant_daily_strata_sizes[stratum_index];
      int relevant_daily_stratum_end = relevant_daily_stratum_pos + curr_relevant_daily_stratum_size - 1;
      
      strata_beta_day1[stratum_index] ~ multi_normal_cholesky(hyper_beta_day1, strata_beta_day1_L_vcov);
      
      for (beta_index in 1:num_all_treatment_coef) {
        strata_beta_dyn_effect[stratum_index, , beta_index] ~ multi_normal_cholesky(hyper_beta_dyn_effect[, beta_index], 
                                                                                    strata_beta_dyn_L_vcov[beta_index]);
        // strata_beta_dyn_effect[stratum_index, , beta_index] ~ multi_normal_cholesky(mu_dyn_zero, strata_beta_dyn_L_vcov[beta_index]);
      }
      
      latent_var[relevant_daily_stratum_pos:relevant_daily_stratum_end] =
        rows_dot_product(treatment_design_matrix_long[relevant_daily_stratum_pos:relevant_daily_stratum_end],
                         strata_beta[stratum_index, dewormed_day_long[relevant_daily_stratum_pos:relevant_daily_stratum_end]]);
      
      relevant_daily_stratum_pos = relevant_daily_stratum_end + 1;
    }
   
    /* //  Non-multilevel version 
    vector[num_relevant_obs_days] latent_var = rows_dot_product(treatment_design_matrix_long, hyper_beta[dewormed_day_long]);
    */
    
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
    
    matrix[num_strata, num_ate_treatments] latent_var_map_day1 = strata_beta_day1 * (treatment_map_design_matrix[ate_treatments[, 1]]');
    // matrix[num_ate_treatments, num_strata] latent_var_map_day1 = treatment_map_design_matrix[ate_treatments[, 1]] * (strata_beta_day1');
    matrix[num_strata, num_deworming_days] latent_var_map[num_ate_treatments];
    // matrix[num_deworming_days, num_strata] latent_var_map[num_ate_treatments];
    
    for (ate_treat_index in 1:num_ate_treatments) {
      latent_var_map[ate_treat_index] = rep_matrix(latent_var_map_day1[, ate_treat_index], num_deworming_days);
      
      for (stratum_index in 1:num_strata) {
        latent_var_map[ate_treat_index, stratum_index, 2:num_deworming_days] = latent_var_map[ate_treat_index, stratum_index, 2:num_deworming_days] 
          + treatment_map_design_matrix[ate_treatments[ate_treat_index, 2]] * (strata_beta_dyn_effect[stratum_index]');
        // latent_var_map[ate_treat_index, 2:num_deworming_days, stratum_index] = latent_var_map[ate_treat_index, 2:num_deworming_days, stratum_index] 
        //   + treatment_map_design_matrix[ate_treatments[, 2]] * (strata_beta_dyn_effect[stratum_index]');
      }
    }

    // vector[num_obs] census_covar_latent_var; 
    // matrix[num_strata, num_deworming_days] latent_var_map_dyn[num_ate_treatments]; 
    // matrix[num_strata, num_ate_treatments] latent_var_map_day1 = treatment_map_design_matrix[ate_treatments[, 1]] * (strata_beta_day1');

    for (stratum_index in 1:num_strata) {
      // int curr_stratum_size = strata_sizes[stratum_index];
      // int stratum_end = stratum_pos + curr_stratum_size - 1;

      // census_covar_latent_var[stratum_pos:stratum_end] = census_covar_dm[stratum_pos:stratum_end] * stratum_census_covar_coef_mat[, stratum_index];

      // stratum_pos = stratum_end + 1;
    }

    est_deworming_days =
      treatment_cell_deworming_day_prop_rng(ate_treatments,
                                            missing_obs_ids,
                                            missing_treatment_stratum_id,
                                            missing_treatment_cluster_id,
                                            missing_treatment_sizes,
                                            observed_treatment_sizes,
                                            // cluster_effects,
                                            // rep_matrix(census_covar_latent_var, num_deworming_days),
                                            latent_var_map,
                                            // hyper_dyn_latent_var,
                                            observed_dewormed_day,
                                            use_logit);

    // est_takeup = 1 - est_deworming_days[, 13];
    // est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }
}
