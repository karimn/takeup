functions {
  matrix cov_exp_quad_ARD(row_vector[] x, real alpha, row_vector rho, real delta, int decompose) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    
    for (i in 1:(N - 1)) {
      K[i,i] = sq_alpha + delta;
      
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp((- 0.5) * dot_self((x[i] - x[j]) ./ rho));
        K[j, i] = K[i, j];
      }
    }
    
    K[N,N] = sq_alpha + delta;
   
    if (decompose) { 
      return cholesky_decompose(K);
    } else {
      return K;
    }
  }
  
  matrix treatment_deworming_day_rng(matrix census_covar_latent_var, 
                                     matrix day_constant_treat_latent_var,
                                     row_vector dynamic_latent_var, 
                                     vector cluster_effects) {
    int curr_treatment_size = rows(census_covar_latent_var); 
    int num_deworming_days = cols(census_covar_latent_var);
    
    matrix[curr_treatment_size, num_deworming_days + 1] deworming_days_mask = rep_matrix(0, curr_treatment_size, num_deworming_days + 1); 
      
    matrix[curr_treatment_size, num_deworming_days] log_lambda_t = 
      census_covar_latent_var +
      day_constant_treat_latent_var +
      rep_matrix(dynamic_latent_var, curr_treatment_size) + rep_matrix(cluster_effects, num_deworming_days);
    
    for (obs_ids_index in 1:curr_treatment_size) {
      vector[num_deworming_days + 1] deworming_day_prob = rep_vector(1, num_deworming_days + 1);
      vector[num_deworming_days + 1] cumul_prod_alphas = rep_vector(1, num_deworming_days + 1);
    
      real prev_not_deworm_prob = 1; 
      
      for (day_index in 1:num_deworming_days) {
        cumul_prod_alphas[day_index] = prev_not_deworm_prob; 
        deworming_day_prob[day_index] = 1 - gumbel_cdf(- log_lambda_t[obs_ids_index, day_index], 0, 1); 
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
  matrix treatment_cell_deworming_day_prop_rng(int[] treatment_ids, 
                                               int[] missing_obs_ids, 
                                               int[] missing_stratum_id,
                                               int[] missing_cluster_id,
                                               int[] all_treat_dyn_treat_id,
                                               int[] missing_treatment_sizes,
                                               int[] observed_treatment_sizes,
                                               vector cluster_effects,
                                               matrix census_covar_latent_var,
                                               matrix day_constant_treat_latent_var,
                                               matrix dynamic_latent_var,
                                               int[] observed_dewormed_any) {
    int num_treatment_ids = size(treatment_ids);
    int num_deworming_days = cols(census_covar_latent_var);
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    matrix[num_treatment_ids, num_deworming_days + 1] cell_deworming_day = rep_matrix(0, num_treatment_ids, num_deworming_days + 1);
    
    matrix[num_deworming_days + 1, num_deworming_days + 1] days_diag = diag_matrix(rep_vector(1, num_deworming_days + 1));
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      int curr_all_treatment_id = treatment_ids[treatment_ids_index];
      int curr_dyn_treatment_id = all_treat_dyn_treat_id[treatment_ids_index];
      
      matrix[curr_missing_treatment_size, num_deworming_days + 1] missing_deworming_days_mask = 
        treatment_deworming_day_rng(census_covar_latent_var[missing_obs_ids[missing_treatment_pos:missing_treatment_end]],
                                    rep_matrix(day_constant_treat_latent_var[treatment_ids_index, missing_stratum_id[missing_treatment_pos:missing_treatment_end]]', num_deworming_days),
                                    dynamic_latent_var[curr_dyn_treatment_id],
                                    cluster_effects[missing_cluster_id[missing_treatment_pos:missing_treatment_end]]);
                                    
      if (curr_observed_treatment_size > 0) {
        int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
        
        matrix[curr_observed_treatment_size, num_deworming_days + 1] observed_deworming_days_mask = days_diag[observed_dewormed_any[observed_treatment_pos:observed_treatment_end]];
      
        cell_deworming_day[treatment_ids_index, 1:(num_deworming_days + 1)] = 
          (diagonal(crossprod(missing_deworming_days_mask)) + diagonal(crossprod(observed_deworming_days_mask)))' / (curr_missing_treatment_size + curr_observed_treatment_size);
          
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
  
  int<lower = 0> num_dynamic_treatments;
  
  int<lower = 1, upper = num_obs> dynamic_treatment_id[num_obs];
  int<lower = 0> num_dynamic_treatment_col;
 
  matrix[num_deworming_days, num_dynamic_treatment_col] dynamic_treatment_map[num_dynamic_treatments];
  int all_treatment_dyn_id[num_all_treatments];
  int relevant_latent_var_map[num_deworming_days + 1, num_deworming_days];
 
  // Counterfactuals and ATE
  
  int<lower = 0, upper = num_all_treatments> num_ate_treatments;
  int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments];
  
  int<lower = 0> num_missing_obs_ids;
  int<lower = 0, upper = num_obs> missing_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> missing_obs_ids[num_missing_obs_ids];
  
  int<lower = 0> num_observed_obs_ids;
  int<lower = 0, upper = num_obs> observed_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> observed_obs_ids[num_observed_obs_ids];
  
  int<lower = 0, upper = num_all_treatments> num_ate_pairs;
  
  int<lower = 1, upper = num_all_treatments> ate_pairs[num_ate_pairs, 2];
  
  int missing_treatment_stratum_id[num_missing_obs_ids]; 
  int missing_treatment_cluster_id[num_missing_obs_ids]; 
  
  matrix[num_missing_obs_ids, num_census_covar_coef] missing_census_covar_dm;
  
  int observed_dewormed_day[num_observed_obs_ids];
  
  // Constants for hyperpriors 
  
  real<lower = 0> scale_df;
  real<lower = 0> scale_sigma;
  
  real<lower = 0> coef_df;
  real<lower = 0> hyper_coef_sigma;
  
  real<lower = 0> lkj_df;
  
  real<lower = 0> hyper_dyn_rho_shape; // alpha
  real<lower = 0> hyper_dyn_rho_scale; // beta
  
  // Configuration
  
  int<lower = 0, upper = 1> estimate_ate;
  int<lower = 0, upper = 1> use_logit; 
  
  int<lower = 1, upper = num_deworming_days> first_dynamics_day;
}

transformed data {
  real delta = 1e-9; 
  
  int<lower = 1, upper = num_obs * num_deworming_days> num_relevant_obs_days = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any]));
  int<lower = 1, upper = num_relevant_obs_days> num_relevant_dyn_obs_days = num_relevant_obs_days - (num_obs * (first_dynamics_day - 1));
  
  int<lower = 1, upper = num_deworming_days> num_dynamics_days = num_deworming_days - first_dynamics_day + 1;
  
  int<lower = 0, upper = num_dynamic_treatments * num_deworming_days> num_dynamic_treatment_days = num_dynamic_treatments * num_dynamics_days;
  
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  row_vector[num_dynamic_treatment_col] nonparam_dyn_treatment_map[num_dynamic_treatment_days];
  
  matrix[num_relevant_obs_days, num_census_covar_coef] census_covar_dm_long;
  matrix[num_relevant_obs_days, num_all_treatment_coef] treatment_design_matrix_long;
  int cluster_id_long[num_relevant_obs_days];
  
  int<lower = 0, upper = 1> relevant_dewormed_any_daily[num_relevant_obs_days] = rep_array(0, num_relevant_obs_days);
  int<lower = 1, upper = num_deworming_days> dewormed_day_long[num_relevant_obs_days];
  int<lower = 1, upper = num_obs * num_deworming_days> relevant_daily_ids[num_relevant_obs_days];
  int<lower = 1, upper = num_relevant_obs_days> relevant_dyn_daily_ids[num_relevant_dyn_obs_days];
  int<lower = 0> relevant_daily_strata_sizes[num_strata];
  
  {
    int dynamic_treatment_dm_pos = 1;
    int non_param_dyn_treat_pos = 1;
    int stratum_pos = 1;
    int relevant_day_obs_pos = 1;
    int relevant_dyn_day_obs_pos = 1;
    
    // vector[num_deworming_days] deworming_day_map = rep_vector(1, num_deworming_days);
    // 
    // deworming_day_map = cumulative_sum(deworming_day_map);
    
    for (dyn_treat_index in 1:num_dynamic_treatments) {
      for (dyn_day_index in 1:num_dynamics_days) {
        nonparam_dyn_treatment_map[non_param_dyn_treat_pos + dyn_day_index - 1] = dynamic_treatment_map[dyn_treat_index, dyn_day_index + first_dynamics_day - 1];
      }
     
      non_param_dyn_treat_pos = non_param_dyn_treat_pos + num_dynamics_days;
    }
   
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
      int relevant_dyn_day_obs_end = relevant_dyn_day_obs_pos + curr_relevant_days - first_dynamics_day;
      
      census_covar_dm_long[relevant_day_obs_pos:relevant_day_obs_end] = rep_matrix(census_covar_dm[obs_index], curr_relevant_days);
      
      treatment_design_matrix_long[relevant_day_obs_pos:relevant_day_obs_end] = rep_matrix(treatment_map_design_matrix[obs_treatment[obs_index]], curr_relevant_days);
        
      cluster_id_long[relevant_day_obs_pos:relevant_day_obs_end] = rep_array(cluster_id[obs_index], curr_relevant_days);
      
      
      for (obs_rel_day_index in 1:curr_relevant_days) {
        int curr_obs_day_pos = relevant_day_obs_pos + obs_rel_day_index - 1;
        
        relevant_daily_ids[curr_obs_day_pos] = ((obs_index - 1) * num_deworming_days) + obs_rel_day_index;
        dewormed_day_long[curr_obs_day_pos] = obs_rel_day_index;
        
        if (obs_rel_day_index >= first_dynamics_day) {
          int curr_dyn_obs_day_pos = relevant_dyn_day_obs_pos + obs_rel_day_index - first_dynamics_day;
          
          relevant_dyn_daily_ids[curr_dyn_obs_day_pos] = curr_obs_day_pos;
        }
      }
     
      if (dewormed_any[obs_index]) {
        relevant_dewormed_any_daily[relevant_day_obs_pos + dewormed_day_any[obs_index] - 1] = 1;
      }
      
      relevant_day_obs_pos = relevant_day_obs_end + 1;
      relevant_dyn_day_obs_pos = relevant_dyn_day_obs_end + 1;
    }
  }
}

parameters {
  real<lower = 0, upper = 1> baseline_takeup;
  row_vector[num_all_treatment_coef - 1] hyper_treat_beta_day1; 
  matrix[num_deworming_days - 1, num_all_treatment_coef] hyper_beta_dyn; 
 
  corr_matrix[num_all_treatment_coef] strata_beta_day1_corr_mat;
  vector<lower = 0>[num_all_treatment_coef] strata_beta_day1_tau;
  corr_matrix[num_deworming_days - 1] strata_beta_dyn_corr_mat[num_all_treatment_coef];
  vector<lower = 0>[num_deworming_days - 1] strata_beta_dyn_tau;
  matrix[num_deworming_days, num_all_treatment_coef] strata_beta[num_strata]; 
}

transformed parameters {
  real hyper_intercept;
 
  row_vector[num_all_treatment_coef] hyper_beta_day1; 
  matrix[num_deworming_days, num_all_treatment_coef] hyper_beta; 
 
  if (use_logit) { 
    hyper_intercept = logit(baseline_takeup);
  } else {
    hyper_intercept = log(- log(1 - baseline_takeup));
  }
  
  hyper_beta_day1 = append_col(hyper_intercept, hyper_treat_beta_day1);
  hyper_beta = append_row(hyper_beta_day1, hyper_beta_dyn); 
}

model {
  hyper_treat_beta_day1 ~ normal(0, hyper_coef_sigma);
  
  for (deworming_day_index in 1:(num_deworming_days - 1)) {
    hyper_beta_dyn[deworming_day_index] ~ normal(hyper_beta_day1, hyper_coef_sigma); 
  }
 
  strata_beta_day1_corr_mat ~ lkj_corr(lkj_df);
  strata_beta_day1_tau ~ normal(0, scale_sigma);
  strata_beta_dyn_tau ~ normal(0, scale_sigma);
  
  {
    matrix[num_all_treatment_coef, num_all_treatment_coef] strata_beta_day1_vcov = quad_form_diag(strata_beta_day1_corr_mat, strata_beta_day1_tau);
    matrix[num_deworming_days - 1, num_deworming_days - 1] strata_beta_dyn_vcov[num_all_treatment_coef];
    
    vector[num_relevant_obs_days] latent_var = rep_vector(0, num_relevant_obs_days);
    
    int relevant_daily_stratum_pos = 1;
    
    for (beta_index in 1:num_all_treatment_coef) {
      strata_beta_dyn_corr_mat[beta_index] ~ lkj_corr(lkj_df);

      strata_beta_dyn_vcov[beta_index] = quad_form_diag(strata_beta_dyn_corr_mat[beta_index], strata_beta_dyn_tau);
    }
    
    for (stratum_index in 1:num_strata) {
      int curr_relevant_daily_stratum_size = relevant_daily_strata_sizes[stratum_index];
      int relevant_daily_stratum_end = relevant_daily_stratum_pos + curr_relevant_daily_stratum_size - 1;
      
      strata_beta[stratum_index, 1] ~ multi_normal(hyper_beta_day1, strata_beta_day1_vcov);
      // strata_beta[stratum_index, 1] ~ multi_normal(hyper_beta_day1, diag_matrix(strata_beta_day1_tau));
      
      for (beta_index in 1:num_all_treatment_coef) {
        vector[num_deworming_days - 1] strata_beta_dyn_mean = 
          rep_vector(strata_beta[stratum_index, 1, beta_index], num_deworming_days - 1) + hyper_beta_dyn[, beta_index];
        
        strata_beta[stratum_index, 2:num_deworming_days, beta_index] ~ multi_normal(strata_beta_dyn_mean, strata_beta_dyn_vcov[beta_index]);
        // strata_beta[stratum_index, 2:num_deworming_days, beta_index] ~ multi_normal(strata_beta_dyn_mean, diag_matrix(strata_beta_dyn_tau));
      }
      
      latent_var[relevant_daily_stratum_pos:relevant_daily_stratum_end] =
        rows_dot_product(treatment_design_matrix_long[relevant_daily_stratum_pos:relevant_daily_stratum_end],
                         strata_beta[stratum_index, dewormed_day_long[relevant_daily_stratum_pos:relevant_daily_stratum_end]]);
      
      relevant_daily_stratum_pos = relevant_daily_stratum_end + 1;
    }
   
    /*  Non-multilevel version 
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
}
