functions {
  matrix L_cov_exp_quad_ARD(row_vector[] x, real alpha, row_vector rho, real delta) {
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
    
    return cholesky_decompose(K);
  }
  
  matrix calculate_treatment_log_lambda_t(vector day_constant_kappa, matrix dyn_treat_dm, matrix dyn_treat_coef, matrix stratum_hazard, vector cluster_frailty) {
    int num_obs = rows(day_constant_kappa);
    int num_days = cols(stratum_hazard);
    
    return rep_matrix(day_constant_kappa, num_days) + (dyn_treat_coef * dyn_treat_dm') + log(stratum_hazard .* rep_matrix(cluster_frailty, num_days));
  } 
 
  matrix treatment_deworming_day_rng(vector log_kappa_census_covar, 
                                     matrix treatment_coef, matrix dyn_treatment_coef, matrix hazard, vector cluster_frailty,
                                     vector treatment_dm, matrix dyn_treatment_dm) {
    int curr_treatment_size = num_elements(log_kappa_census_covar); 
    int num_deworming_days = cols(hazard);
    
    matrix[curr_treatment_size, num_deworming_days + 1] deworming_days_mask; 
      
    matrix[curr_treatment_size, num_deworming_days] log_lambda_t = 
      calculate_treatment_log_lambda_t(treatment_coef * treatment_dm + log_kappa_census_covar, dyn_treatment_dm, dyn_treatment_coef, hazard, cluster_frailty);
    
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
      
      deworming_days_mask = rep_matrix(0, curr_treatment_size, num_deworming_days + 1);
      deworming_days_mask[obs_ids_index, categorical_rng(deworming_day_prob)] = 1;
    }
  
    return deworming_days_mask; 
  }
  
  /**
   * @return A (number of treatments) X (number of deworming days + 2) matrix. The last column has the total number of observed and imputed observations used to 
   * calculate the mean daily take-up
   */
  matrix treatment_cell_deworming_day_prop_rng(int[] treatment_ids, 
                                               int[] missing_obs_ids, 
                                               int[] missing_stratum_id,
                                               int[] missing_cluster_id,
                                               matrix missing_census_covar_dm,
                                               matrix treatment_map_dm,
                                               matrix[] dyn_treatment_map,
                                               int[] all_treat_dyn_treat_id,
                                               int[] missing_treatment_sizes,
                                               int[] observed_treatment_sizes,
                                               matrix stratum_hazard,
                                               vector cluster_frailty,
                                               matrix stratum_census_covar_coef,
                                               matrix stratum_treatment_coef,
                                               matrix stratum_dyn_treatment_coef,
                                               int[] observed_dewormed_any) {
    int num_treatment_ids = size(treatment_ids);
    int num_deworming_days = cols(stratum_hazard);
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    matrix[num_treatment_ids, num_deworming_days + 1] cell_deworming_day = rep_matrix(0, num_treatment_ids, num_deworming_days + 1);
    
    matrix[num_deworming_days + 1, num_deworming_days + 1] days_diag = diag_matrix(rep_vector(1, num_deworming_days + 1));
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
      int curr_all_treatment_id = all_treat_dyn_treat_id[treatment_ids_index];
      
      matrix[curr_missing_treatment_size, num_deworming_days + 1] missing_deworming_days_mask = 
        treatment_deworming_day_rng(rows_dot_product(stratum_census_covar_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                                     missing_census_covar_dm[missing_treatment_pos:missing_treatment_end]),
                                    stratum_treatment_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                    stratum_dyn_treatment_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                    stratum_hazard[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                    cluster_frailty[missing_cluster_id[missing_treatment_pos:missing_treatment_end]],
        
                                    treatment_map_dm[treatment_ids[treatment_ids_index], ]',
                                    dyn_treatment_map[curr_all_treatment_id]); 
      
      matrix[curr_observed_treatment_size, num_deworming_days + 1] observed_deworming_days_mask = days_diag[observed_dewormed_any[observed_treatment_pos:observed_treatment_end]];
      
      cell_deworming_day[treatment_ids_index, 1:(num_deworming_days + 1)] = 
        (diagonal(crossprod(missing_deworming_days_mask)) + diagonal(crossprod(observed_deworming_days_mask)))' / (curr_missing_treatment_size + curr_observed_treatment_size);
        
      missing_treatment_pos = missing_treatment_end + 1;
      observed_treatment_pos = observed_treatment_end + 1;
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
  
  int relevant_latent_var_map[num_deworming_days + 1, num_deworming_days];
  
  // Covariates
  
  int<lower = 1> num_census_covar_coef;
  matrix[num_obs, num_census_covar_coef] census_covar_dm;
  
  // Dynamics
  
  int<lower = 0> num_dynamic_treatments;
  
  int<lower = 1, upper = num_obs> dynamic_treatment_id[num_obs];
  int<lower = 0> num_dynamic_treatment_col;
 
  matrix[num_deworming_days, num_dynamic_treatment_col] dynamic_treatment_map[num_dynamic_treatments];
  int all_treatment_dyn_id[num_all_treatments];
 
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
  real<lower = 0> scale_sigma_dynamic;
  
  real<lower = 0> coef_df;
  real<lower = 0> coef_sigma;
  real<lower = 0> coef_sigma_dynamic;
  
  // Configuration
  
  int estimate_ate;
}

transformed data {
  real delta = 1e-9; 
  
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  matrix[num_deworming_days + 1, num_deworming_days] relevant_latent_var_map_mat = to_matrix(relevant_latent_var_map);
  row_vector[num_dynamic_treatment_col] nonparam_dyn_treatment_map[num_dynamic_treatments * num_deworming_days]; 
  
  matrix[num_obs * num_deworming_days, num_census_covar_coef] census_covar_dm_long;
  matrix[num_obs * num_deworming_days, num_all_treatment_coef] treatment_design_matrix_long;
  int cluster_id_long[num_obs * num_deworming_days];
  
  int<lower = 1> dewormed_days_ids[num_dewormed];
  int<lower = 1> not_dewormed_days_ids[sum(dewormed_day_any) - num_obs];
  
  {
    int dynamic_treatment_dm_pos = 1;
    int dewormed_days_ids_pos = 1;
    int not_dewormed_days_ids_pos = 1;
    int non_param_dyn_treat_pos = 1;
    
    for (dyn_treat_index in 1:num_dynamic_treatments) {
      for (dyn_day_index in 1:num_deworming_days) {
        nonparam_dyn_treatment_map[non_param_dyn_treat_pos + dyn_day_index - 1] = dynamic_treatment_map[dyn_treat_index, dyn_day_index];
      }
     
      non_param_dyn_treat_pos = non_param_dyn_treat_pos + num_deworming_days;
    }
    
    for (obs_index in 1:num_obs) {
      int not_dewormed_days = dewormed_day_any[obs_index] - 1;
      int dynamic_treatment_dm_end = dynamic_treatment_dm_pos + num_deworming_days - 1;
      
      census_covar_dm_long[dynamic_treatment_dm_pos:dynamic_treatment_dm_end] = - rep_matrix(census_covar_dm[obs_index], num_deworming_days);
      
      treatment_design_matrix_long[dynamic_treatment_dm_pos:dynamic_treatment_dm_end] = - rep_matrix(treatment_map_design_matrix[obs_treatment[obs_index]], num_deworming_days);
        
      cluster_id_long[dynamic_treatment_dm_pos:dynamic_treatment_dm_end] = rep_array(cluster_id[obs_index], num_deworming_days);
     
      if (dewormed_any[obs_index]) {
        dewormed_days_ids[dewormed_days_ids_pos] = dynamic_treatment_dm_pos + not_dewormed_days;
        dewormed_days_ids_pos = dewormed_days_ids_pos + 1;
      }
      
      if (not_dewormed_days > 0) {
        for (not_dewormed_day_index in 1:not_dewormed_days) {
          not_dewormed_days_ids[not_dewormed_days_ids_pos] = dynamic_treatment_dm_pos + not_dewormed_day_index - 1; 
          not_dewormed_days_ids_pos = not_dewormed_days_ids_pos + 1;
        }
      }
      
      dynamic_treatment_dm_pos = dynamic_treatment_dm_end + 1;
    }
  }
}

parameters {
  row_vector<lower = 0, upper = 1>[num_deworming_days] hyper_baseline_cond_takeup; // Uniform[0, 1] prior
  
  real<lower = 0> stratum_hazard_effect[num_strata];
  real<lower = 0> stratum_hazard_frailty_var;
  
  vector[num_clusters] cluster_effects;
  real<lower = 0> tau_cluster_effect;
  
  vector[num_all_treatment_coef] hyper_beta; 
  vector[num_all_treatment_coef] stratum_beta[num_strata];
  vector<lower = 0>[num_all_treatment_coef] stratum_tau_treatment;
  
  vector[num_census_covar_coef] hyper_census_covar_coef;
  vector[num_census_covar_coef] stratum_census_covar_coef[num_strata];
  vector<lower = 0>[num_census_covar_coef] stratum_tau_census_covar;
  
  //  Dynamics
  
  row_vector<lower=0>[num_dynamic_treatment_col] hyper_dyn_rho;
  real<lower=0> hyper_dyn_alpha;
  vector[num_dynamic_treatments * num_deworming_days] hyper_dyn_eta; 
}

transformed parameters {
  row_vector<lower = 0>[num_deworming_days] hyper_baseline_hazard = - log(1 - hyper_baseline_cond_takeup);
  
  matrix[num_strata, num_all_treatment_coef] stratum_beta_mat;
  matrix[num_strata, num_dynamic_treatment_col] stratum_dyn_treatment_mat;
  matrix[num_strata, num_census_covar_coef] stratum_census_covar_coef_mat;
  matrix<lower = 0>[num_strata, num_deworming_days] stratum_hazard_mat;
  
  matrix[num_dynamic_treatments, num_deworming_days] hyper_dyn_latent_var;
  
  {
    matrix[num_dynamic_treatments * num_deworming_days, num_dynamic_treatments * num_deworming_days] L_hyper_dyn_K = 
      L_cov_exp_quad_ARD(nonparam_dyn_treatment_map, hyper_dyn_alpha, hyper_dyn_rho, delta);
    
    hyper_dyn_latent_var = to_matrix(L_hyper_dyn_K * hyper_dyn_eta, num_dynamic_treatments, num_deworming_days, 0); // Row-major order
  }
 
  for (strata_index in 1:num_strata) {
    stratum_hazard_mat[strata_index] = stratum_hazard_effect[strata_index] * hyper_baseline_hazard;
    
    stratum_beta_mat[strata_index] = stratum_beta[strata_index]';
    stratum_census_covar_coef_mat[strata_index] = stratum_census_covar_coef[strata_index]';
  } 
}

model {
  stratum_hazard_frailty_var ~ student_t(scale_df, 0, scale_sigma);
  stratum_hazard_effect ~ gamma(1 / stratum_hazard_frailty_var, 1 / stratum_hazard_frailty_var);
  
  hyper_beta ~ student_t(coef_df, 0, coef_sigma); 
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  stratum_beta ~ multi_student_t(coef_df, hyper_beta, diag_matrix(stratum_tau_treatment));
  
  hyper_census_covar_coef ~ student_t(coef_df, 0, coef_sigma);
  stratum_tau_census_covar ~ student_t(scale_df, 0, scale_sigma);
  stratum_census_covar_coef ~ multi_student_t(coef_df, hyper_census_covar_coef, diag_matrix(stratum_tau_census_covar));
  
  tau_cluster_effect ~ student_t(scale_df, 0, scale_sigma);
  cluster_effects ~ student_t(coef_df, 0, tau_cluster_effect); 
 
  hyper_dyn_rho ~ gamma(4, 4); 
  hyper_dyn_alpha ~ normal(0, 1);
  hyper_dyn_eta ~ normal(0, 1);
  
  {
    vector[num_obs * num_deworming_days] latent_var;
  
    int stratum_pos = 1;
    int dynamic_stratum_pos = 1;
    
    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;

      int curr_dynamic_stratum_size = curr_stratum_size * num_deworming_days;
      int dynamic_stratum_end = dynamic_stratum_pos + curr_dynamic_stratum_size - 1;
      
      latent_var[dynamic_stratum_pos:dynamic_stratum_end] = 
        census_covar_dm_long[dynamic_stratum_pos:dynamic_stratum_end] * stratum_census_covar_coef[strata_index] + 
        treatment_design_matrix_long[dynamic_stratum_pos:dynamic_stratum_end] * stratum_beta[strata_index] -
        cluster_effects[cluster_id_long[dynamic_stratum_pos:dynamic_stratum_end]] -
        to_vector(log(rep_matrix(stratum_hazard_mat[strata_index], curr_stratum_size))) -
        to_vector(hyper_dyn_latent_var[dynamic_treatment_id[stratum_pos:stratum_end]]);
          
      stratum_pos = stratum_end + 1;
      dynamic_stratum_pos = dynamic_stratum_end + 1;
    }
    
    target += gumbel_lcdf(latent_var[not_dewormed_days_ids] | 0, 1) + gumbel_lccdf(latent_var[dewormed_days_ids] | 0, 1);
  }
}

generated quantities {
  matrix<lower = 0, upper = 1>[num_strata, num_deworming_days] stratum_baseline_cond_takeup = 1 - exp(- stratum_hazard_mat);
  
  matrix<lower = 0>[num_ate_treatments, num_deworming_days + 1] est_deworming_days = rep_matrix(0, num_ate_treatments, num_deworming_days + 1);
  vector<lower = 0, upper = 1>[num_ate_treatments] est_takeup = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] est_takeup_ate = rep_vector(0, num_ate_pairs);
  
  if (estimate_ate) {
    // est_deworming_days =
    //   treatment_cell_deworming_day_prop_rng(ate_treatments,
    //                                         missing_obs_ids,
    //                                         missing_treatment_stratum_id,
    //                                         missing_treatment_cluster_id,
    //                                         missing_census_covar_dm,
    //                                         treatment_map_design_matrix,
    //                                         dynamic_treatment_map,
    //                                         all_treatment_dyn_id,
    //                                         missing_treatment_sizes,
    //                                         observed_treatment_sizes,
    //                                         stratum_hazard_mat,
    //                                         cluster_hazard_effect,
    //                                         stratum_census_covar_coef_mat,
    //                                         stratum_beta_mat,
    //                                         stratum_dyn_treatment_mat,
    //                                         observed_dewormed_day);
    // 
    // est_takeup = 1 - est_deworming_days[, 13]; 
    // est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }
}
