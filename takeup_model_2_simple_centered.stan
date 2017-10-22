functions {
  matrix calculate_stratum_log_lambda_t(vector day_constant_kappa, matrix dyn_treat_dm, vector dyn_treat_coef, row_vector stratum_hazard, vector cluster_frailty) {
    int num_obs = rows(day_constant_kappa);
    int num_days = cols(stratum_hazard);
    
    return rep_matrix(day_constant_kappa, num_days)  
      + to_matrix(dyn_treat_dm * dyn_treat_coef, num_obs, num_days, 0)  
      + log(rep_matrix(stratum_hazard, num_obs) .* rep_matrix(cluster_frailty, num_days));
  } 
  
  matrix calculate_treatment_log_lambda_t(vector day_constant_kappa, matrix dyn_treat_dm, matrix dyn_treat_coef, matrix stratum_hazard, vector cluster_frailty) {
    int num_obs = rows(day_constant_kappa);
    int num_days = cols(stratum_hazard);
    
    return rep_matrix(day_constant_kappa, num_days) + (dyn_treat_coef * dyn_treat_dm') + log(stratum_hazard .* rep_matrix(cluster_frailty, num_days));
  } 
 
  matrix treatment_deworming_day_rng(vector log_kappa_census_covar, 
                                     matrix treatment_coef, matrix dyn_treatment_coef, matrix hazard, vector cluster_frailty,
                                     vector treatment_dm, matrix dyn_treatment_dm) {
                                     // int[] private_value_calendar_coef, int[] private_value_bracelet_coef) {
    int curr_treatment_size = num_elements(log_kappa_census_covar); // rows(census_covar_coef); 
    int num_deworming_days = cols(hazard);
    
    matrix[curr_treatment_size, num_deworming_days + 1] deworming_days_mask; //[num_treatments]; 
    
    // vector[curr_treatment_size] log_kappa = treatment_coef * treatment_dm + log_kappa_census_covar; 
      // treatment_coef[, private_value_calendar_coef] * treatment_dm[private_value_bracelet_coef]; 
      
    matrix[curr_treatment_size, num_deworming_days] log_lambda_t = 
      calculate_treatment_log_lambda_t(treatment_coef * treatment_dm + log_kappa_census_covar, dyn_treatment_dm, dyn_treatment_coef, hazard, cluster_frailty);
    
    for (obs_ids_index in 1:curr_treatment_size) {
      vector[num_deworming_days + 1] deworming_day_prob = rep_vector(1, num_deworming_days + 1);
      // vector[num_deworming_days] one_m_alphas = 1 - gumbel_cdf(log_lambda_t[obs_ids_index], 0, 1)'; //inv_cloglog(log_lambda_t[obs_ids_index])';
      vector[num_deworming_days + 1] cumul_prod_alphas = rep_vector(1, num_deworming_days + 1);
    
      real prev_not_deworm_prob = 1; 
      // deworming_day_prob[1:num_deworming_days] = one_m_alphas;
      
      // for (day_index in 2:(num_deworming_days + 1)) {
      for (day_index in 1:num_deworming_days) {
        cumul_prod_alphas[day_index] = prev_not_deworm_prob; // prod(1 - one_m_alphas[1:(day_index - 1)])';
        deworming_day_prob[day_index] = 1 - gumbel_cdf(log_lambda_t[obs_ids_index, day_index], 0, 1); 
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
                                               // int[] private_value_calendar_coef,
                                               // int[] private_value_bracelet_coef,
                                               matrix missing_census_covar_dm,
                                               matrix treatment_map_dm,
                                               matrix[] dyn_treatment_map,
                                               int[] all_treat_dyn_treat_id,
                                               // int[] all_treat_signal_observed_days,
                                               // int[] all_treat_reminder_info_days,
                                               int[] missing_treatment_sizes,
                                               int[] observed_treatment_sizes,
                                               matrix stratum_hazard,
                                               vector cluster_frailty,
                                               matrix stratum_census_covar_coef,
                                               matrix stratum_treatment_coef,
                                               matrix stratum_dyn_treatment_coef,
                                               // vector observed_log_kappa_census_covar,
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
        // treatment_deworming_day_rng(observed_log_kappa_census_covar[missing_obs_ids[missing_treatment_pos:missing_treatment_end]],
        treatment_deworming_day_rng(rows_dot_product(stratum_census_covar_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                                     missing_census_covar_dm[missing_treatment_pos:missing_treatment_end]),
                                    stratum_treatment_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                    stratum_dyn_treatment_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                    stratum_hazard[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                    cluster_frailty[missing_cluster_id[missing_treatment_pos:missing_treatment_end]],
        
                                    treatment_map_dm[treatment_ids[treatment_ids_index], ]',
                                    dyn_treatment_map[curr_all_treatment_id]); 
                                                              // all_treat_signal_observed_days[curr_all_treatment_id] + 1, all_treat_reminder_info_days[curr_all_treatment_id] + 1],
                                    // private_value_calendar_coef, private_value_bracelet_coef);
      
      matrix[curr_observed_treatment_size, num_deworming_days + 1] observed_deworming_days_mask = days_diag[observed_dewormed_any[observed_treatment_pos:observed_treatment_end]];
      
      cell_deworming_day[treatment_ids_index, 1:(num_deworming_days + 1)] = 
        (diagonal(crossprod(missing_deworming_days_mask)) + diagonal(crossprod(observed_deworming_days_mask)))' / (curr_missing_treatment_size + curr_observed_treatment_size);
        
      // cell_deworming_day[treatment_ids_index, num_deworming_days + 2] = curr_missing_treatment_size + curr_observed_treatment_size;
        
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
  
  // int num_private_value_calendar_coef;
  // int num_private_value_bracelet_coef;
  // int private_value_calendar_coef[num_private_value_calendar_coef];
  // int private_value_bracelet_coef[num_private_value_bracelet_coef];
  // int not_private_value_bracelet_coef[num_all_treatment_coef - num_private_value_bracelet_coef];
  
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
  
  matrix[num_deworming_days, num_deworming_days] hazard_day_map;
  matrix[num_deworming_days + 1, num_deworming_days] hazard_day_triangle_map;
  
  // Covariates
  
  int<lower = 1> num_census_covar_coef;
  matrix[num_obs, num_census_covar_coef] census_covar_dm;
  
  // Dynamics
  
  int<lower = 0> num_dynamic_treatments;
  
  int<lower = 1, upper = num_obs> dynamic_treatment_id[num_obs];
  int<lower = 0> num_dynamic_treatment_col;
 
  matrix[num_deworming_days, num_dynamic_treatment_col] dynamic_treatment_map[num_dynamic_treatments];
  int all_treatment_dyn_id[num_all_treatments];
  
  // int<lower = 1, upper = num_dynamic_treatment_col> signal_observed_coef[num_signal_observed_coef];
  // int<lower = 1, upper = num_dynamic_treatment_col> reminder_info_coef[num_reminder_info_coef];
  
  // int<lower = 0, upper = num_deworming_days> all_treatment_signal_observed_days[num_all_treatments];
  // int<lower = 0, upper = num_deworming_days> all_treatment_reminder_info_days[num_all_treatments];
  // 
  // matrix<lower = 0>[num_deworming_days, num_dynamic_treatment_col] dynamic_treatment_map_dm[num_dynamic_treatments, num_deworming_days, num_deworming_days];
 
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
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  
  // int num_not_private_value_bracelet_coef = num_all_treatment_coef - num_private_value_bracelet_coef;
  
  matrix<lower = 0>[num_obs * num_deworming_days, num_dynamic_treatment_col] dynamic_treatment_dm;
  
  {
    int dynamic_treatment_dm_pos = 1;
    
    for (obs_index in 1:num_obs) {
      dynamic_treatment_dm[dynamic_treatment_dm_pos:(dynamic_treatment_dm_pos + num_deworming_days - 1)] = dynamic_treatment_map[dynamic_treatment_id[obs_index]];
      
      dynamic_treatment_dm_pos = dynamic_treatment_dm_pos + num_deworming_days;
    }
  }
}

parameters {
  row_vector<lower = 0, upper = 1>[num_deworming_days] hyper_baseline_cond_takeup; // Uniform[0, 1] prior
  
  real<lower = 0> stratum_hazard_effect[num_strata];
  real<lower = 0> stratum_hazard_frailty_var;
  
  vector<lower = 0>[num_clusters] cluster_hazard_effect;
  real<lower = 0> cluster_hazard_frailty_var;
  
  // vector[num_not_private_value_bracelet_coef] hyper_beta; 
  vector[num_all_treatment_coef] hyper_beta; 
  // vector[num_not_private_value_bracelet_coef] stratum_beta[num_strata];
  vector[num_all_treatment_coef] stratum_beta[num_strata];
  // vector<lower = 0>[num_not_private_value_bracelet_coef] stratum_tau_treatment;
  vector<lower = 0>[num_all_treatment_coef] stratum_tau_treatment;
   
  row_vector[num_dynamic_treatment_col] hyper_dynamic_treatment_coef;
  row_vector[num_dynamic_treatment_col] stratum_dynamic_treatment_coef[num_strata];
  vector<lower = 0>[num_dynamic_treatment_col] stratum_tau_dynamic_treatment;
  
  vector[num_census_covar_coef] hyper_census_covar_coef;
  vector[num_census_covar_coef] stratum_census_covar_coef[num_strata];
  vector<lower = 0>[num_census_covar_coef] stratum_tau_census_covar;
}

transformed parameters {
  row_vector<lower = 0>[num_deworming_days] hyper_baseline_hazard = - log(1 - hyper_baseline_cond_takeup);
  
  // vector[num_strata] stratum_lp;
  
  matrix[num_strata, num_all_treatment_coef] stratum_beta_mat;
  matrix[num_strata, num_dynamic_treatment_col] stratum_dyn_treatment_mat;
  matrix[num_strata, num_census_covar_coef] stratum_census_covar_coef_mat;
  matrix[num_strata, num_deworming_days] stratum_hazard_mat;
 
  for (strata_index in 1:num_strata) {
    stratum_beta_mat[strata_index] = stratum_beta[strata_index]';
    stratum_dyn_treatment_mat[strata_index] = stratum_dynamic_treatment_coef[strata_index];
    stratum_census_covar_coef_mat[strata_index] = stratum_census_covar_coef[strata_index]';
    stratum_hazard_mat[strata_index] = stratum_hazard_effect[strata_index] * hyper_baseline_hazard;
  } 
}

model {
  cluster_hazard_frailty_var ~ student_t(scale_df, 0, scale_sigma);
  cluster_hazard_effect ~ gamma(1 / cluster_hazard_frailty_var, 1 / cluster_hazard_frailty_var);
  
  stratum_hazard_frailty_var ~ student_t(scale_df, 0, scale_sigma);
  stratum_hazard_effect ~ gamma(1 / stratum_hazard_frailty_var, 1 / stratum_hazard_frailty_var); 
  
  hyper_beta ~ student_t(coef_df, 0, coef_sigma); 
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  stratum_beta ~ multi_student_t(coef_df, hyper_beta, diag_matrix(stratum_tau_treatment));
  
  hyper_dynamic_treatment_coef ~ student_t(coef_df, 0, coef_sigma_dynamic);
  stratum_tau_dynamic_treatment ~ student_t(scale_df, 0, scale_sigma_dynamic);
  stratum_dynamic_treatment_coef ~ multi_student_t(coef_df, hyper_dynamic_treatment_coef, diag_matrix(stratum_tau_dynamic_treatment));
  
  hyper_census_covar_coef ~ student_t(coef_df, 0, coef_sigma);
  stratum_tau_census_covar ~ student_t(scale_df, 0, scale_sigma);
  stratum_census_covar_coef ~ multi_student_t(coef_df, hyper_census_covar_coef, diag_matrix(stratum_tau_census_covar));
  
  {
    int stratum_pos = 1;
    int dewormed_stratum_pos = 1;
    int dynamic_stratum_pos = 1;
  
    // vector[num_obs] observed_log_kappa = rep_vector(0, num_obs);
    // matrix[num_obs, num_deworming_days] observed_log_lambda_t = rep_matrix(0, num_obs, num_deworming_days);
    // vector[num_obs] observed_log_kappa_census_covar = rep_vector(0, num_obs);
    
    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      int curr_dewormed_stratum_size = strata_dewormed_sizes[strata_index];
      int dewormed_stratum_end = dewormed_stratum_pos + curr_dewormed_stratum_size - 1;
      
      int curr_dynamic_stratum_size = curr_stratum_size * num_deworming_days;  
      int dynamic_stratum_end = dynamic_stratum_pos + curr_dynamic_stratum_size - 1;
      
      // vector[num_all_treatment_coef] local_stratum_beta = rep_vector(0, num_all_treatment_coef);
        
      int stratum_dewormed_ids[curr_dewormed_stratum_size] = dewormed_ids[dewormed_stratum_pos:dewormed_stratum_end];
      int stratum_dewormed_day_all[curr_stratum_size] = dewormed_day_any[stratum_pos:stratum_end]; 
      int stratum_dewormed_day_dewormed[curr_dewormed_stratum_size] = dewormed_day_any[stratum_dewormed_ids];
      
      matrix[curr_stratum_size, num_deworming_days] stratum_hazard_day_triangle_map = hazard_day_triangle_map[stratum_dewormed_day_all];
      matrix[curr_dewormed_stratum_size, num_deworming_days] stratum_hazard_day_map = hazard_day_map[stratum_dewormed_day_dewormed];
      
      // row_vector[num_deworming_days] stratum_baseline_hazard = stratum_hazard_effect[strata_index] * hyper_baseline_hazard;
      
      // vector[curr_stratum_size] stratum_log_kappa; // = rep_vector(0, num_obs);
      // matrix[curr_stratum_size, num_deworming_days] stratum_log_lambda_t; // = rep_matrix(0, num_obs, num_deworming_days);
      
      // local_stratum_beta[not_private_value_bracelet_coef] = stratum_beta[strata_index];
      // local_stratum_beta[private_value_bracelet_coef] = local_stratum_beta[private_value_calendar_coef];
      
      // observed_log_kappa_census_covar[stratum_pos:stratum_end] = census_covar_dm[stratum_pos:stratum_end] * stratum_census_covar_coef[strata_index];
              
      // stratum_log_kappa =
      //   census_covar_dm[stratum_pos:stratum_end] * stratum_census_covar_coef[strata_index] +
      //     // observed_log_kappa_census_covar[stratum_pos:stratum_end] +
      //     treatment_design_matrix[stratum_pos:stratum_end] * local_stratum_beta +
      //     treatment_design_matrix[stratum_pos:stratum_end, private_value_bracelet_coef] * local_stratum_beta[private_value_calendar_coef];
      
      matrix[curr_stratum_size, num_deworming_days] stratum_log_lambda_t =
        calculate_stratum_log_lambda_t(census_covar_dm[stratum_pos:stratum_end] * stratum_census_covar_coef[strata_index] +
                                         treatment_design_matrix[stratum_pos:stratum_end] * stratum_beta[strata_index],
                                       dynamic_treatment_dm[dynamic_stratum_pos:dynamic_stratum_end],
                                       stratum_dynamic_treatment_coef[strata_index]',
                                       stratum_hazard_mat[strata_index],
                                       cluster_hazard_effect[cluster_id[stratum_pos:stratum_end]]);
       //  // This is the time-invariant part of the kappa function
       //  rep_matrix(census_covar_dm[stratum_pos:stratum_end] * stratum_census_covar_coef[strata_index]  
       //              + treatment_design_matrix[stratum_pos:stratum_end] * local_stratum_beta, num_deworming_days) 
       //  
       //  // This is the time varying part of the kappa function              
       //  + to_matrix(dynamic_treatment_dm[dynamic_stratum_pos:dynamic_stratum_end] * stratum_dynamic_treatment_coef[strata_index]', curr_stratum_size, num_deworming_days, 0)  
       // 
       // // Hazard and cluster frailty (basically cluster effects) 
       //  + log(rep_matrix(stratum_hazard_effect[strata_index] * hyper_baseline_hazard, curr_stratum_size) .* 
       //          rep_matrix(cluster_hazard_effect[cluster_id[stratum_pos:stratum_end]], num_deworming_days));
          
        
      // stratum_lp[strata_index] = 
      target += 
        gumbel_lcdf(to_vector(stratum_log_lambda_t .* stratum_hazard_day_triangle_map) | 0, 1) + 
        gumbel_lccdf(to_vector(stratum_log_lambda_t[stratum_dewormed_index[dewormed_stratum_pos:dewormed_stratum_end]] .* stratum_hazard_day_map) | 0, 1); 
        // - (exp(log_sum_exp(stratum_log_lambda_t .* stratum_hazard_day_triangle_map)) - sum(1 - stratum_hazard_day_triangle_map)) 
        // + sum(log(inv_cloglog(stratum_log_lambda_t[stratum_dewormed_index[dewormed_stratum_pos:dewormed_stratum_end]] .* stratum_hazard_day_map))) 
        // - curr_dewormed_stratum_size * (num_deworming_days - 1) * log1m_exp(-1);
        
      // if (is_nan(stratum_lp[strata_index]) || is_inf(stratum_lp[strata_index])) {
      //   reject("Stratum ", strata_index, ": log probability is ", stratum_lp[strata_index]);
      // }
        
          
      stratum_pos = stratum_end + 1;
      dewormed_stratum_pos = dewormed_stratum_end + 1;
      dynamic_stratum_pos = dynamic_stratum_end + 1;
    }
  }
  
  // target += stratum_lp;
}

generated quantities {
  row_vector<lower = 0, upper = 1>[num_deworming_days] stratum_baseline_cond_takeup[num_strata];
  
  matrix<lower = 0>[num_ate_treatments, num_deworming_days + 1] est_deworming_days = rep_matrix(0, num_ate_treatments, num_deworming_days + 1);
  vector<lower = 0, upper = 1>[num_ate_treatments] est_takeup = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] est_takeup_ate = rep_vector(0, num_ate_pairs);
  
  for (stratum_index in 1:num_strata) {
    stratum_baseline_cond_takeup[stratum_index] = exp(- hyper_baseline_hazard * stratum_hazard_effect[stratum_index]);
  }
  
  if (estimate_ate) {
    est_deworming_days =
      treatment_cell_deworming_day_prop_rng(ate_treatments,
                                            missing_obs_ids,
                                            missing_treatment_stratum_id,
                                            missing_treatment_cluster_id,
                                       //   private_value_calendar_coef,
                                       // private_value_bracelet_coef,
                                            missing_census_covar_dm,
                                            treatment_map_design_matrix,
                                            dynamic_treatment_map,
                                            all_treatment_dyn_id,
                                       // all_treatment_signal_observed_days,
                                       // all_treatment_reminder_info_days,
                                            missing_treatment_sizes,
                                            observed_treatment_sizes,
                                            stratum_hazard_mat,
                                            cluster_hazard_effect,
                                            stratum_census_covar_coef_mat,
                                            stratum_beta_mat,
                                            stratum_dyn_treatment_mat,
                                       // observed_log_kappa_census_covar,
                                            observed_dewormed_day);
   
    est_takeup = 1 - est_deworming_days[, 13]; 
    est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }
    
}
