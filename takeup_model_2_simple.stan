functions {
}

data {
  int<lower = 0> num_obs;
  int<lower = 1> num_all_treatments; // # of treatment cells
  int<lower = 1> num_all_treatment_coef; // # of linear model coefficients minus intercept
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  int num_private_value_calendar_coef;
  int num_private_value_bracelet_coef;
  int private_value_calendar_coef[num_private_value_calendar_coef];
  int private_value_bracelet_coef[num_private_value_bracelet_coef];
  int not_private_value_bracelet_coef[num_all_treatment_coef - num_private_value_bracelet_coef];
  
  int<lower = 0> num_bracelet_treated;
  int<lower = 1, upper = num_obs> bracelet_treated_id[num_bracelet_treated];
  int<lower = 0, upper = num_obs> strata_bracelet_sizes[num_strata];
  
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
  
  matrix[num_deworming_days, num_deworming_days] hazard_day_map;
  matrix[num_deworming_days + 1, num_deworming_days] hazard_day_triangle_map;
  
  // Dynamics
  
  int<lower = 0> num_dynamic_treatment_col;
 
  matrix[num_deworming_days, num_dynamic_treatment_col] dynamic_treatment_map;
  matrix<lower = 0, upper = 1>[num_obs, num_dynamic_treatment_col] dynamic_treatment_col;
  matrix<lower = 0>[num_obs * num_deworming_days, num_dynamic_treatment_col] dynamic_treatment_dm;
  
  // Constants for hyperpriors 
  real<lower = 0> scale_df;
  real<lower = 0> scale_sigma;
  real<lower = 0> scale_sigma_dynamic;
  
  real<lower = 0> coef_df;
  real<lower = 0> coef_sigma;
  real<lower = 0> coef_sigma_dynamic;
}

transformed data {
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  
  int num_not_private_value_bracelet_coef = num_all_treatment_coef - num_private_value_bracelet_coef;
}

parameters {
  row_vector<lower = 0, upper = 1>[num_deworming_days] hyper_baseline_cond_takeup; // Uniform[0, 1]
  
  vector[num_not_private_value_bracelet_coef] hyper_beta_raw; // No tau for hyper parameter; coef_sigma is the SD
  vector[num_not_private_value_bracelet_coef] stratum_beta_raw[num_strata];
  vector<lower = 0>[num_not_private_value_bracelet_coef] stratum_tau_treatment;
  
  vector[num_dynamic_treatment_col] hyper_dynamic_treatment_coef_raw;
  vector[num_dynamic_treatment_col] stratum_dynamic_treatment_coef_raw[num_strata];
  vector<lower = 0>[num_dynamic_treatment_col] stratum_tau_dynamic_treatment;
}

transformed parameters {
  row_vector<lower = 0>[num_deworming_days] hyper_baseline_hazard = - log(1 - hyper_baseline_cond_takeup);
  
  vector[num_all_treatment_coef] hyper_beta = rep_vector(0, num_all_treatment_coef); 
  vector[num_dynamic_treatment_col] hyper_dynamic_treatment_coef = hyper_dynamic_treatment_coef_raw * coef_sigma_dynamic;
  
  vector[num_strata] stratum_lp;
  
  hyper_beta[not_private_value_bracelet_coef] = hyper_beta_raw * coef_sigma;
  
  {
    int stratum_pos = 1;
    int dewormed_stratum_pos = 1;
    int dynamic_stratum_pos = 1;
    
    vector[num_obs] latent_utility = rep_vector(0, num_obs);
    matrix[num_obs, num_deworming_days] log_lambda_t;
    
    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      int curr_dewormed_stratum_size = strata_dewormed_sizes[strata_index];
      int dewormed_stratum_end = dewormed_stratum_pos + curr_dewormed_stratum_size - 1;
      
      int curr_dynamic_stratum_size = curr_stratum_size * num_deworming_days;  
      int dynamic_stratum_end = dynamic_stratum_pos + curr_dynamic_stratum_size - 1;
      
      vector[num_all_treatment_coef] stratum_beta = rep_vector(0, num_all_treatment_coef);
      
      vector[num_dynamic_treatment_col] stratum_dynamic_treatment_coef = 
        hyper_dynamic_treatment_coef + stratum_dynamic_treatment_coef_raw[strata_index] .* stratum_tau_dynamic_treatment; 
        
      int stratum_dewormed_ids[curr_dewormed_stratum_size] = dewormed_ids[dewormed_stratum_pos:dewormed_stratum_end]; 
      int stratum_dewormed_day_all[curr_stratum_size] = dewormed_day_any[stratum_pos:stratum_end]; 
      int stratum_dewormed_day_dewormed[curr_dewormed_stratum_size] = dewormed_day_any[stratum_dewormed_ids];
      
      matrix[curr_stratum_size, num_deworming_days] stratum_hazard_day_triangle_map = hazard_day_triangle_map[stratum_dewormed_day_all];
      matrix[curr_dewormed_stratum_size, num_deworming_days] stratum_hazard_day_map = hazard_day_map[stratum_dewormed_day_dewormed];
      
      stratum_beta[not_private_value_bracelet_coef] = 
        hyper_beta[not_private_value_bracelet_coef] + stratum_beta_raw[strata_index] .* stratum_tau_treatment;
              
      latent_utility[stratum_pos:stratum_end] =
          treatment_design_matrix[stratum_pos:stratum_end] * stratum_beta +
          treatment_design_matrix[stratum_pos:stratum_end, private_value_bracelet_coef] * stratum_beta[private_value_calendar_coef];
         
      // print("stratum = ", strata_index, ", latent_utility[1] = ", latent_utility[stratum_pos]);
      // print("stratum_beta = ", stratum_beta);
      // print("stratum_beta_raw = ", stratum_beta_raw);
      // print("hyper_beta = ", hyper_beta);
      
      log_lambda_t[stratum_pos:stratum_end] = 
         (rep_matrix(latent_utility[stratum_pos:stratum_end], num_deworming_days) +
          to_matrix(dynamic_treatment_dm[dynamic_stratum_pos:dynamic_stratum_end] * stratum_dynamic_treatment_coef, curr_stratum_size, num_deworming_days, 0) +
          rep_matrix(log(hyper_baseline_hazard), curr_stratum_size)); 
        
      stratum_lp[strata_index] = 
        - (exp(log_sum_exp(log_lambda_t[stratum_pos:stratum_end] .* stratum_hazard_day_triangle_map)) - sum(1 - stratum_hazard_day_triangle_map)) +
        sum(log(inv_cloglog(log_lambda_t[stratum_dewormed_ids] .* stratum_hazard_day_map))) - sum(1 - stratum_hazard_day_map) * log1m_exp(-1);
        
      if (is_nan(stratum_lp[strata_index]) || is_inf(stratum_lp[strata_index])) {
        reject("Stratum ", strata_index, ": log probability is ", stratum_lp[strata_index]);
      }
      
      // print("stratum ", strata_index, ": latent_utility[stratum_pos] = ", latent_utility[stratum_pos]);
      // print("stratum ", strata_index, ": latent_utility[stratum_end] = ", latent_utility[stratum_end]);
      // print("stratum ", strata_index, ": loglog_lambda_t[stratum_pos, ] = ", loglog_lambda_t[stratum_pos]);
      // print("stratum ", strata_index, ": loglog_lambda_t[stratum_end, ] = ", loglog_lambda_t[stratum_end]);
      // print("stratum ", strata_index, ": stratum_dewormed_day_all[1:100] = ", stratum_dewormed_day_all[1:100]);
      // print("stratum ", strata_index, ": stratum_dewormed_day_dewormed[1:10] = ", stratum_dewormed_day_dewormed[1:10]);
      // print("stratum ", strata_index, ": stratum_hazard_day_triangle_map[1:100, ] = ", stratum_hazard_day_triangle_map[1:100]);
      // print("stratum ", strata_index, ": stratum_hazard_day_map[1:10, ] = ", stratum_hazard_day_map[1:10]);
      // print("Debug value: ",
      //   hyper_baseline_cond_takeup,
      //   hyper_baseline_hazard,
      //   rep_matrix(log(hyper_baseline_hazard), 10));
        // stratum_dynamic_treatment_coef,
        // dynamic_treatment_dm[dynamic_stratum_pos:dynamic_stratum_end][1:(12 * 2)],
        // to_matrix(dynamic_treatment_dm[dynamic_stratum_pos:dynamic_stratum_end][1:(12 * 2)] * stratum_dynamic_treatment_coef, 2, num_deworming_days, 0));
        // dynamic_treatment_dm[dynamic_stratum_pos:dynamic_stratum_end][1:(12 * 10), ]);
        // to_matrix(dynamic_treatment_dm[dynamic_stratum_pos:dynamic_stratum_end] * stratum_dynamic_treatment_coef, curr_stratum_size, num_deworming_days, 0)[1:10, ]);
        // - (exp(log_sum_exp(loglog_lambda_t[stratum_pos:stratum_end] .* stratum_hazard_day_triangle_map)) - sum(1 - stratum_hazard_day_triangle_map))); 
          
      stratum_pos = stratum_end + 1;
      dewormed_stratum_pos = dewormed_stratum_end + 1;
      dynamic_stratum_pos = dynamic_stratum_end + 1;
    }
    
    // print("stratum_lp = ", stratum_lp);
  }
}

model {
  hyper_beta_raw ~ student_t(coef_df, 0, 1); 
  stratum_beta_raw ~ multi_student_t(coef_df, 
                                     rep_vector(0.0, num_not_private_value_bracelet_coef), 
                                     diag_matrix(rep_vector(1, num_not_private_value_bracelet_coef)));
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  
  hyper_dynamic_treatment_coef_raw ~ student_t(coef_df, 0, 1);
  stratum_dynamic_treatment_coef_raw ~ multi_student_t(coef_df, rep_vector(0.0, num_dynamic_treatment_col), diag_matrix(rep_vector(1, num_dynamic_treatment_col)));
  stratum_tau_dynamic_treatment ~ student_t(scale_df, 0, scale_sigma_dynamic);

  target += stratum_lp;
}

generated quantities {
  // matrix[num_all_treatments, num_deworming_days] all_treatment_cond_takeup = 1 - 
  //   exp(- exp(treatment_map_design_matrix[, not_private_value_bracelet_coef] * hyper_beta[not_private_value_bracelet_coef] + 
  //             treatment_map_design_matrix[, private_value_bracelet_coef] * hyper_beta[private_value_calendar_coef]) * 
  //         hyper_baseline_hazard');
}
