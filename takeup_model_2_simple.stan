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
  
  // Constants for hyperpriors 
  real<lower = 0> scale_df;
  real<lower = 0> scale_sigma;
  
  real<lower = 0> coef_df;
  real<lower = 0> coef_sigma;
}

transformed data {
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  
  int num_not_private_value_bracelet_coef = num_all_treatment_coef - num_private_value_bracelet_coef;
}

parameters {
  vector<lower = 0, upper = 1>[num_deworming_days] hyper_baseline_cond_takeup; // Uniform[0, 1]
  
  vector[num_not_private_value_bracelet_coef] hyper_beta_raw; // No tau for hyper parameter; coef_sigma is the SD
  vector[num_not_private_value_bracelet_coef] stratum_beta_raw[num_strata];
  vector<lower = 0>[num_not_private_value_bracelet_coef] stratum_tau_treatment;
}

transformed parameters {
  vector<lower = 0>[num_deworming_days] hyper_baseline_hazard = - log(1 - hyper_baseline_cond_takeup);
  
  vector[num_all_treatment_coef] hyper_beta = rep_vector(0, num_all_treatment_coef); 
  
  vector[num_obs] latent_utility = rep_vector(0, num_obs);
  vector<upper = 0>[num_obs] hetero_kappa = rep_vector(0, num_obs);
  
  vector[num_strata] stratum_lp;
  
  hyper_beta[not_private_value_bracelet_coef] = hyper_beta_raw * coef_sigma;
  
  {
    int stratum_pos = 1;
    int dewormed_stratum_pos = 1;
    
    vector[num_deworming_days + 1] stratum_triangle_sum_lambda = hazard_day_triangle_map * hyper_baseline_hazard;
    vector[num_deworming_days] stratum_dewormed_day_lambda = hazard_day_map * hyper_baseline_hazard;

    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      int curr_dewormed_stratum_size = strata_dewormed_sizes[strata_index];
      int dewormed_stratum_end = dewormed_stratum_pos + curr_dewormed_stratum_size - 1;
      
      vector[num_all_treatment_coef] stratum_beta = rep_vector(0, num_all_treatment_coef);
      
      stratum_beta[not_private_value_bracelet_coef] = hyper_beta[not_private_value_bracelet_coef] + stratum_beta_raw[strata_index] .* stratum_tau_treatment;
              
      latent_utility[stratum_pos:stratum_end] =
          treatment_design_matrix[stratum_pos:stratum_end] * stratum_beta + 
          treatment_design_matrix[stratum_pos:stratum_end, private_value_bracelet_coef] * stratum_beta[private_value_calendar_coef];
         
      hetero_kappa[stratum_pos:stratum_end] = 
        (- exp(latent_utility[stratum_pos:stratum_end])); 
    
      stratum_lp[strata_index] = 
        sum(hetero_kappa[stratum_pos:stratum_end] .* stratum_triangle_sum_lambda[dewormed_day_any[stratum_pos:stratum_end]]) +
        sum(log1m_exp(hetero_kappa[dewormed_ids[dewormed_stratum_pos:dewormed_stratum_end]] .*
                        stratum_dewormed_day_lambda[dewormed_day_any[dewormed_ids[dewormed_stratum_pos:dewormed_stratum_end]]]));
      
      stratum_pos = stratum_end + 1;
      dewormed_stratum_pos = dewormed_stratum_end + 1;
    }
  }
}

model {
  hyper_beta_raw ~ student_t(coef_df, 0, 1); 
  stratum_beta_raw ~ multi_student_t(coef_df, 
                                     rep_vector(0.0, num_not_private_value_bracelet_coef), 
                                     diag_matrix(rep_vector(1, num_not_private_value_bracelet_coef)));
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);

  target += sum(stratum_lp);
}

generated quantities {
  matrix[num_all_treatments, num_deworming_days] all_treatment_cond_takeup = 1 - 
    exp(- exp(treatment_map_design_matrix[, not_private_value_bracelet_coef] * hyper_beta[not_private_value_bracelet_coef] + 
              treatment_map_design_matrix[, private_value_bracelet_coef] * hyper_beta[private_value_calendar_coef]) * 
          hyper_baseline_hazard');
  }
