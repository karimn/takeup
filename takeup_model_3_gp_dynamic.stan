// Non parametric Guassian process for take-up dynamics

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
  
  matrix[num_deworming_days + 1, num_deworming_days] hazard_day_map;
  matrix[num_deworming_days + 1, num_deworming_days] hazard_day_triangle_map;
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
  
  // int<lower = 0, upper = num_all_treatments> num_ate_treatments;
  // int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments];
  // 
  // int<lower = 0> num_missing_obs_ids;
  // int<lower = 0, upper = num_obs> missing_treatment_sizes[num_ate_treatments];
  // int<lower = 1, upper = num_obs> missing_obs_ids[num_missing_obs_ids];
  // 
  // int<lower = 0> num_observed_obs_ids;
  // int<lower = 0, upper = num_obs> observed_treatment_sizes[num_ate_treatments];
  // int<lower = 1, upper = num_obs> observed_obs_ids[num_observed_obs_ids];
  // 
  // int<lower = 0, upper = num_all_treatments> num_ate_pairs;
  // 
  // int<lower = 1, upper = num_all_treatments> ate_pairs[num_ate_pairs, 2];
  // 
  // int missing_treatment_stratum_id[num_missing_obs_ids]; 
  // int missing_treatment_cluster_id[num_missing_obs_ids]; 
  // 
  // matrix[num_missing_obs_ids, num_census_covar_coef] missing_census_covar_dm;
  // 
  // int observed_dewormed_day[num_observed_obs_ids];
  
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
  
  matrix[num_deworming_days + 1, num_deworming_days] relevant_latent_var_map_mat = to_matrix(relevant_latent_var_map);
  int dewormed_any_daily[num_obs * num_deworming_days] = to_array_1d(relevant_latent_var_map[dewormed_day_any]);
  
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  row_vector[num_dynamic_treatment_col] nonparam_dyn_treatment_map[num_dynamic_treatments * num_deworming_days]; 
  
  {
    int non_param_dyn_treat_pos = 1;
    
    for (dyn_treat_index in 1:num_dynamic_treatments) {
      for (dyn_day_index in 1:num_deworming_days) {
        nonparam_dyn_treatment_map[non_param_dyn_treat_pos + dyn_day_index - 1] = dynamic_treatment_map[dyn_treat_index, dyn_day_index];
      }
     
      non_param_dyn_treat_pos = non_param_dyn_treat_pos + num_deworming_days;
    }
  }
}

parameters {
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
  matrix[num_dynamic_treatments, num_deworming_days] hyper_dyn_latent_var;
  
  matrix[num_strata, num_all_treatment_coef] stratum_beta_mat;
  matrix[num_strata, num_census_covar_coef] stratum_census_covar_coef_mat;
  matrix[num_strata, num_deworming_days] stratum_hazard_mat;
  
  {
    matrix[num_dynamic_treatments * num_deworming_days, num_dynamic_treatments * num_deworming_days] L_hyper_dyn_K = 
      L_cov_exp_quad_ARD(nonparam_dyn_treatment_map, hyper_dyn_alpha, hyper_dyn_rho, delta);
    
    hyper_dyn_latent_var = to_matrix(L_hyper_dyn_K * hyper_dyn_eta, num_dynamic_treatments, num_deworming_days, 0); // Row-major order
  }
 
  for (strata_index in 1:num_strata) {
    stratum_beta_mat[strata_index] = stratum_beta[strata_index]';
    stratum_census_covar_coef_mat[strata_index] = stratum_census_covar_coef[strata_index]';
  } 
}

model {
  hyper_beta ~ student_t(coef_df, 0, coef_sigma); 
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  stratum_beta ~ multi_student_t(coef_df, hyper_beta, diag_matrix(stratum_tau_treatment));
  
  hyper_census_covar_coef ~ student_t(coef_df, 0, coef_sigma);
  stratum_tau_census_covar ~ student_t(scale_df, 0, scale_sigma);
  stratum_census_covar_coef ~ multi_student_t(coef_df, hyper_census_covar_coef, diag_matrix(stratum_tau_census_covar));
 
  hyper_dyn_rho ~ gamma(4, 4); 
  hyper_dyn_alpha ~ normal(0, 1);
  hyper_dyn_eta ~ normal(0, 1);
  
  {
    int stratum_pos = 1;
    
    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      int stratum_dewormed_day_all[curr_stratum_size] = dewormed_day_any[stratum_pos:stratum_end];

      matrix[curr_stratum_size, num_deworming_days] stratum_relevant_latent_var_map = relevant_latent_var_map_mat[stratum_dewormed_day_all];
      
      matrix[curr_stratum_size, num_deworming_days] day_constant_latent_var = 
        rep_matrix(census_covar_dm[stratum_pos:stratum_end] * stratum_census_covar_coef[strata_index] + 
                     treatment_design_matrix[stratum_pos:stratum_end] * stratum_beta[strata_index], num_deworming_days); 
      
      matrix[curr_stratum_size, num_deworming_days] day_varying_latent_var = hyper_dyn_latent_var[dynamic_treatment_id[stratum_pos:stratum_end]];
      
      int stratum_daily_pos = ((stratum_pos - 1) * num_deworming_days) + 1; 
      int stratum_daily_end = stratum_daily_pos + (curr_stratum_size * num_deworming_days) - 1;
      
      dewormed_any_daily[stratum_daily_pos:stratum_daily_end] ~ bernoulli_logit(to_vector((day_constant_latent_var + day_varying_latent_var) .* stratum_relevant_latent_var_map));
          
      stratum_pos = stratum_end + 1;
    }
  }
}

generated quantities {
  /*
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
                                            missing_census_covar_dm,
                                            treatment_map_design_matrix,
                                            dynamic_treatment_map,
                                            all_treatment_dyn_id,
                                            missing_treatment_sizes,
                                            observed_treatment_sizes,
                                            stratum_hazard_mat,
                                            cluster_hazard_effect,
                                            stratum_census_covar_coef_mat,
                                            stratum_beta_mat,
                                            stratum_dyn_treatment_mat,
                                            observed_dewormed_day);
   
    est_takeup = 1 - est_deworming_days[, 13]; 
    est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }*/
}
