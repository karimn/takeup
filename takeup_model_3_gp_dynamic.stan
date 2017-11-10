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
  
  matrix treatment_deworming_day_rng(matrix full_missing_latent_var) {
    int curr_treatment_size = rows(full_missing_latent_var);
    int num_deworming_days = cols(full_missing_latent_var);
    matrix[curr_treatment_size, num_deworming_days + 1] deworming_days_mask;
    
    for (obs_ids_index in 1:curr_treatment_size) {
      vector[num_deworming_days + 1] deworming_day_prob = rep_vector(1, num_deworming_days + 1);
      vector[num_deworming_days + 1] cumul_prod_alphas = rep_vector(1, num_deworming_days + 1);
    
      real prev_not_deworm_prob = 1; 
      
      for (day_index in 1:num_deworming_days) {
        cumul_prod_alphas[day_index] = prev_not_deworm_prob; 
        deworming_day_prob[day_index] = inv_logit(full_missing_latent_var[obs_ids_index, day_index]); 
        prev_not_deworm_prob = prev_not_deworm_prob * (1 - deworming_day_prob[day_index]);
      }
      
      cumul_prod_alphas[num_deworming_days + 1] = prev_not_deworm_prob;
      
      deworming_day_prob = deworming_day_prob .* cumul_prod_alphas; 
      
      deworming_days_mask = rep_matrix(0, curr_treatment_size, num_deworming_days + 1);
      deworming_days_mask[obs_ids_index, categorical_rng(deworming_day_prob)] = 1;
    }
  
    return deworming_days_mask; 
  }
  
  matrix treatment_cell_deworming_day_prop_rng(matrix day_constant_latent_var,
                                               matrix day_varying_latent_var,
                                               int[] missing_obs_ids, 
                                               int[] missing_stratum_id,
                                               int[] missing_cluster_id,
                                               matrix missing_census_covar_dm,
                                               int[] missing_treatment_sizes,
                                               int[] observed_treatment_sizes,
                                               vector cluster_effects,
                                               matrix stratum_census_covar_coef,
                                               int[] observed_dewormed_any) {
    int num_treatment_ids = rows(day_constant_latent_var); 
    int num_deworming_days = cols(day_varying_latent_var); 
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    matrix[num_treatment_ids, num_deworming_days + 1] cell_deworming_day = rep_matrix(0, num_treatment_ids, num_deworming_days + 1);
    
    matrix[num_deworming_days + 1, num_deworming_days + 1] days_diag = diag_matrix(rep_vector(1, num_deworming_days + 1));
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
      
      matrix[curr_missing_treatment_size, num_deworming_days] full_missing_latent_var =
        rep_matrix(day_constant_latent_var[treatment_ids_index, missing_stratum_id[missing_treatment_pos:missing_treatment_end]]' +
                     rows_dot_product(stratum_census_covar_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]],
                                      missing_census_covar_dm[missing_treatment_pos:missing_treatment_end]) + 
                                      cluster_effects[missing_cluster_id[missing_treatment_pos:missing_treatment_end]],
                   num_deworming_days) +
        rep_matrix(day_varying_latent_var[treatment_ids_index], curr_missing_treatment_size);
      
      matrix[curr_missing_treatment_size, num_deworming_days + 1] missing_deworming_days_mask = treatment_deworming_day_rng(full_missing_latent_var);
      
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
  real<lower = 0, upper = 1> hyper_baseline_takeup; // Completely uninformed prior
  vector[num_strata] stratum_intercept;
  real<lower = 0> tau_stratum_intercept;
  
  vector[num_all_treatment_coef] hyper_beta; 
  vector[num_all_treatment_coef] stratum_beta[num_strata];
  vector<lower = 0>[num_all_treatment_coef] stratum_tau_treatment;
  
  vector[num_census_covar_coef] hyper_census_covar_coef;
  vector[num_census_covar_coef] stratum_census_covar_coef[num_strata];
  vector<lower = 0>[num_census_covar_coef] stratum_tau_census_covar;
  
  vector[num_clusters] cluster_effects;
  real<lower = 0> tau_cluster_effect;
  
  //  Dynamics
  
  row_vector<lower=0>[num_dynamic_treatment_col] hyper_dyn_rho;
  real<lower=0> hyper_dyn_alpha;
  vector[num_dynamic_treatments * num_deworming_days] hyper_dyn_eta; 
}

transformed parameters {
  real hyper_intercept = logit(hyper_baseline_takeup);
  matrix[num_dynamic_treatments, num_deworming_days] hyper_dyn_latent_var;
  
  matrix[num_all_treatment_coef, num_strata] stratum_beta_mat;
  matrix[num_census_covar_coef, num_strata] stratum_census_covar_coef_mat;
  
  {
    matrix[num_dynamic_treatments * num_deworming_days, num_dynamic_treatments * num_deworming_days] L_hyper_dyn_K = 
      L_cov_exp_quad_ARD(nonparam_dyn_treatment_map, hyper_dyn_alpha, hyper_dyn_rho, delta);
    
    hyper_dyn_latent_var = to_matrix(L_hyper_dyn_K * hyper_dyn_eta, num_dynamic_treatments, num_deworming_days, 0); // Row-major order
  }
 
  for (strata_index in 1:num_strata) {
    stratum_beta_mat[, strata_index] = stratum_beta[strata_index];
    stratum_census_covar_coef_mat[, strata_index] = stratum_census_covar_coef[strata_index];
  } 
}

model {
  tau_stratum_intercept ~ student_t(scale_df, 0, scale_sigma);
  stratum_intercept ~ student_t(coef_df, hyper_intercept, tau_stratum_intercept); 
  
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
    int stratum_pos = 1;
    
    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      int stratum_dewormed_day_all[curr_stratum_size] = dewormed_day_any[stratum_pos:stratum_end];

      matrix[curr_stratum_size, num_deworming_days] stratum_relevant_latent_var_map = relevant_latent_var_map_mat[stratum_dewormed_day_all];
      
      matrix[curr_stratum_size, num_deworming_days] day_constant_latent_var = 
        stratum_intercept[strata_index] +
        rep_matrix(census_covar_dm[stratum_pos:stratum_end] * stratum_census_covar_coef[strata_index] + 
                     treatment_design_matrix[stratum_pos:stratum_end] * stratum_beta[strata_index] +
                     cluster_effects[cluster_id[stratum_pos:stratum_end]], 
                   num_deworming_days); 
      
      matrix[curr_stratum_size, num_deworming_days] day_varying_latent_var = hyper_dyn_latent_var[dynamic_treatment_id[stratum_pos:stratum_end]];
      
      int stratum_daily_pos = ((stratum_pos - 1) * num_deworming_days) + 1; 
      int stratum_daily_end = stratum_daily_pos + (curr_stratum_size * num_deworming_days) - 1;
      
      dewormed_any_daily[stratum_daily_pos:stratum_daily_end] ~ bernoulli_logit(to_vector((day_constant_latent_var + day_varying_latent_var) .* stratum_relevant_latent_var_map));
          
      stratum_pos = stratum_end + 1;
    }
  }
}

generated quantities {
  matrix[num_ate_treatments, num_strata] day_constant_latent_var = rep_matrix(0, num_ate_treatments, num_strata);
  matrix[num_ate_treatments, num_deworming_days] day_varying_latent_var = rep_matrix(0, num_ate_treatments, num_deworming_days);
  
  matrix<lower = 0>[num_ate_treatments, num_deworming_days + 1] est_deworming_days = rep_matrix(0, num_ate_treatments, num_deworming_days + 1);
  vector<lower = 0, upper = 1>[num_ate_treatments] est_takeup = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] est_takeup_ate = rep_vector(0, num_ate_pairs);
  
  if (estimate_ate) {
    day_constant_latent_var = treatment_map_design_matrix[ate_treatments] * stratum_beta_mat + rep_matrix(stratum_intercept', num_ate_treatments);
    day_varying_latent_var = hyper_dyn_latent_var[all_treatment_dyn_id[ate_treatments]];
    
    est_deworming_days =
      treatment_cell_deworming_day_prop_rng(day_constant_latent_var,
                                            day_varying_latent_var,
                                            missing_obs_ids,
                                            missing_treatment_stratum_id,
                                            missing_treatment_cluster_id,
                                            missing_census_covar_dm,
                                            missing_treatment_sizes,
                                            observed_treatment_sizes,
                                            cluster_effects,
                                            stratum_census_covar_coef_mat,
                                            observed_dewormed_day);
   
    est_takeup = 1 - est_deworming_days[, 13]; 
    est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }
}
