functions {
}

data {
  int<lower = 0> num_obs;
  int<lower = 1> num_missing; // # of missing counterfactuals
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  int<lower = 1, upper = num_obs> stratum_covar_id[num_obs]; // Observation indices ordered by stratum and endline covariate missingness
  int<lower = 1> stratum_missing_covar_sizes[num_strata]; // Number of covariate missingness per stratum (ordered by stratum)
  
  // Below references to "maps" refer to a finite set of treatment cells that has a one-to-many relationship to actual observations
  
  // # of linear model of covariates coefficients minus intercept
  int<lower = 1> num_census_covar_coef; 
  int<lower = 1> num_endline_covar_coef; 
  
  // Design matrix for covariates
  matrix[num_obs, num_census_covar_coef] census_covar_dm;
  //matrix[num_obs - sum(stratum_missing_sizes), num_endline_covar_coef] endline_covar_dm;
  vector[num_endline_covar_coef] endline_covar_dm[num_obs - sum(stratum_missing_covar_sizes)];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1> strata_sizes[num_strata]; 
}

transformed data {
  int<lower = 1, upper = num_obs> num_missing_endline_covar = sum(stratum_missing_covar_sizes);
  int<lower = 1, upper = num_obs> num_obs_endline_covar = num_obs - num_missing_endline_covar;
}

parameters {
  // matrix[num_endline_covar_coef, num_census_covar_coef] beta_census_covar[num_strata];
  matrix[num_census_covar_coef, num_endline_covar_coef] beta_census_covar[num_strata];
  cholesky_factor_corr[num_endline_covar_coef] L_Omega_census_covar[num_strata];
  vector<lower = 0>[num_endline_covar_coef] L_sigma_census_covar[num_strata];
}

transformed parameters {
}

model {
  matrix[num_endline_covar_coef, num_endline_covar_coef] L_Sigma_census_covar[num_strata];

  int strata_pos = 1;
  int stratum_obs_start = 1;

  for (stratum_index in 1:num_strata) {
    int curr_stratum_size = strata_sizes[stratum_index];
    int curr_missing_endline_size = stratum_missing_covar_sizes[stratum_index];
    int curr_endline_size = curr_stratum_size - curr_missing_endline_size;
    int stratum_end = strata_pos + curr_endline_size - 1;
    int stratum_obs_end = stratum_obs_start + curr_endline_size - 1;
    int curr_obs_index[curr_endline_size] = stratum_covar_id[strata_pos:stratum_end];
    
    row_vector[num_endline_covar_coef] mu[curr_endline_size];
    
    to_vector(beta_census_covar[stratum_index]) ~ normal(0, 2);
    
    L_Omega_census_covar[stratum_index] ~ lkj_corr_cholesky(4);
    L_sigma_census_covar[stratum_index] ~ cauchy(0, 2.5);

    L_Sigma_census_covar[stratum_index] = diag_pre_multiply(L_sigma_census_covar[stratum_index], L_Omega_census_covar[stratum_index]);
   
    // print(beta_census_covar[stratum_index]); 
    
    for (i in 1:curr_endline_size) {
      mu[i] = census_covar_dm[curr_obs_index[i]] * beta_census_covar[stratum_index];
    }
    // 
    endline_covar_dm[stratum_obs_start:stratum_obs_end] ~ multi_normal_cholesky(mu, L_Sigma_census_covar[stratum_index]);

    strata_pos = strata_pos + curr_stratum_size;
    stratum_obs_start = stratum_obs_end + 1;
  }
}

generated quantities {
}
