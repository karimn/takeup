functions {
  real[] takeup_proportion_rng(int[] eval_treatment_ids, // the treatment IDs to evaluate proportions for  
                               matrix treatment_map_design_matrix,
                               matrix name_match_interact_map_design_matrix,
                               int[] dewormed_any,
                               vector[] beta,
                               vector stratum_intercept,
                               vector cluster_effects,
                               int[] treatment_id,
                               int[] missing_treatment_id,
                               int[] treatment_sizes,
                               int[] missing_treatment_sizes,
                               int[] missing_stratum_id_individ,
                               int[] missing_cluster_id) { 
    int num_eval_treatment_prop = size(eval_treatment_ids);
    int total_num_coef = cols(treatment_map_design_matrix) + cols(name_match_interact_map_design_matrix);
    real takeup_proportion[num_eval_treatment_prop]; 
    int num_strata = num_elements(stratum_intercept);
    real missing_link_model[num_strata, num_eval_treatment_prop];
    
    for (strata_index in 1:num_strata) {
      matrix[num_eval_treatment_prop, total_num_coef] all_treat_map_design_matrix = 
        append_col(treatment_map_design_matrix[eval_treatment_ids], name_match_interact_map_design_matrix[eval_treatment_ids]); 
        
      missing_link_model[strata_index, ] = to_array_1d(all_treat_map_design_matrix * beta[strata_index] + stratum_intercept[strata_index]);
    }
    
    for (treatment_index_index in 1:num_eval_treatment_prop) {
      int treatment_index = eval_treatment_ids[treatment_index_index]; 
      int treatment_pos;  
      int missing_treatment_pos; 
      int curr_treatment_size = treatment_sizes[treatment_index];
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_index];
      int treatment_end;
      int missing_treatment_end;
      int missing_stratum_id[curr_missing_treatment_size];
      int treatment_missing_cluster_id[curr_missing_treatment_size];
      real missing_treatment_link[curr_missing_treatment_size];
      real missing_dewormed[curr_missing_treatment_size];
      int curr_missing_id[curr_missing_treatment_size];
      
      if (treatment_index > 1) {
        treatment_pos = sum(treatment_sizes[1:(treatment_index - 1)]) + 1;
        missing_treatment_pos = sum(missing_treatment_sizes[1:(treatment_index - 1)]) + 1;
      } else {
        treatment_pos = 1;
        missing_treatment_pos = 1;
      }
      
      treatment_end = treatment_pos + curr_treatment_size - 1;
      missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      curr_missing_id = missing_treatment_id[missing_treatment_pos:missing_treatment_end];
      missing_stratum_id = missing_stratum_id_individ[curr_missing_id];
      treatment_missing_cluster_id = missing_cluster_id[curr_missing_id];
      missing_treatment_link = missing_link_model[missing_stratum_id, treatment_index_index];
      
      for(i in 1:curr_missing_treatment_size) {
        missing_dewormed[i] = bernoulli_rng(Phi(missing_treatment_link[i] + cluster_effects[treatment_missing_cluster_id[i]]));
      }
      
      takeup_proportion[treatment_index_index] = 
        (sum(missing_dewormed) + sum(dewormed_any[treatment_id[treatment_pos:treatment_end]])) /
        (curr_missing_treatment_size + curr_treatment_size);
    }
    
    return(takeup_proportion);
  }

  vector row_sums(matrix to_sum) {
    return(to_sum * rep_vector(1, rows(to_sum)));
  }
}

data {
  int<lower = 0> num_obs;
  int<lower = 1> num_all_treatments; // # of treatment cells
  int<lower = 1> num_all_treatment_coef; // # of linear model coefficients minus intercept
  int<lower = 0> num_experiment_coef;
  int<lower = 1> num_name_match_interact_coef; // # of linear model for name matched indicator interacted with other treatments
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  // Below references to "maps" refer to a finite set of treatment cells that has a one-to-many relationship to actual observations
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix; // Design matrix generated from treatment map
  matrix[num_all_treatments, num_name_match_interact_coef] name_match_interact_map_design_matrix; // Design matrix generated from name-match interaction map
  int<lower = 1, upper = num_all_treatment_coef> experiment_coef[num_experiment_coef];
  int<lower = 1, upper = num_all_treatments> obs_treatment[num_obs]; // ID of observed treatment (from treatment map)
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1> strata_sizes[num_strata]; 
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  int<lower = 1> treatment_sizes[num_all_treatments]; // Number of observations per treatment ID, in ascending order
 
  // Counterfactual information for finite sample analysis
  
  int<lower = 1> num_missing; // # of missing counterfactuals
  int<lower = 1> num_eval_treatment_prop; // # of treatment cells to be imputed
  
  int<lower = 1, upper = num_obs> stratum_covar_id[num_obs]; // Observation indices ordered by stratum and endline covariate missingness
  int<lower = 1> stratum_missing_covar_sizes[num_strata]; // Number of covariate missingness per stratum (ordered by stratum)
  int<lower = 1, upper = num_all_treatments> eval_treatment_prop_id[num_eval_treatment_prop]; // Which treatments to impute
  int<lower = 1, upper = num_missing> missing_treatment_id[num_missing]; // Observation indices with missing treatment counterfactural, ordered by missing treatment ID
  int<lower = 1> missing_treatment_sizes[num_all_treatments]; // Number of obserations missing treatment, in ascending order of missing treatment ID
  
  // Covariates
  
  int<lower = 1> num_census_covar_coef;
  int<lower = 1> num_endline_covar_coef; 
  int<lower = 1> num_distinct_census_covar;
  
  matrix[num_distinct_census_covar, num_census_covar_coef] census_covar_map_dm; 
  int census_covar_id[num_obs];
  // matrix[num_obs - sum(stratum_missing_covar_sizes), num_endline_covar_coef] endline_covar_dm;
 
  // Binary deworming outcome 
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
}

transformed data {
  // int<lower = 1, upper = num_strata> missing_stratum_id[num_missing] = stratum_id[missing_treatment_id];
  // int<lower = 1, upper = num_clusters> missing_cluster_id[num_missing] = cluster_id[missing_treatment_id];
  // 
  // int<lower = 1, upper = num_obs> num_missing_endline_covar = sum(stratum_missing_covar_sizes);
  
  // Unmodeled parameters for priors and hyperpriors
  
  vector[num_all_treatment_coef] mu = rep_vector(0, num_all_treatment_coef);
  vector[num_census_covar_coef] mu_census_covar = rep_vector(0, num_census_covar_coef);
  // vector[num_endline_covar_coef] mu_endline_covar = rep_vector(0, num_endline_covar_coef);
  
  vector[num_all_treatment_coef] tau_treatment = rep_vector(1, num_all_treatment_coef);
  // vector[num_name_match_interact_coef] tau_name_match_interact = rep_vector(sqrt(0.25), num_name_match_interact_coef);
  // vector[num_treat_name_match_coef] tau = append_row(tau_treatment, tau_name_match_interact);
  cov_matrix[num_all_treatment_coef] Sigma_beta = diag_matrix(tau_treatment);
  
  vector[num_census_covar_coef] tau_census_covar = rep_vector(1, num_census_covar_coef);
  cov_matrix[num_census_covar_coef] Sigma_census_covar = diag_matrix(tau_census_covar);
  // vector[num_endline_covar_coef] tau_endline_covar = rep_vector(1, num_endline_covar_coef);
  // cov_matrix[num_endline_covar_coef] Sigma_endline_covar = diag_matrix(tau_endline_covar);
}

parameters {
  // Modelled parameters
  
  real<lower=-3, upper=3> mu_strata; // intercept hyperparameter, uniformly distributed prior
  vector[num_strata] stratum_intercept; // Stratum variation in intercept
  
  vector[num_clusters] cluster_effects;
  
  vector[num_all_treatment_coef] hyper_beta;
  vector[num_all_treatment_coef] beta[num_strata];
  vector[num_name_match_interact_coef] beta_name_match_interact;
  
  vector[num_census_covar_coef] hyper_census_covar_coef;
  vector[num_census_covar_coef] census_covar_coef[num_strata];
  
  matrix[num_census_covar_coef, num_experiment_coef] hyper_beta_census_covar_interact; // census characteristics' effect on treatment effects
  // matrix[num_census_covar_coef, num_experiment_coef] beta_census_covar_interact[num_strata]; 
  // vector[num_endline_covar_coef] hyper_endline_covar_coef;
  // vector[num_endline_covar_coef] endline_covar_coef[num_strata];
  
  // Scale hyperparameters for betas 
  
  // vector[num_strata] stratum_intercept_covar;
  
  // corr_matrix[num_treatments] Omega_beta;
  // vector<lower = 0>[num_all_treatment_coef] tau_beta;
  // real<lower = 0> tau_stratum_effect;
  // real<lower = 0> tau_cluster_effect;
  
  // Parameters used to predict missing endline covariates
  
  // Imputed data
  
  // matrix[num_missing_endline_covar, num_endline_covar_coef] missing_endline_covar; 
}

transformed parameters {
}

model {
  int strata_pos = 1;
  matrix[num_distinct_census_covar, num_experiment_coef] census_covar_interact_map_dm; 
  
  hyper_census_covar_coef ~ multi_normal(mu_census_covar, Sigma_census_covar);
  census_covar_coef ~ multi_normal(hyper_census_covar_coef, Sigma_census_covar);
  
  to_vector(hyper_beta_census_covar_interact) ~ normal(0, 1);
  
  // hyper_endline_covar_coef ~ multi_normal(mu_endline_covar, Sigma_endline_covar);
  // endline_covar_coef ~ multi_normal(hyper_endline_covar_coef, Sigma_endline_covar);
  
  stratum_intercept ~ normal(mu_strata, 1); # tau_stratum_effect); 
  cluster_effects ~ normal(0, 1); # tau_cluster_effect);
  
  beta_name_match_interact ~ normal(0, 1);
  
  hyper_beta ~ multi_normal(mu, Sigma_beta); // For now assuming no correlation between effects
  beta ~ multi_normal(hyper_beta, Sigma_beta);
  
  census_covar_interact_map_dm = census_covar_map_dm * hyper_beta_census_covar_interact;
  
  for (strata_index in 1:num_strata) {
    int curr_stratum_size = strata_sizes[strata_index];
    int stratum_end = strata_pos + curr_stratum_size - 1;
    
    // matrix[curr_stratum_size, num_all_treatment_coef] stratum_treat_map_dm = treatment_map_design_matrix[obs_treatment[strata_pos:stratum_end]];
    // matrix[curr_stratum_size, num_census_covar_coef] stratum_census_covar_map_dm = census_covar_map_dm[census_covar_id[strata_pos:stratum_end]];
    
    vector[curr_stratum_size] test = rows_dot_product(census_covar_interact_map_dm[census_covar_id[strata_pos:stratum_end]], 
                                                      treatment_map_design_matrix[obs_treatment[strata_pos:stratum_end], experiment_coef]);
  
    // Probit
    dewormed_any[strata_pos:stratum_end] ~ bernoulli(Phi(
      stratum_intercept[strata_index] + 
      cluster_effects[cluster_id[strata_pos:stratum_end]] + 
      treatment_map_design_matrix[obs_treatment[strata_pos:stratum_end]] * beta[strata_index] +
      // stratum_treat_map_dm * beta[strata_index] + 
      name_match_interact_map_design_matrix[obs_treatment[strata_pos:stratum_end]] * beta_name_match_interact +
      census_covar_map_dm[census_covar_id[strata_pos:stratum_end]] * census_covar_coef[strata_index] 
    )); 
    
    strata_pos = strata_pos + curr_stratum_size;
  }
}

generated quantities {
  // real<lower = 0, upper = 1> takeup_proportion[num_eval_treatment_prop] =
  //   takeup_proportion_rng(
  //     eval_treatment_prop_id, treatment_map_design_matrix, name_match_interact_map_design_matrix, dewormed_any, beta, stratum_intercept, cluster_effects,
  //     treatment_id, missing_treatment_id, treatment_sizes, missing_treatment_sizes, missing_stratum_id, missing_cluster_id
  //   );
}
