functions {/*
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
    // int total_num_coef = cols(treatment_map_design_matrix) + cols(name_match_interact_map_design_matrix);
    real takeup_proportion[num_eval_treatment_prop]; 
    int num_strata = num_elements(stratum_intercept);
    // real missing_link_model[num_strata, num_eval_treatment_prop];
    
    // for (strata_index in 1:num_strata) {
    //   matrix[num_eval_treatment_prop, total_num_coef] all_treat_map_design_matrix = 
    //     append_col(treatment_map_design_matrix[eval_treatment_ids], name_match_interact_map_design_matrix[eval_treatment_ids]); 
    //     
    //   missing_link_model[strata_index, ] = to_array_1d(all_treat_map_design_matrix * beta[strata_index] + stratum_intercept[strata_index]);
    // }
    
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
      } else
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
*/
  real dewormed_matched_lpmf(int[] matched, int[] dewormed, real prob_false_pos, real prob_false_neg) {
    return(bernoulli_lpmf(matched | prob_false_pos + (1 - prob_false_pos - prob_false_neg) * to_vector(dewormed)));
  }

  real dewormed_monitored_probit_lpmf(int[] dewormed, real stratum_intercept, vector cluster_effects, vector treatment_effects, vector census_covar_effects) {
    return(bernoulli_lpmf(dewormed | Phi_approx(stratum_intercept + cluster_effects + treatment_effects + census_covar_effects)));
  }  
}

data {
  int<lower = 0> num_obs;
  int<lower = 1> num_all_treatments; // # of treatment cells
  int<lower = 1> num_all_treatment_coef; // # of linear model coefficients minus intercept
  int<lower = 0> num_experiment_coef;
  // int<lower = 1> num_name_match_interact_coef; // # of linear model for name matched indicator interacted with other treatments
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  // Constants for priors
  
  real<lower = 0> name_match_false_pos_alpha;
  real<lower = 0> name_match_false_pos_beta;
  
  // Below references to "maps" refer to a finite set of treatment cells that has a one-to-many relationship to actual observations
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix; // Design matrix generated from treatment map
  // matrix[num_all_treatments, num_name_match_interact_coef] name_match_interact_map_design_matrix; // Design matrix generated from name-match interaction map
  int<lower = 1, upper = num_all_treatment_coef> experiment_coef[num_experiment_coef];
  int<lower = 1, upper = num_all_treatments> obs_treatment[num_obs]; // ID of observed treatment (from treatment map)
  
  int <lower = 0, upper = num_obs> num_name_matched;
  int <lower = 0, upper = num_obs> num_name_matching_errors_ids;
  vector<lower = 0, upper = 1>[num_obs] name_matched;
  int<lower = 1, upper = num_obs> name_matched_id[num_name_matched];
  int<lower = 1, upper = num_obs> name_matched_strata_sizes[num_strata];
  int<lower = 1, upper = num_obs> name_matched_dewormed_strata_sizes[num_strata];
  int<lower = 1, upper = num_obs> name_matching_error_ids[num_name_matching_errors_ids];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  int<lower = 1, upper = num_obs> treatment_sizes[num_all_treatments]; // Number of observations per treatment ID, in ascending order
  // int<lower = 1, upper = num_all_treatments> num_non_phone_owner_treatments;
 
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
  
  vector<lower = 0, upper = 1>[num_obs] monitored = 1 - name_matched;
  
  // Unmodeled parameters for priors and hyperpriors
  
  vector[num_all_treatment_coef] mu = rep_vector(0, num_all_treatment_coef);
  vector[num_census_covar_coef] mu_census_covar = rep_vector(0, num_census_covar_coef);
  // vector[num_endline_covar_coef] mu_endline_covar = rep_vector(0, num_endline_covar_coef);
  
  // vector[num_all_treatment_coef] tau_treatment = rep_vector(1, num_all_treatment_coef);
  // vector[num_name_match_interact_coef] tau_name_match_interact = rep_vector(sqrt(0.25), num_name_match_interact_coef);
  // vector[num_treat_name_match_coef] tau = append_row(tau_treatment, tau_name_match_interact);
  // cov_matrix[num_all_treatment_coef] Sigma_beta = diag_matrix(tau_treatment);
  
  vector[num_census_covar_coef] tau_census_covar = rep_vector(1, num_census_covar_coef);
  cov_matrix[num_census_covar_coef] Sigma_census_covar = diag_matrix(tau_census_covar);
  // vector[num_endline_covar_coef] tau_endline_covar = rep_vector(1, num_endline_covar_coef);
  // cov_matrix[num_endline_covar_coef] Sigma_endline_covar = diag_matrix(tau_endline_covar);
}

parameters {
  // Modelled parameters
  
  // real<lower=-3, upper=3> mu_strata; // intercept hyperparameter, uniformly distributed prior
  // vector<lower = 0, upper = 1>[num_strata] stratum_baseline_takeup; // Stratum variation in intercept
  real<lower = 0, upper = 1> hyper_baseline_takeup; 
  vector[num_strata] stratum_intercept;
  real<lower = 0> tau_stratum_intercept;
  
  vector[num_clusters] cluster_effects;
  real<lower = 0> tau_cluster_effect;
  
  vector[num_all_treatment_coef] hyper_beta;
  vector[num_all_treatment_coef] beta[num_strata];
  vector<lower = 0>[num_all_treatment_coef] tau_treatment;
  
  // vector[num_name_match_interact_coef] beta_name_match_interact;
  
  vector[num_census_covar_coef] hyper_census_covar_coef;
  vector[num_census_covar_coef] census_covar_coef[num_strata];
  
  // matrix[num_census_covar_coef, num_experiment_coef] hyper_beta_census_covar_interact; // census characteristics' effect on treatment effects
  // matrix[num_census_covar_coef, num_experiment_coef] beta_census_covar_interact[num_strata]; 
  // vector[num_endline_covar_coef] hyper_endline_covar_coef;
  // vector[num_endline_covar_coef] endline_covar_coef[num_strata];
  
  // Scale hyperparameters for betas 
  
  // vector[num_strata] stratum_intercept_covar;
  
  // corr_matrix[num_treatments] Omega_beta;
  // vector<lower = 0>[num_all_treatment_coef] tau_beta;
  // real<lower = 0> tau_stratum_effect;
  
  // Parameters used to predict missing endline covariates
  
  // Imputed data
  
  // matrix[num_missing_endline_covar, num_endline_covar_coef] missing_endline_covar; 
  
  // Name matching model parameters
  
  // real<lower = 0, upper = 0.5> name_match_false_pos; 
  // real<lower = 0, upper = 1> name_match_false_neg; // uninformative prior on this
  // 
  // real<lower = 0, upper = 1> name_match_error_intercept;
}

transformed parameters {
  // real<lower = 0, upper = 1> name_match_error_diff = name_match_error_intercept * (name_match_false_neg - name_match_false_pos);
  // vector[num_strata] stratum_intercept = inv_Phi(stratum_baseline_takeup);
  real hyper_intercept = inv_Phi(hyper_baseline_takeup);
}

model {
  int strata_pos = 1;
  
  stratum_intercept ~ normal(0, 1); 
  tau_stratum_intercept ~ normal(0, 1);
  
  cluster_effects ~ normal(0, 1); 
  tau_cluster_effect ~ normal(0, 1);
  
  hyper_beta ~ normal(0, 10); # diag_matrix(rep_vector(2, num_all_treatment_coef))); // For now assuming no correlation between effects
  // hyper_beta ~ multi_normal(mu, diag_matrix(rep_vector(2, num_all_treatment_coef))); // For now assuming no correlation between effects
  beta ~ multi_normal(rep_vector(0.0, num_all_treatment_coef), diag_matrix(rep_vector(1, num_all_treatment_coef)));
  // beta ~ multi_normal(hyper_beta, rep_vector(1, num_all_treatment_coef));
  tau_treatment ~ normal(0, 20);
  
  // matrix[num_distinct_census_covar, num_experiment_coef] census_covar_interact_map_dm; 
  // real name_matched_lp[2, num_strata]; // binary (latent) deworming outcome
  
  // name_match_false_pos ~ beta(name_match_false_pos_alpha, name_match_false_pos_beta);
  
  // hyper_census_covar_coef ~ multi_normal(mu_census_covar, Sigma_census_covar);
  hyper_census_covar_coef ~ normal(0, 10);
  census_covar_coef ~ multi_normal(rep_vector(0.0, num_census_covar_coef), Sigma_census_covar);
  
  // to_vector(hyper_beta_census_covar_interact) ~ normal(0, 1);
  
  // hyper_endline_covar_coef ~ multi_normal(mu_endline_covar, Sigma_endline_covar);
  // endline_covar_coef ~ multi_normal(hyper_endline_covar_coef, Sigma_endline_covar);
  
  // census_covar_interact_map_dm = census_covar_map_dm * hyper_beta_census_covar_interact;
  
  for (strata_index in 1:num_strata) {
    int curr_stratum_size = strata_sizes[strata_index];
    int curr_matched_stratum_size = name_matched_strata_sizes[strata_index];
    int curr_monitored_stratum_size = curr_stratum_size - curr_matched_stratum_size;
    int curr_matched_dewormed_stratum_size = name_matched_dewormed_strata_sizes[strata_index];
    int stratum_end = strata_pos + curr_stratum_size - 1;
    int strata_matched_pos = strata_pos + curr_monitored_stratum_size;
    int strata_matched_not_dewormed_end = stratum_end - curr_matched_dewormed_stratum_size;
    int strata_matched_dewormed_pos = strata_matched_not_dewormed_end + 1;
    int monitored_stratum_end = strata_pos + curr_monitored_stratum_size - 1;
    
    // print("stratum[", strata_index, "]: stratum size = ", curr_stratum_size, ", strata_pos = ", strata_pos, 
    //       ", monitored_stratum_end = ", monitored_stratum_end, ", strata_matched_pos = ", strata_matched_pos, ", stratum_end = ", stratum_end);
    // print("sum monitored:", sum(name_matched[strata_pos:monitored_stratum_end]));
    // print("sum matched:", sum(name_matched[strata_matched_pos:stratum_end]));
    // print("size of monitored: ", curr_monitored_stratum_size);
    // print("size of matched: ", curr_matched_stratum_size);
    
    // vector[curr_stratum_size] test = rows_dot_product(census_covar_interact_map_dm[census_covar_id[strata_pos:stratum_end]], 
    //                                                   treatment_map_design_matrix[obs_treatment[strata_pos:stratum_end], experiment_coef]);
  
    target += dewormed_monitored_probit_lpmf(
      dewormed_any[strata_pos:monitored_stratum_end] | hyper_intercept + stratum_intercept[strata_index] * tau_stratum_intercept, 
      
                                                       cluster_effects[cluster_id[strata_pos:monitored_stratum_end]] * tau_cluster_effect,
                                                       
                                                       // treatment_map_design_matrix[obs_treatment[strata_pos:monitored_stratum_end]] * (hyper_beta + beta[strata_index]),
                                                       treatment_map_design_matrix[obs_treatment[strata_pos:monitored_stratum_end]] * (hyper_beta + beta[strata_index] .* tau_treatment),
                                                         
                                                       census_covar_map_dm[census_covar_id[strata_pos:monitored_stratum_end]] * (hyper_census_covar_coef + census_covar_coef[strata_index]));
                                                       
    // for (true_dewormed in 0:1) {
    // name_matched_lp[true_dewormed + 1, strata_index] =
    //   dewormed_monitored_probit_lpmf(rep_array(true_dewormed, curr_matched_stratum_size) | hyper_intercept + stratum_intercept[strata_index] * tau_stratum_intercept, 
    //   
    //                                                                                        cluster_effects[cluster_id[strata_matched_pos:stratum_end]] * tau_cluster_effect,
    //                                                                                        
    //                                                                                        treatment_map_design_matrix[obs_treatment[strata_matched_pos:stratum_end]] * beta[strata_index],
    //                                                                                        
    //                                                                                        census_covar_map_dm[census_covar_id[strata_matched_pos:stratum_end]] * census_covar_coef[strata_index]) +
    //                                                                                        
    //   dewormed_matched_lpmf(dewormed_any[strata_matched_pos:stratum_end] | rep_array(true_dewormed, curr_matched_stratum_size), name_match_false_pos, name_match_false_neg);
    // }
               
    // Marginalizing out the imputed deworming outcome (for the name matched)           
    // name_matched_lp[1, strata_index] = 
    //   dewormed_monitored_probit_lpmf(rep_array(0, curr_matched_stratum_size) | stratum_intercept[strata_index], 
    //                                                                            cluster_effects[cluster_id[strata_matched_pos:stratum_end]],
    //                                                                            treatment_map_design_matrix[obs_treatment[strata_matched_pos:stratum_end]] * beta[strata_index],
    //                                                                            census_covar_map_dm[census_covar_id[strata_matched_pos:stratum_end]] * census_covar_coef[strata_index]) +
    //   dewormed_matched_lpmf(dewormed_any[strata_matched_pos:stratum_end] | rep_array(0, curr_matched_stratum_size), name_match_false_pos, name_match_false_neg);
    //   
    // name_matched_lp[2, strata_index] = 
    //   dewormed_monitored_probit_lpmf(rep_array(1, curr_matched_stratum_size) | stratum_intercept[strata_index], 
    //                                                                            cluster_effects[cluster_id[strata_matched_pos:stratum_end]],
    //                                                                            treatment_map_design_matrix[obs_treatment[strata_matched_pos:stratum_end]] * beta[strata_index],
    //                                                                            census_covar_map_dm[census_covar_id[strata_matched_pos:stratum_end]] * census_covar_coef[strata_index]) +
    //   dewormed_matched_lpmf(dewormed_any[strata_matched_pos:stratum_end] | rep_array(1, curr_matched_stratum_size), name_match_false_pos, name_match_false_neg);
    
    strata_pos = stratum_end + 1;
  }
  
  // target += log_sum_exp(sum(name_matched_lp[1]), sum(name_matched_lp[2]));
  // 
  // target += bernoulli_lpmf(dewormed_any[name_matching_error_ids] | name_match_error_intercept + 
  //                                                                  name_match_error_diff * (monitored[name_matching_error_ids]));
}

generated quantities {
  // real<lower = 0, upper = 1> takeup_proportion[num_eval_treatment_prop] =
  //   takeup_proportion_rng(
  //     eval_treatment_prop_id, treatment_map_design_matrix, name_match_interact_map_design_matrix, dewormed_any, beta, stratum_intercept, cluster_effects,
  //     treatment_id, missing_treatment_id, treatment_sizes, missing_treatment_sizes, missing_stratum_id, missing_cluster_id
  //   );
}
