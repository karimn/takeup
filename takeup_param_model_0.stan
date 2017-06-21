functions {
  real[] takeup_proportion_rng(int[] eval_treatment_ids, // the treatment IDs to evaluate proportions for  
                               matrix treatment_map_design_matrix,
                               matrix census_covar_dm,
                               int[] dewormed_any,
                               vector[] beta,
                               vector[] census_covar_coef,
                               vector stratum_intercept,
                               vector cluster_effects,
                               int[] treatment_id,
                               int[] missing_treatment_id,
                               int[] treatment_sizes,
                               int[] missing_treatment_sizes,
                               int[,] missing_treatment_strata_sizes,
                               int[] missing_stratum_id_individ,
                               int[] missing_cluster_id) { 
                                 
    int num_eval_treatment_prop = size(eval_treatment_ids);
    real takeup_proportion[num_eval_treatment_prop]; 
    int num_strata = num_elements(stratum_intercept);
    
    for (treatment_index_index in 1:num_eval_treatment_prop) {
      int treatment_index = eval_treatment_ids[treatment_index_index]; 
      int curr_treatment_size = treatment_sizes[treatment_index];
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_index];
      int treatment_pos;  
      int missing_treatment_pos; 
      int treatment_end;
      int missing_treatment_end;
      int missing_treatment_stratum_pos;
      int missing_stratum_id[curr_missing_treatment_size];
      int treatment_missing_cluster_id[curr_missing_treatment_size];
      real missing_treatment_link[curr_missing_treatment_size];
      real missing_dewormed[curr_missing_treatment_size];
      int curr_missing_id[curr_missing_treatment_size];
      vector[curr_missing_treatment_size] missing_latent_utility = rep_vector(0, curr_missing_treatment_size);
      int missing_latent_util_pos = 1;
      
      if (treatment_index > 1) {
        treatment_pos = sum(treatment_sizes[1:(treatment_index - 1)]) + 1;
        missing_treatment_pos = sum(missing_treatment_sizes[1:(treatment_index - 1)]) + 1;
      } else {
        treatment_pos = 1;
        missing_treatment_pos = 1;
      }
      
      treatment_end = treatment_pos + curr_treatment_size - 1;
      missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      missing_treatment_stratum_pos = missing_treatment_pos;
      curr_missing_id = missing_treatment_id[missing_treatment_pos:missing_treatment_end];
      missing_stratum_id = missing_stratum_id_individ[missing_treatment_pos:missing_treatment_end]; // curr_missing_id];
      treatment_missing_cluster_id = missing_cluster_id[missing_treatment_pos:missing_treatment_end]; // curr_missing_id];
      // missing_treatment_link = missing_link_model[missing_stratum_id, treatment_index_index];
      
      for (stratum_index in 1:num_strata) {
        int curr_missing_treatment_stratum_size = missing_treatment_strata_sizes[treatment_index, stratum_index];
        int missing_treatment_stratum_end = missing_treatment_stratum_pos + curr_missing_treatment_stratum_size - 1;
        int missing_latent_util_end = missing_latent_util_pos + curr_missing_treatment_stratum_size - 1;
        
        missing_latent_utility[missing_latent_util_pos:missing_latent_util_end] = 
          stratum_intercept[missing_stratum_id[missing_latent_util_pos:missing_latent_util_end]] +
          cluster_effects[treatment_missing_cluster_id[missing_latent_util_pos:missing_latent_util_end]] +
          census_covar_dm[curr_missing_id[missing_latent_util_pos:missing_latent_util_end]] * census_covar_coef[stratum_index] +
          treatment_map_design_matrix[treatment_index] * beta[stratum_index]; 
          
        missing_treatment_stratum_pos = missing_treatment_stratum_end + 1;
        missing_latent_util_pos = missing_latent_util_end + 1;
      }
      
      for(missing_index in 1:curr_missing_treatment_size) {
        missing_dewormed[missing_index] = bernoulli_logit_rng(missing_latent_utility[missing_index]);
      }
      
      takeup_proportion[treatment_index_index] = 
        (sum(missing_dewormed) + sum(dewormed_any[treatment_id[treatment_pos:treatment_end]])) /
        (curr_missing_treatment_size + curr_treatment_size);
    }
    
    return(takeup_proportion);
  }

  real dewormed_matched_lpmf(int matched, int dewormed, real prob_false_pos, real prob_false_neg) {
    return(bernoulli_lpmf(matched | prob_false_pos + (1 - prob_false_pos - prob_false_neg) * dewormed));
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
  int<lower = 1, upper = num_obs> monitored_id[num_obs - num_name_matched];
  int<lower = 1, upper = num_obs> name_matched_strata_sizes[num_strata];
  int<lower = 1, upper = num_obs> name_matched_dewormed_strata_sizes[num_strata];
  int<lower = 1, upper = num_obs> name_matching_error_ids[num_name_matching_errors_ids];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  int<lower = 1, upper = num_obs> treatment_sizes[num_all_treatments]; // Number of observations per treatment ID, in ascending order
  int<lower = 1, upper = num_all_treatments> num_non_phone_owner_treatments;
  int<lower = 1, upper = num_all_treatments> non_phone_owner_treatments[num_non_phone_owner_treatments];
  int<lower = 1, upper = num_all_treatments> num_phone_owner_treatments;
  int<lower = 1, upper = num_all_treatments> phone_owner_treatments[num_phone_owner_treatments];
 
  // Counterfactual information for finite sample analysis
  
  int<lower = 1> num_missing; // # of missing counterfactuals
  int<lower = 1> num_eval_treatment_prop; // # of treatment cells to be imputed
  
  int<lower = 1, upper = num_obs> stratum_covar_id[num_obs]; // Observation indices ordered by stratum and endline covariate missingness
  int<lower = 1> stratum_missing_covar_sizes[num_strata]; // Number of covariate missingness per stratum (ordered by stratum)
  int<lower = 1, upper = num_all_treatments> eval_treatment_prop_id[num_eval_treatment_prop]; // Which treatments to impute
  int<lower = 1, upper = num_missing> missing_treatment_id[num_missing]; // Observation indices with missing treatment counterfactural, ordered by missing treatment ID
  int<lower = 1, upper = num_obs> missing_treatment_strata_sizes[num_all_treatments, num_strata];
  
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
  int<lower = 1, upper = num_strata> missing_stratum_id[num_missing] = stratum_id[missing_treatment_id];
  int<lower = 1, upper = num_clusters> missing_cluster_id[num_missing] = cluster_id[missing_treatment_id];
  // 
  // int<lower = 1, upper = num_obs> num_missing_endline_covar = sum(stratum_missing_covar_sizes);
  
  vector<lower = 0, upper = 1>[num_obs] monitored = 1 - name_matched;
  
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  matrix[num_obs, num_census_covar_coef] census_covar_dm = census_covar_map_dm[census_covar_id];
  
  // Constants for hyperpriors 
  real scale_df = 3;
  real scale_sigma = 1;
  
  real coef_df = 7;
  real coef_sigma = 10;
  
  int<lower = 1, upper = num_obs> missing_treatment_sizes[num_all_treatments]; 
  
  for (treatment_index in 1:num_all_treatments) {
    missing_treatment_sizes[treatment_index] = sum(missing_treatment_strata_sizes[treatment_index]);
  }
 
}

parameters {
  // Modelled parameters
  
  real<lower = 0, upper = 1> hyper_baseline_takeup; 
  vector[num_strata] stratum_intercept_raw;
  real<lower = 0> tau_stratum_intercept;
  
  vector[num_clusters] cluster_effects_raw;
  real<lower = 0> tau_cluster_effect;
  
  vector[num_all_treatment_coef] hyper_beta_raw;
  vector[num_all_treatment_coef] beta_raw[num_strata];
  // matrix[num_all_treatment_coef, num_strata] beta_raw;
  vector<lower = 0>[num_all_treatment_coef] tau_treatment;
  
  // vector[num_name_match_interact_coef] beta_name_match_interact;
  
  vector[num_census_covar_coef] hyper_census_covar_coef_raw;
  vector[num_census_covar_coef] census_covar_coef_raw[num_strata];
  // matrix[num_census_covar_coef, num_strata] census_covar_coef_raw;
  vector<lower = 0>[num_census_covar_coef] tau_census_covar;
  
  // matrix[num_census_covar_coef, num_experiment_coef] hyper_beta_census_covar_interact_raw; // census characteristics' effect on treatment effects
  // matrix[num_census_covar_coef, num_experiment_coef] beta_census_covar_interact[num_strata];
  // matrix<lower = 0>[num_census_covar_coef, num_experiment_coef] tau_beta_census_covar_interact[num_strata];
  
  /* Name matching model parameters *****************************************************/
  
  // real<lower = 0, upper = 1> name_match_false_pos;
  // real<lower = 0, upper = 1> name_match_false_neg;
  // 
  // real<lower = 0, upper = 1> name_match_error_intercept;
}

transformed parameters {
  real hyper_intercept = inv_logit(hyper_baseline_takeup);
  
  // real<lower = 0, upper = 1> name_match_error_diff =
  //   name_match_error_intercept * (name_match_false_neg - name_match_false_pos) / (1 - name_match_false_neg + name_match_false_pos);
    
  vector[num_all_treatment_coef] hyper_beta = hyper_beta_raw * coef_sigma; 
  vector[num_census_covar_coef] hyper_census_covar_coef = hyper_census_covar_coef_raw * coef_sigma;
  // matrix[num_census_covar_coef, num_experiment_coef] hyper_beta_census_covar_interact = hyper_beta_census_covar_interact_raw * coef_sigma;
  
  vector[num_strata] stratum_intercept = hyper_intercept + stratum_intercept_raw * tau_stratum_intercept;
  vector[num_clusters] cluster_effects = cluster_effects_raw * tau_cluster_effect;
  
  vector[num_all_treatment_coef] beta[num_strata];
  vector[num_census_covar_coef] census_covar_coef[num_strata];
  vector[num_obs] latent_utility = rep_vector(0, num_obs);
 
  {
    int stratum_pos = 1;

    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int curr_matched_stratum_size = name_matched_strata_sizes[strata_index];
      int curr_monitored_stratum_size = curr_stratum_size - curr_matched_stratum_size;
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      int monitored_stratum_end = stratum_pos + curr_monitored_stratum_size - 1;
      // matrix[curr_stratum_size, num_experiment_coef] census_covar_interact_dm;
      
      beta[strata_index] = hyper_beta + beta_raw[strata_index] .* tau_treatment;
      census_covar_coef[strata_index] = hyper_census_covar_coef + census_covar_coef_raw[strata_index] .* tau_census_covar;

      latent_utility[stratum_pos:stratum_end] =
          stratum_intercept[strata_index] +
          cluster_effects[cluster_id[stratum_pos:stratum_end]] +
          census_covar_dm[stratum_pos:stratum_end] * census_covar_coef[strata_index] +
          treatment_design_matrix[stratum_pos:stratum_end] * beta[strata_index];

      stratum_pos = stratum_end + 1;
    }
  }
}

model {
  
  // name_match_false_pos ~ beta(name_match_false_pos_alpha, name_match_false_pos_beta);
  
  stratum_intercept_raw ~ student_t(coef_df, 0, 1); 
  tau_stratum_intercept ~ student_t(scale_df, 0, scale_sigma);
  
  cluster_effects_raw ~ student_t(coef_df, 0, 1); 
  tau_cluster_effect ~ student_t(scale_df, 0, scale_sigma);
  
  hyper_beta_raw ~ student_t(coef_df, 0, 1); 
  beta_raw ~ multi_student_t(coef_df, rep_vector(0.0, num_all_treatment_coef), diag_matrix(rep_vector(1, num_all_treatment_coef)));
  // to_vector(beta_raw) ~ student_t(coef_df, 0, 1);
  tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  
  hyper_census_covar_coef_raw ~ student_t(coef_df, 0, 1);
  census_covar_coef_raw ~ multi_student_t(coef_df, rep_vector(0.0, num_census_covar_coef), diag_matrix(rep_vector(1, num_census_covar_coef)));
  // to_vector(census_covar_coef_raw) ~ student_t(coef_df, 0, 1);
  tau_census_covar ~ student_t(scale_df, 0, scale_sigma);
  
  dewormed_any ~ bernoulli_logit(latent_utility);
}

generated quantities {
  real<lower = 0, upper = 1> takeup_proportion[num_eval_treatment_prop] =
    takeup_proportion_rng(
      eval_treatment_prop_id,
      treatment_map_design_matrix,
      census_covar_dm,
      dewormed_any,
      beta, census_covar_coef, stratum_intercept, cluster_effects,
      treatment_id, missing_treatment_id,
      treatment_sizes, missing_treatment_sizes, missing_treatment_strata_sizes,
      missing_stratum_id, missing_cluster_id
    );
}