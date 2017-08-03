functions {
  matrix treatment_cell_takeup_rng(int[] treatment_ids, 
                                   int[] missing_obs_ids, 
                                   int[] missing_stratum_id,
                                   int[] missing_cluster_id,
                                   matrix missing_census_covar_dm,
                                   vector missing_name_matched,
                                   vector observed_name_matched,
                                   matrix treatment_map_dm,
                                   matrix diag_treatment_map_dm,
                                   int[] missing_treatment_sizes,
                                   int[] observed_treatment_sizes,
                                   vector stratum_intercept,
                                   vector cluster_effects,
                                   vector census_covar_coef,
                                   vector stratum_name_matching_effect,
                                   matrix stratum_treatment_name_matching_interact,
                                   matrix stratum_treatment_coef,
                                   vector observed_dewormed_any) {
    int num_treatment_ids = size(treatment_ids);
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    matrix[num_treatment_ids, 2] cell_takeup_prop;
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
      
      vector[curr_missing_treatment_size] missing_latent_utility =
        stratum_intercept[missing_stratum_id[missing_treatment_pos:missing_treatment_end]] +
        cluster_effects[missing_cluster_id[missing_treatment_pos:missing_treatment_end]] +
        missing_census_covar_dm[missing_treatment_pos:missing_treatment_end] * census_covar_coef +
        stratum_treatment_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]] * treatment_map_dm[treatment_ids[treatment_ids_index]]' +
        missing_name_matched[missing_treatment_pos:missing_treatment_end] .* stratum_name_matching_effect[missing_stratum_id[missing_treatment_pos:missing_treatment_end]] +
        (diag_treatment_map_dm[treatment_ids[treatment_ids_index]] * stratum_treatment_name_matching_interact[, missing_stratum_id[missing_treatment_pos:missing_treatment_end]])' .*
          missing_name_matched[missing_treatment_pos:missing_treatment_end];
      
      vector[curr_missing_treatment_size] all_takeup;
      real all_observed_takeup_total = sum(observed_dewormed_any[observed_treatment_pos:observed_treatment_end]); 
      real monitored_observed_takeup_total = 
          sum(observed_dewormed_any[observed_treatment_pos:observed_treatment_end] .* (1 - observed_name_matched[observed_treatment_pos:observed_treatment_end]));
          
      real all_cell_num_obs = curr_missing_treatment_size + curr_observed_treatment_size;
      real monitored_cell_num_obs = all_cell_num_obs - 
        sum(missing_name_matched[missing_treatment_pos:missing_treatment_end]) - sum(observed_name_matched[missing_treatment_pos:missing_treatment_end]);
      
      for (missing_obs_ids_index in 1:curr_missing_treatment_size) {
        all_takeup[missing_obs_ids_index] = bernoulli_logit_rng(missing_latent_utility[missing_obs_ids_index]);  
      }
     
      cell_takeup_prop[treatment_ids_index, 1] = (sum(all_takeup) + all_observed_takeup_total) / all_cell_num_obs;
      cell_takeup_prop[treatment_ids_index, 2] = 
        (sum(all_takeup .* (1 - observed_name_matched[observed_treatment_pos:observed_treatment_end])) + monitored_observed_takeup_total) / monitored_cell_num_obs;
        
      missing_treatment_pos = missing_treatment_end + 1;
      observed_treatment_pos = observed_treatment_end + 1;
    }
    
    return(cell_takeup_prop);
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
  
  // int <lower = 0, upper = num_obs> num_name_matched;
  // int <lower = 0, upper = num_obs> num_name_matching_errors_ids;
  vector<lower = 0, upper = 1>[num_obs] name_matched;
  // int<lower = 1, upper = num_obs> name_matched_id[num_name_matched];
  // int<lower = 1, upper = num_obs> monitored_id[num_obs - num_name_matched];
  // int<lower = 1, upper = num_obs> name_matched_strata_sizes[num_strata];
  // int<lower = 1, upper = num_obs> name_matched_dewormed_strata_sizes[num_strata];
  // int<lower = 1, upper = num_obs> name_matching_error_ids[num_name_matching_errors_ids];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 
  int<lower = 1, upper = num_obs> cluster_sizes[num_clusters]; 
  int<lower = 1, upper = num_clusters> strata_num_clusters[num_strata];
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  int<lower = 1, upper = num_obs> treatment_sizes[num_all_treatments]; // Number of observations per treatment ID, in ascending order
 
  // Counterfactual information for finite sample analysis
  
  int<lower = 1, upper = num_obs> stratum_covar_id[num_obs]; // Observation indices ordered by stratum and endline covariate missingness
  int<lower = 1> stratum_missing_covar_sizes[num_strata]; // Number of covariate missingness per stratum (ordered by stratum)
  
  int<lower = 0, upper = num_all_treatments> num_non_phone_owner_treatments;
  int<lower = 0, upper = num_all_treatments> num_phone_owner_treatments;
  int<lower = 1, upper = num_all_treatments> non_phone_owner_treatments[num_non_phone_owner_treatments];
  int<lower = 1, upper = num_all_treatments> phone_owner_treatments[num_phone_owner_treatments];
  
  int<lower = 0> num_missing_non_phone_owner_obs_ids;
  int<lower = 0, upper = num_obs> missing_non_phone_owner_treatment_sizes[num_non_phone_owner_treatments];
  int<lower = 1, upper = num_obs> missing_non_phone_owner_obs_ids[num_missing_non_phone_owner_obs_ids];
  int<lower = 0> num_missing_phone_owner_obs_ids;
  int<lower = 0, upper = num_obs> missing_phone_owner_treatment_sizes[num_phone_owner_treatments];
  int<lower = 1, upper = num_obs> missing_phone_owner_obs_ids[num_missing_phone_owner_obs_ids];
  
  int<lower = 0> num_observed_non_phone_owner_obs_ids;
  int<lower = 0, upper = num_obs> observed_non_phone_owner_treatment_sizes[num_non_phone_owner_treatments];
  int<lower = 1, upper = num_obs> observed_non_phone_owner_obs_ids[num_observed_non_phone_owner_obs_ids];
  int<lower = 0> num_observed_phone_owner_obs_ids;
  int<lower = 0, upper = num_obs> observed_phone_owner_treatment_sizes[num_phone_owner_treatments];
  int<lower = 1, upper = num_obs> observed_phone_owner_obs_ids[num_observed_phone_owner_obs_ids];
  
  int<lower = 0, upper = num_all_treatments> num_non_phone_owner_ate_pairs;
  int<lower = 0, upper = num_all_treatments> num_phone_owner_ate_pairs;
  
  int<lower = 1, upper = num_all_treatments> non_phone_owner_ate_pairs[num_non_phone_owner_ate_pairs, 2];
  int<lower = 1, upper = num_all_treatments> phone_owner_ate_pairs[num_phone_owner_ate_pairs, 2];
  
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
  // int<lower = 1, upper = num_obs> num_missing_endline_covar = sum(stratum_missing_covar_sizes);
  
  // vector<lower = 0, upper = 1>[num_obs] monitored = 1 - name_matched;
  
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  matrix[num_obs, num_census_covar_coef] census_covar_dm = census_covar_map_dm[census_covar_id];
  matrix[num_all_treatments, num_all_treatment_coef] diag_treatment_map_dm;
  matrix[num_obs, num_all_treatment_coef] diag_treatment_dm;
  
  int non_phone_missing_treatment_stratum_id[num_missing_non_phone_owner_obs_ids] = stratum_id[missing_non_phone_owner_obs_ids];
  int phone_missing_treatment_stratum_id[num_missing_phone_owner_obs_ids] = stratum_id[missing_phone_owner_obs_ids];
  int non_phone_missing_treatment_cluster_id[num_missing_non_phone_owner_obs_ids] = cluster_id[missing_non_phone_owner_obs_ids];
  int phone_missing_treatment_cluster_id[num_missing_phone_owner_obs_ids] = cluster_id[missing_phone_owner_obs_ids];
  vector<lower = 0, upper = 1>[num_missing_non_phone_owner_obs_ids] non_phone_missing_name_matched = name_matched[missing_non_phone_owner_obs_ids];
  vector<lower = 0, upper = 1>[num_missing_phone_owner_obs_ids] phone_missing_name_matched = name_matched[missing_phone_owner_obs_ids];
  vector<lower = 0, upper = 1>[num_observed_non_phone_owner_obs_ids] non_phone_observed_name_matched = name_matched[observed_non_phone_owner_obs_ids];
  vector<lower = 0, upper = 1>[num_observed_phone_owner_obs_ids] phone_observed_name_matched = name_matched[observed_phone_owner_obs_ids];
  
  matrix[num_missing_non_phone_owner_obs_ids, num_census_covar_coef] non_phone_missing_census_covar_dm = census_covar_dm[missing_non_phone_owner_obs_ids];
  matrix[num_missing_phone_owner_obs_ids, num_census_covar_coef] phone_missing_census_covar_dm = census_covar_dm[missing_phone_owner_obs_ids];
  
  int observed_non_phone_dewormed_any[num_observed_non_phone_owner_obs_ids] = dewormed_any[observed_non_phone_owner_obs_ids];
  int observed_phone_dewormed_any[num_observed_phone_owner_obs_ids] = dewormed_any[observed_phone_owner_obs_ids];
  
  // Constants for hyperpriors 
  real scale_df = 3;
  real scale_sigma = 1;
  
  real coef_df = 7;
  real coef_sigma = 10;

  diag_treatment_map_dm[1] = rep_row_vector(0, num_all_treatment_coef);
  diag_treatment_map_dm[2:num_all_treatments] = diag_matrix(rep_vector(1, num_all_treatment_coef));
  
  diag_treatment_dm = diag_treatment_map_dm[obs_treatment];
}

parameters {
  // Modelled parameters
  
  real<lower = 0, upper = 1> hyper_baseline_takeup; // Completely uninformed prior
  vector[num_strata] stratum_intercept_raw;
  real<lower = 0> tau_stratum_intercept;
  
  vector[num_clusters] cluster_effects_raw;
  real<lower = 0> tau_cluster_effect;
  
  vector[num_all_treatment_coef] hyper_beta_raw;
  vector[num_all_treatment_coef] stratum_beta_raw[num_strata];
  vector<lower = 0>[num_all_treatment_coef] stratum_tau_treatment;
  
  vector[num_strata] stratum_name_matching_effect_raw;
  real<lower = 0> tau_stratum_name_matching;
  
  matrix[num_all_treatment_coef, num_strata] stratum_treatment_name_matching_interact_raw;
  real<lower = 0> tau_stratum_treatment_name_matching_interact;
  
  // vector[num_name_match_interact_coef] beta_name_match_interact;
  
  vector[num_census_covar_coef] hyper_census_covar_coef_raw;
  // vector[num_census_covar_coef] census_covar_coef_raw[num_strata];
  // matrix[num_census_covar_coef, num_strata] census_covar_coef_raw;
  // vector<lower = 0>[num_census_covar_coef] tau_census_covar;
  
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
  real hyper_intercept = logit(hyper_baseline_takeup);
  
  // real<lower = 0, upper = 1> name_match_error_diff =
  //   name_match_error_intercept * (name_match_false_neg - name_match_false_pos) / (1 - name_match_false_neg + name_match_false_pos);
    
  vector[num_all_treatment_coef] hyper_beta = hyper_beta_raw * coef_sigma; 
  vector[num_census_covar_coef] hyper_census_covar_coef = hyper_census_covar_coef_raw * coef_sigma;
  // matrix[num_census_covar_coef, num_experiment_coef] hyper_beta_census_covar_interact = hyper_beta_census_covar_interact_raw * coef_sigma;
  
  vector[num_strata] stratum_intercept = hyper_intercept + stratum_intercept_raw * tau_stratum_intercept;
  vector[num_clusters] cluster_effects = cluster_effects_raw * tau_cluster_effect;
  
  vector[num_all_treatment_coef] stratum_beta[num_strata];
  vector[num_obs] latent_utility = rep_vector(0, num_obs);
  
  vector[num_strata] stratum_name_matching_effect = stratum_name_matching_effect_raw * tau_stratum_name_matching;
  matrix[num_all_treatment_coef, num_strata] stratum_treatment_name_matching_interact = 
    stratum_treatment_name_matching_interact_raw * tau_stratum_treatment_name_matching_interact;
 
  {
    int stratum_pos = 1;

    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      // int curr_matched_stratum_size = name_matched_strata_sizes[strata_index];
      // int curr_monitored_stratum_size = curr_stratum_size - curr_matched_stratum_size;
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      // int monitored_stratum_end = stratum_pos + curr_monitored_stratum_size - 1;
      // matrix[curr_stratum_size, num_experiment_coef] census_covar_interact_dm;
      
      stratum_beta[strata_index] = hyper_beta + stratum_beta_raw[strata_index] .* stratum_tau_treatment;
              
      // census_covar_coef[strata_index] = hyper_census_covar_coef + census_covar_coef_raw[strata_index] .* tau_census_covar;
      
      latent_utility[stratum_pos:stratum_end] =
          stratum_intercept[strata_index] +
          cluster_effects[cluster_id[stratum_pos:stratum_end]] +
          census_covar_dm[stratum_pos:stratum_end] * hyper_census_covar_coef +
          treatment_design_matrix[stratum_pos:stratum_end] * stratum_beta[strata_index] +
          name_matched[stratum_pos:stratum_end] * stratum_name_matching_effect[strata_index] +
          (diag_treatment_dm[stratum_pos:stratum_end] * stratum_treatment_name_matching_interact[, strata_index]) .* name_matched[stratum_pos:stratum_end];
      
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
  stratum_beta_raw ~ multi_student_t(coef_df, rep_vector(0.0, num_all_treatment_coef), diag_matrix(rep_vector(1, num_all_treatment_coef)));
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  
  stratum_name_matching_effect_raw ~ student_t(coef_df, 0, 1);
  tau_stratum_name_matching ~ student_t(scale_df, 0, scale_sigma);
  to_vector(stratum_treatment_name_matching_interact_raw) ~ student_t(coef_df, 0, 1);
  tau_stratum_treatment_name_matching_interact ~ student_t(scale_df, 0, scale_sigma);
  
  hyper_census_covar_coef_raw ~ student_t(coef_df, 0, 1);
  // census_covar_coef_raw ~ multi_student_t(coef_df, rep_vector(0.0, num_census_covar_coef), diag_matrix(rep_vector(1, num_census_covar_coef)));
  // tau_census_covar ~ student_t(scale_df, 0, scale_sigma);
  
  dewormed_any ~ bernoulli_logit(latent_utility);
}

generated quantities {
  matrix<lower = 0, upper = 1>[num_non_phone_owner_treatments, 2] non_phone_takeup; 
  matrix<lower = 0, upper = 1>[num_phone_owner_treatments, 2] phone_takeup; 
  vector<lower = -1, upper = 1>[num_non_phone_owner_ate_pairs] non_phone_takeup_ate;
  vector<lower = -1, upper = 1>[num_phone_owner_ate_pairs] phone_takeup_ate;
  vector[num_non_phone_owner_ate_pairs] non_phone_takeup_ate_percent;
  vector[num_phone_owner_ate_pairs] phone_takeup_ate_percent;
  
  matrix[num_strata, num_all_treatment_coef] stratum_beta_mat;
  
  for (stratum_index in 1:num_strata) {
    stratum_beta_mat[stratum_index, ] = stratum_beta[stratum_index]';
  }
  
  non_phone_takeup = treatment_cell_takeup_rng(non_phone_owner_treatments,
                                           missing_non_phone_owner_obs_ids,
                                           non_phone_missing_treatment_stratum_id,
                                           non_phone_missing_treatment_cluster_id,
                                           non_phone_missing_census_covar_dm,
                                           non_phone_missing_name_matched,
                                           non_phone_observed_name_matched,
                                           treatment_map_design_matrix,
                                           diag_treatment_dm,
                                           missing_non_phone_owner_treatment_sizes,
                                           observed_non_phone_owner_treatment_sizes,
                                           stratum_intercept,
                                           cluster_effects,
                                           hyper_census_covar_coef,
                                           stratum_name_matching_effect,
                                           stratum_treatment_name_matching_interact_raw,
                                           stratum_beta_mat,
                                           to_vector(observed_non_phone_dewormed_any));
                                           
  phone_takeup = treatment_cell_takeup_rng(phone_owner_treatments,
                                           missing_phone_owner_obs_ids,
                                           phone_missing_treatment_stratum_id,
                                           phone_missing_treatment_cluster_id,
                                           phone_missing_census_covar_dm,
                                           phone_missing_name_matched,
                                           phone_observed_name_matched,
                                           treatment_map_design_matrix,
                                           diag_treatment_dm,
                                           missing_phone_owner_treatment_sizes,
                                           observed_phone_owner_treatment_sizes,
                                           stratum_intercept,
                                           cluster_effects,
                                           hyper_census_covar_coef,
                                           stratum_name_matching_effect,
                                           stratum_treatment_name_matching_interact_raw,
                                           stratum_beta_mat,
                                           to_vector(observed_phone_dewormed_any));
                                           
  non_phone_takeup_ate = non_phone_takeup[non_phone_owner_ate_pairs[, 1], 1] - non_phone_takeup[non_phone_owner_ate_pairs[, 2], 1];
  phone_takeup_ate = phone_takeup[phone_owner_ate_pairs[, 1], 1] - phone_takeup[phone_owner_ate_pairs[, 2], 1];
  non_phone_takeup_ate_percent = non_phone_takeup_ate ./ non_phone_takeup[non_phone_owner_ate_pairs[, 2], 2];
  phone_takeup_ate_percent = phone_takeup_ate ./ phone_takeup[phone_owner_ate_pairs[, 2], 2];
}
