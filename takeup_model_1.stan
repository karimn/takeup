functions {
  vector treatment_cell_takeup_rng(int[] treatment_ids, 
                                   int[] missing_obs_ids, 
                                   int[] observed_obs_ids, 
                                   int[] missing_stratum_id,
                                   int[] missing_cluster_id,
                                   int[] private_value_calendar_coef,
                                   int[] private_value_bracelet_coef,
                                   matrix missing_census_covar_dm,
                                   matrix treatment_map_dm,
                                   int[] missing_treatment_sizes,
                                   int[] observed_treatment_sizes,
                                   vector stratum_intercept,
                                   vector cluster_effects,
                                   vector census_covar_coef,
                                   matrix stratum_treatment_coef,
                                   int[] observed_takeup_total) { 
    int num_treatment_ids = size(treatment_ids);
    vector[num_treatment_ids] cell_takeup_prop;
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
      
      vector[curr_missing_treatment_size] missing_latent_utility =
        stratum_intercept[missing_stratum_id[missing_treatment_pos:missing_treatment_end]] +
        cluster_effects[missing_cluster_id[missing_treatment_pos:missing_treatment_end]] +
        missing_census_covar_dm[missing_treatment_pos:missing_treatment_end] * census_covar_coef +
        stratum_treatment_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end]] * treatment_map_dm[treatment_ids[treatment_ids_index]]';
      
      vector[curr_missing_treatment_size] all_takeup;
      // real all_observed_takeup_total = sum(observed_dewormed_any[observed_obs_ids[observed_treatment_pos:observed_treatment_end]]); 
      // real all_cell_num_obs = curr_missing_treatment_size + curr_observed_treatment_size;
      
      for (missing_obs_ids_index in 1:curr_missing_treatment_size) {
        all_takeup[missing_obs_ids_index] = bernoulli_logit_rng(missing_latent_utility[missing_obs_ids_index]);  
      }
     
      cell_takeup_prop[treatment_ids_index] = (sum(all_takeup) + observed_takeup_total[treatment_ids_index]) / (curr_missing_treatment_size + curr_observed_treatment_size);
        
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
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  int num_private_value_calendar_coef;
  int num_private_value_bracelet_coef;
  int private_value_calendar_coef[num_private_value_calendar_coef];
  int private_value_bracelet_coef[num_private_value_bracelet_coef];
  int private_value_bracelet_indicator;
  int not_private_value_bracelet_coef[num_all_treatment_coef - num_private_value_bracelet_coef];
 
  // Below references to "maps" refer to a finite set of treatment cells that has a one-to-many relationship to actual observations
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix; // Design matrix generated from treatment map
  int<lower = 1, upper = num_all_treatments> obs_treatment[num_obs]; // ID of observed treatment (from treatment map)
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  
  // Covariates
  
  int<lower = 1> num_census_covar_coef;
  int<lower = 1> num_endline_covar_coef; 
  int<lower = 1> num_distinct_census_covar;
  
  matrix[num_distinct_census_covar, num_census_covar_coef] census_covar_map_dm; 
  int census_covar_id[num_obs];
  
  // Counterfactual information for finite sample analysis
  
  int<lower = 0, upper = num_all_treatments> num_ate_treatments;
  int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments];
  int<lower = 0, upper = num_obs> observed_takeup_total[num_ate_treatments];
  
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
  
  // Hyperprior settings 
  
  real<lower = 0> scale_df;
  real<lower = 0> scale_sigma;
  
  real<lower = 0> coef_df;
  real<lower = 0> coef_sigma;
 
  // Binary deworming outcome 
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
  
  // Configuration
  
  int estimate_ate;
}

transformed data {
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  matrix[num_obs, num_census_covar_coef] census_covar_dm = census_covar_map_dm[census_covar_id];
  matrix[num_all_treatments, num_all_treatments] diag_treatment_map_dm = diag_matrix(rep_vector(1, num_all_treatments));
  
  int num_not_private_value_bracelet_coef = num_all_treatment_coef - num_private_value_bracelet_coef;
}

parameters {
  // Modelled parameters
  
  real<lower = 0, upper = 1> hyper_baseline_takeup; // Completely uninformed prior
  vector[num_strata] stratum_intercept;
  real<lower = 0> tau_stratum_intercept;
  
  vector[num_clusters] cluster_effects;
  real<lower = 0> tau_cluster_effect;
  
  vector[num_not_private_value_bracelet_coef] hyper_beta; // No tau for hyper parameter; coef_sigma is the SD
  vector[num_not_private_value_bracelet_coef] stratum_beta[num_strata];
  vector<lower = 0>[num_not_private_value_bracelet_coef] stratum_tau_treatment;
  
  vector[num_census_covar_coef] hyper_census_covar_coef;
}

transformed parameters {
  real hyper_intercept = logit(hyper_baseline_takeup);
  matrix[num_strata, num_all_treatment_coef] stratum_beta_mat;
  vector[num_obs] latent_utility = rep_vector(0, num_obs);
  
  {
    int stratum_pos = 1;
    int bracelet_val_stratum_pos = 1;

    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      vector[num_all_treatment_coef] local_stratum_beta = rep_vector(0, num_all_treatment_coef);
      
      local_stratum_beta[not_private_value_bracelet_coef] = stratum_beta[strata_index];
      local_stratum_beta[private_value_bracelet_coef] = local_stratum_beta[private_value_calendar_coef];
              
      latent_utility[stratum_pos:stratum_end] =
          stratum_intercept[strata_index] +
          cluster_effects[cluster_id[stratum_pos:stratum_end]] +
          census_covar_dm[stratum_pos:stratum_end] * hyper_census_covar_coef +
          treatment_design_matrix[stratum_pos:stratum_end] * local_stratum_beta; 
          
      stratum_pos = stratum_end + 1;
      stratum_beta_mat[strata_index] = local_stratum_beta';
    }
  }
}

model {
  tau_stratum_intercept ~ student_t(scale_df, 0, scale_sigma);
  stratum_intercept ~ student_t(coef_df, hyper_intercept, tau_stratum_intercept); 
  
  tau_cluster_effect ~ student_t(scale_df, 0, scale_sigma);
  cluster_effects ~ student_t(coef_df, 0, tau_cluster_effect); 
  
  hyper_beta ~ student_t(coef_df, 0, coef_sigma); 
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  stratum_beta ~ multi_student_t(coef_df, hyper_beta, diag_matrix(stratum_tau_treatment));
  
  hyper_census_covar_coef ~ student_t(coef_df, 0, coef_sigma);
  
  dewormed_any ~ bernoulli_logit(latent_utility);
}

generated quantities {
  vector<lower = 0, upper = 1>[num_ate_treatments] treatment_cell_takeup;
  vector<lower = -1, upper = 1>[num_ate_pairs] treatment_cell_ate;
  
  if (estimate_ate) {
    treatment_cell_takeup = treatment_cell_takeup_rng(ate_treatments,
                                                      missing_obs_ids,
                                                      observed_obs_ids,
                                                      missing_treatment_stratum_id,
                                                      missing_treatment_cluster_id,
                                                      private_value_calendar_coef,
                                                      private_value_bracelet_coef,
                                                      missing_census_covar_dm,
                                                      treatment_map_design_matrix,
                                                      missing_treatment_sizes,
                                                      observed_treatment_sizes,
                                                      stratum_intercept,
                                                      cluster_effects,
                                                      hyper_census_covar_coef,
                                                      stratum_beta_mat,
                                                      observed_takeup_total);
                                                      
    treatment_cell_ate = treatment_cell_takeup[ate_pairs[, 1]] - treatment_cell_takeup[ate_pairs[, 2]]; 
  }
}
