functions {
  vector treatment_cell_takeup_rng(int[] treatment_ids, 
                                   int[] missing_obs_ids, 
                                   int[] missing_cluster_id,
                                   int[] missing_treatment_sizes,
                                   int[] observed_treatment_sizes,
                                   matrix cluster_treat_latent_utility,
                                   vector census_covar_latent_utility,
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
        cluster_treat_latent_utility[treatment_ids_index, missing_cluster_id[missing_treatment_pos:missing_treatment_end]]' +
        census_covar_latent_utility[missing_obs_ids[missing_treatment_pos:missing_treatment_end]];
      
      vector[curr_missing_treatment_size] all_takeup;
      
      for (missing_obs_ids_index in 1:curr_missing_treatment_size) {
        all_takeup[missing_obs_ids_index] = bernoulli_logit_rng(missing_latent_utility[missing_obs_ids_index]);  
      }
     
      cell_takeup_prop[treatment_ids_index] = (sum(all_takeup) + observed_takeup_total[treatment_ids_index]) / (curr_missing_treatment_size + curr_observed_treatment_size);
        
      missing_treatment_pos = missing_treatment_end + 1;
      observed_treatment_pos = observed_treatment_end + 1;
    }
    
    return(cell_takeup_prop);
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
  int<lower = 1, upper = num_strata> stratum_id[num_obs]; // Observations' stratum IDs
  int<lower = 1, upper = num_strata> cluster_stratum_ids[num_clusters]; // Clusters' stratum IDs
  int<lower = 1, upper = num_clusters> strata_cluster_ids[num_clusters]; // Cluster IDs ordered by stratum
  int<lower = 1, upper = num_obs> cluster_obs_ids[num_obs]; // Observation IDs ordered by cluster
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; // Observations in strata 
  int<lower = 1, upper = num_clusters> strata_num_clusters[num_strata]; // Clusters in strata
  int<lower = 1, upper = num_obs> cluster_sizes[num_clusters]; // Observations in clusters
  
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
  
  real<lower = 0> lkj_df;
 
  // Binary deworming outcome 
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
  
  // Configuration
  
  int estimate_ate;
}

transformed data {
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
  matrix[num_obs, num_census_covar_coef] census_covar_dm = census_covar_map_dm[census_covar_id];
  matrix[num_all_treatments, num_all_treatments] diag_treatment_map_dm = diag_matrix(rep_vector(1, num_all_treatments));
}

parameters {
  // Modelled parameters
  real<lower = 0, upper = 1> hyper_baseline_takeup; // Completely uninformed prior

  vector[num_all_treatment_coef - 1] hyper_beta; // Non intercept hyperparameters
  vector[num_census_covar_coef] hyper_census_covar_coef;

  vector[num_all_treatment_coef] stratum_beta[num_strata];
  vector<lower = 0>[num_all_treatment_coef] stratum_tau_treatment;
  corr_matrix[num_all_treatment_coef] stratum_beta_corr_mat;

  vector[num_all_treatment_coef] cluster_beta[num_clusters];
  vector<lower = 0>[num_all_treatment_coef] cluster_tau_treatment;
  corr_matrix[num_all_treatment_coef] cluster_beta_corr_mat;

  vector[num_census_covar_coef] stratum_census_covar_coef[num_strata];
  vector<lower = 0>[num_census_covar_coef] stratum_tau_census_covar;
  corr_matrix[num_census_covar_coef] stratum_census_covar_corr_mat;
  
  vector[num_census_covar_coef] cluster_census_covar_coef[num_clusters];
  vector<lower = 0>[num_census_covar_coef] cluster_tau_census_covar;
  corr_matrix[num_census_covar_coef] cluster_census_covar_corr_mat;
}

transformed parameters {
  real hyper_intercept = logit(hyper_baseline_takeup);

  matrix[num_all_treatment_coef, num_clusters] cluster_beta_mat;
  vector[num_obs] census_covar_latent_utility = rep_vector(0, num_obs);
  vector[num_obs] latent_utility = rep_vector(0, num_obs);
  
  {
    int cluster_pos = 1;

    for (cluster_index in 1:num_clusters) {
      int curr_cluster_size = cluster_sizes[cluster_index];
      int cluster_end = cluster_pos + curr_cluster_size - 1;

      vector[num_all_treatment_coef] combined_cluster_beta = 
        stratum_beta[cluster_stratum_ids[cluster_index]] + cluster_beta[cluster_index];
        
      vector[num_census_covar_coef] combined_cluster_census_covar_coef = 
        stratum_census_covar_coef[cluster_stratum_ids[cluster_index]] + cluster_census_covar_coef[cluster_index];

      census_covar_latent_utility[cluster_obs_ids[cluster_pos:cluster_end]] =
        census_covar_dm[cluster_obs_ids[cluster_pos:cluster_end]] * combined_cluster_census_covar_coef;
        
      latent_utility[cluster_obs_ids[cluster_pos:cluster_end]] =
        census_covar_latent_utility[cluster_obs_ids[cluster_pos:cluster_end]] +
        census_covar_dm[cluster_obs_ids[cluster_pos:cluster_end]] * combined_cluster_census_covar_coef + 
        treatment_design_matrix[cluster_obs_ids[cluster_pos:cluster_end]] * combined_cluster_beta; 

      cluster_pos = cluster_end + 1;
      
      if (estimate_ate) {
        cluster_beta_mat[, cluster_index] = combined_cluster_beta;
      }
    }
  }
}

model {
  hyper_beta ~ student_t(coef_df, 0, coef_sigma); 
  hyper_census_covar_coef ~ student_t(coef_df, 0, coef_sigma);
  
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  stratum_beta_corr_mat ~ lkj_corr(lkj_df);
  
  cluster_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  cluster_beta_corr_mat ~ lkj_corr(lkj_df);
  
  stratum_tau_census_covar ~ student_t(scale_df, 0, scale_sigma);
  stratum_census_covar_corr_mat ~ lkj_corr(lkj_df);
  
  cluster_tau_census_covar ~ student_t(scale_df, 0, scale_sigma);
  cluster_census_covar_corr_mat ~ lkj_corr(lkj_df);

  {
    matrix[num_all_treatment_coef, num_all_treatment_coef] stratum_beta_vcov_mat = quad_form_diag(stratum_beta_corr_mat, stratum_tau_treatment);
    matrix[num_all_treatment_coef, num_all_treatment_coef] cluster_beta_vcov_mat = quad_form_diag(cluster_beta_corr_mat, cluster_tau_treatment);
    
    matrix[num_census_covar_coef, num_census_covar_coef] stratum_census_covar_vcov_mat = quad_form_diag(stratum_census_covar_corr_mat, stratum_tau_census_covar);
    matrix[num_census_covar_coef, num_census_covar_coef] cluster_census_covar_vcov_mat = quad_form_diag(cluster_census_covar_corr_mat, cluster_tau_census_covar);

    stratum_beta ~ multi_student_t(coef_df, append_row(hyper_intercept, hyper_beta), stratum_beta_vcov_mat);
    cluster_beta ~ multi_student_t(coef_df, rep_vector(0, num_all_treatment_coef), cluster_beta_vcov_mat);
    
    stratum_census_covar_coef ~ multi_student_t(coef_df, hyper_census_covar_coef, stratum_census_covar_vcov_mat);
    cluster_census_covar_coef ~ multi_student_t(coef_df, rep_vector(0, num_census_covar_coef), cluster_census_covar_vcov_mat);
  }
  
  dewormed_any ~ bernoulli_logit(latent_utility);
}

generated quantities {
  vector<lower = 0, upper = 1>[num_ate_treatments] treatment_cell_takeup;
  vector<lower = -1, upper = 1>[num_ate_pairs] treatment_cell_ate;
  
  if (estimate_ate) {
    matrix[num_ate_treatments, num_clusters] cluster_treat_latent_utility = treatment_map_design_matrix[ate_treatments] * cluster_beta_mat;
    
    treatment_cell_takeup = treatment_cell_takeup_rng(ate_treatments,
                                                      missing_obs_ids,
                                                      missing_treatment_cluster_id,
                                                      missing_treatment_sizes,
                                                      observed_treatment_sizes,
                                                      cluster_treat_latent_utility,
                                                      census_covar_latent_utility,
                                                      observed_takeup_total);

    treatment_cell_ate = treatment_cell_takeup[ate_pairs[, 1]] - treatment_cell_takeup[ate_pairs[, 2]];
  }
}
