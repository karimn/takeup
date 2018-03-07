functions {
  vector treatment_deworming_rng(vector missing_latent_var) {
    int curr_treatment_size = rows(missing_latent_var); 
    vector[curr_treatment_size] deworming_mask = rep_vector(0, curr_treatment_size); 
    
    for (obs_ids_index in 1:curr_treatment_size) {
      deworming_mask[obs_ids_index] = bernoulli_logit_rng(missing_latent_var[obs_ids_index]);
    }
  
    return deworming_mask; 
  }
  
  vector treatment_cell_deworming_prop_rng(int num_treatment_ids,
                                           int[] missing_obs_ids, 
                                           int[] missing_level_id, // Either stratum or cluster
                                           int[] missing_treatment_sizes,
                                           int[] observed_treatment_sizes,
                                           matrix latent_var_map, 
                                           int[] observed_dewormed_any) {
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    vector[num_treatment_ids] cell_deworming_takeup = rep_vector(0, num_treatment_ids);
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      
      vector[curr_missing_treatment_size] missing_deworming_mask = 
        treatment_deworming_rng(latent_var_map[missing_level_id[missing_treatment_pos:missing_treatment_end], treatment_ids_index]);
                                    
      if (curr_observed_treatment_size > 0) {
        int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
          
        cell_deworming_takeup[treatment_ids_index] = 
          (sum(missing_deworming_mask) + sum(observed_dewormed_any[observed_treatment_pos:observed_treatment_end]))  
          / (curr_missing_treatment_size + curr_observed_treatment_size);
          
        observed_treatment_pos = observed_treatment_end + 1;
      } else {
        cell_deworming_takeup[treatment_ids_index] = mean(missing_deworming_mask);
      }
        
      missing_treatment_pos = missing_treatment_end + 1;
    }
    
    return cell_deworming_takeup;
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
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix;
  matrix[num_obs, num_all_treatment_coef] Q_treatment_design_matrix;
  matrix[num_all_treatment_coef, num_all_treatment_coef] R_inv_treatment_design_matrix;
  int<lower = 1, upper = num_all_treatments> obs_treatment[num_obs]; // ID of observed treatment (from treatment map)
  
  int<lower = 1, upper = num_all_treatments> num_cluster_level_treatments;
  int<lower = 1> within_cluster_treatment_sizes[num_cluster_level_treatments];
  int<lower = 0, upper = num_cluster_level_treatments> unique_within_cluster_treatment_sizes[2]; // Control and treated clusters have different sizes
  matrix<lower = 0, upper = 1>[sum(within_cluster_treatment_sizes), num_all_treatment_coef] within_cluster_treatment_map;
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs]; // Observations' stratum IDs
  int<lower = 1, upper = num_strata> cluster_stratum_ids[num_clusters]; // Clusters' stratum IDs
  int<lower = 1, upper = num_obs> cluster_obs_ids[num_obs]; // Observation IDs ordered by cluster
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; // Observations in strata 
  int<lower = 1, upper = num_clusters> strata_num_clusters[num_strata]; // Clusters in strata
  int<lower = 1, upper = num_obs> cluster_sizes[num_clusters]; // Observations in clusters
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  
  // Deworming outcomes
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
 
  int num_dewormed; 
  int strata_dewormed_sizes[num_strata];
  int dewormed_ids[num_dewormed];
  int<lower = 1> stratum_dewormed_index[num_dewormed];
  
  // Covariates
  
  int<lower = 1> num_census_covar_coef;
  matrix[num_obs, num_census_covar_coef] census_covar_dm;
 
  // Counterfactuals and ATE
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_treatments;
  int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments]; 
  
  int<lower = 0> num_missing_obs_ids;
  int<lower = 0, upper = num_obs> missing_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> missing_obs_ids[num_missing_obs_ids];
  
  int<lower = 0> num_observed_obs_ids;
  int<lower = 0, upper = num_obs> observed_treatment_sizes[num_ate_treatments];
  int<lower = 1, upper = num_obs> observed_obs_ids[num_observed_obs_ids];
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_pairs;
  
  int<lower = 1, upper = num_all_treatments * num_all_treatments> ate_pairs[num_ate_pairs, 2];
  
  int missing_treatment_stratum_id[num_missing_obs_ids]; 
  int missing_treatment_cluster_id[num_missing_obs_ids]; 
  
  // matrix[num_missing_obs_ids, num_census_covar_coef] missing_census_covar_dm;
  
  int observed_dewormed_any[num_observed_obs_ids];
  
  // Constants for hyperpriors 
  
  real<lower = 0> scale_sigma;
  real<lower = 0> cluster_scale_sigma;
  
  real<lower = 0> hyper_coef_sigma;
  real<lower = 0> hyper_intercept_sigma;
  
  real<lower = 0> lkj_df;
  
  // Configuration
  
  int<lower = 0, upper = 1> estimate_ate;
  
  // This allows us to specify how many levels the model has:
  // 1: Only the toplevel parameters
  // 2: Strata level (counties)
  // 3: Cluster level (villages)
  int<lower = 1, upper = 3> model_levels;
}

transformed data {
  int<lower = 0> strata_num_all_treatment_coef = model_levels > 1 ? num_all_treatment_coef : 0;
  int<lower = 0> cluster_num_all_treatment_coef = model_levels > 2 ? num_all_treatment_coef : 0;
  
  int<lower = 0, upper = num_all_treatments> num_within_cluster_rows = model_levels > 2 ? sum(within_cluster_treatment_sizes) : 0;
  int<lower = 0, upper = num_cluster_level_treatments> num_control_cluster_treatments = 0;
  int<lower = 0, upper = num_cluster_level_treatments> num_treated_cluster_treatments = num_cluster_level_treatments;
 
  if (model_levels > 2) { 
    for (cluster_level_treat_index in 1:num_cluster_level_treatments) {
      if (within_cluster_treatment_sizes[cluster_level_treat_index] == unique_within_cluster_treatment_sizes[1]) {
        num_control_cluster_treatments += 1;
        num_treated_cluster_treatments -= 1;
      }
    }
  } else {
    num_treated_cluster_treatments = 0;
  }
}

parameters {
  vector[num_all_treatment_coef] hyper_beta;
  
  real<lower = 2, upper = 8> stratum_student_df;
  real<lower = 2, upper = 8> cluster_student_df;
  
  matrix[strata_num_all_treatment_coef, num_strata] strata_beta_raw;
  vector<lower = 0>[strata_num_all_treatment_coef] strata_beta_tau;
  cholesky_factor_corr[strata_num_all_treatment_coef] strata_beta_L_corr_mat;
  
  matrix[cluster_num_all_treatment_coef, num_clusters] cluster_beta;
  matrix<lower = 0>[num_within_cluster_rows, num_strata] cluster_beta_tau;
  cholesky_factor_corr[unique_within_cluster_treatment_sizes[1]] control_cluster_L_corr_mat[num_strata, num_control_cluster_treatments];
  cholesky_factor_corr[unique_within_cluster_treatment_sizes[2]] treated_cluster_L_corr_mat[num_strata, num_treated_cluster_treatments];
}

transformed parameters {
  // These are always specified even if the model is not multilevel; they'll just have the same values one level up.
  matrix[num_all_treatment_coef, num_strata] strata_beta; 
  matrix[num_all_treatment_coef, num_clusters] effective_cluster_beta;
  matrix[num_within_cluster_rows, num_clusters] within_cluster_beta_raw; 
  
  {
    int cluster_pos = 1;
    matrix[strata_num_all_treatment_coef, strata_num_all_treatment_coef] strata_beta_L_vcov = diag_pre_multiply(strata_beta_tau, strata_beta_L_corr_mat);
    
    if (model_levels > 1) {
      strata_beta = rep_matrix(hyper_beta, num_strata) + strata_beta_L_vcov * strata_beta_raw;
      
      if (model_levels > 2) {
        effective_cluster_beta = cluster_beta;
      }  
    } else {
      strata_beta = rep_matrix(hyper_beta, num_strata);
    }

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      if (model_levels > 2) {
        matrix[num_within_cluster_rows, curr_num_clusters] within_cluster_beta;
        
        matrix[num_within_cluster_rows, num_within_cluster_rows] within_cluster_beta_L_corr_mat; // = rep_matrix(0, num_within_cluster_rows, num_within_cluster_rows);
        matrix[num_within_cluster_rows, num_within_cluster_rows] within_cluster_beta_L_vcov_mat;  
        vector[num_within_cluster_rows] stratum_cluster_level_mu;
        
        int control_pos = 1;
        int treated_pos = 1;
        int corr_pos = 1;
        
        within_cluster_beta = within_cluster_treatment_map * cluster_beta[, cluster_pos:cluster_end];
        stratum_cluster_level_mu = within_cluster_treatment_map * strata_beta[, stratum_index];
      
        for (cluster_level_treat_index in 1:num_cluster_level_treatments) {
          int corr_end; 
          
          if (within_cluster_treatment_sizes[cluster_level_treat_index] == unique_within_cluster_treatment_sizes[1]) {
            corr_end = corr_pos + unique_within_cluster_treatment_sizes[1] - 1;
            within_cluster_beta_L_corr_mat[corr_pos:corr_end, corr_pos:corr_end] = control_cluster_L_corr_mat[stratum_index, control_pos];
            control_pos += 1; 
          } else {
            corr_end = corr_pos + unique_within_cluster_treatment_sizes[2] - 1;
            within_cluster_beta_L_corr_mat[corr_pos:corr_end, corr_pos:corr_end] = treated_cluster_L_corr_mat[stratum_index, treated_pos];
            treated_pos += 1; 
          }
          
          corr_pos = corr_end + 1;
        }
      
        within_cluster_beta_L_vcov_mat = diag_pre_multiply(cluster_beta_tau[, stratum_index], within_cluster_beta_L_corr_mat);
        
        within_cluster_beta_raw[, cluster_pos:cluster_end] = 
          within_cluster_beta_L_vcov_mat \ (within_cluster_beta - rep_matrix(stratum_cluster_level_mu, curr_num_clusters));
      } else {
        effective_cluster_beta[, cluster_pos:cluster_end] = rep_matrix(strata_beta[, stratum_index], curr_num_clusters);
      }
      
      cluster_pos = cluster_end + 1;
    }
  }
}

model {
  hyper_beta[1] ~ normal(0, hyper_intercept_sigma);
  hyper_beta[2:num_all_treatment_coef] ~ normal(0, hyper_coef_sigma);
 
  if (model_levels > 1) { 
    to_vector(strata_beta_raw) ~ student_t(stratum_student_df, 0, 1);
    strata_beta_tau ~ normal(0, scale_sigma);
    strata_beta_L_corr_mat ~ lkj_corr_cholesky(lkj_df);
  }
 
  if (model_levels > 2) { 
    to_vector(within_cluster_beta_raw) ~ student_t(cluster_student_df, 0, 1);
    to_vector(cluster_beta_tau) ~ normal(0, cluster_scale_sigma);
  }
  
  {
    vector[num_obs] latent_var = rep_vector(0, num_obs);
   
    int obs_pos = 1; 
    int cluster_pos = 1;

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      int control_pos = 1;
      int treated_pos = 1;
      
      if (model_levels > 2) {
        for (cluster_level_treat_index in 1:num_cluster_level_treatments) {
          if (within_cluster_treatment_sizes[cluster_level_treat_index] == unique_within_cluster_treatment_sizes[1]) {
            control_cluster_L_corr_mat[stratum_index, control_pos] ~ lkj_corr_cholesky(lkj_df);
            control_pos += 1; 
          } else {
            treated_cluster_L_corr_mat[stratum_index, treated_pos] ~ lkj_corr_cholesky(lkj_df);
            treated_pos += 1; 
          }
        }
      }

      for (cluster_pos_index in cluster_pos:cluster_end) {
        int curr_cluster_num_obs = cluster_sizes[cluster_pos_index];
        int obs_end = obs_pos + curr_cluster_num_obs - 1;
          
        latent_var[obs_pos:obs_end] = treatment_design_matrix[obs_pos:obs_end] * effective_cluster_beta[, cluster_pos_index];
                           
        obs_pos = obs_end + 1;
      }

      cluster_pos = cluster_end + 1;
    }
    
    dewormed_any ~ bernoulli_logit(latent_var);
  }
}

generated quantities { 
  vector<lower = 0, upper = 1>[num_ate_treatments] est_takeup = rep_vector(0, num_ate_treatments);
  vector<lower = 0, upper = 1>[num_ate_treatments] sp_est_takeup = rep_vector(0, num_ate_treatments);
  matrix<lower = 0, upper = 1>[num_strata, num_ate_treatments] stratum_sp_est_takeup = rep_matrix(0, num_strata, num_ate_treatments);
  matrix<lower = 0, upper = 1>[num_clusters, num_ate_treatments] cluster_sp_est_takeup = rep_matrix(0, num_clusters, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] est_takeup_ate = rep_vector(0, num_ate_pairs);
  vector<lower = -1, upper = 1>[num_ate_pairs] sp_est_takeup_ate = rep_vector(0, num_ate_pairs);
  matrix<lower = -1, upper = 1>[num_strata, num_ate_pairs] stratum_sp_est_takeup_ate = rep_matrix(0, num_strata, num_ate_pairs);
  matrix<lower = -1, upper = 1>[num_clusters, num_ate_pairs] cluster_sp_est_takeup_ate = rep_matrix(0, num_clusters, num_ate_pairs);
  
  vector[num_ate_treatments] hyper_latent_var_map; 
  matrix[num_strata, num_ate_treatments] stratum_latent_var_map; 

  if (estimate_ate) {
    int stratum_pos = 1;
    matrix[num_clusters, num_ate_treatments] cluster_latent_var_map; 
  
    hyper_latent_var_map = treatment_map_design_matrix[ate_treatments] * hyper_beta;
    stratum_latent_var_map = (treatment_map_design_matrix[ate_treatments] * strata_beta)';
    cluster_latent_var_map = (treatment_map_design_matrix[ate_treatments] * effective_cluster_beta)';

    est_takeup =
      treatment_cell_deworming_prop_rng(num_ate_treatments,
                                        missing_obs_ids,
                                        missing_treatment_cluster_id,
                                        missing_treatment_sizes,
                                        observed_treatment_sizes,
                                        cluster_latent_var_map,
                                        observed_dewormed_any);

    est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
    sp_est_takeup = inv_logit(hyper_latent_var_map);  
    sp_est_takeup_ate = sp_est_takeup[ate_pairs[, 1]] - sp_est_takeup[ate_pairs[, 2]];
    
    stratum_sp_est_takeup = inv_logit(stratum_latent_var_map); 
    stratum_sp_est_takeup_ate = stratum_sp_est_takeup[, ate_pairs[, 1]] - stratum_sp_est_takeup[, ate_pairs[, 2]];
    
    cluster_sp_est_takeup = inv_logit(cluster_latent_var_map); 
    cluster_sp_est_takeup_ate = cluster_sp_est_takeup[, ate_pairs[, 1]] - cluster_sp_est_takeup[, ate_pairs[, 2]];
  }
}
