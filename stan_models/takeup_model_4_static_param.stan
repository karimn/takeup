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
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs]; // Observations' stratum IDs
  int<lower = 1, upper = num_strata> cluster_stratum_ids[num_clusters]; // Clusters' stratum IDs
  int<lower = 1, upper = num_clusters> strata_cluster_ids[num_clusters]; // Cluster IDs ordered by stratum
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
  
  // real<lower = 0> stratum_student_df;
  real<lower = 0> scale_sigma;
  // real<lower = 0> cluster_student_df;
  real<lower = 0> cluster_scale_sigma;
  
  real<lower = 0> hyper_coef_sigma;
  real<lower = 0> hyper_intercept_sigma;
  
  real<lower = 0> lkj_df;
  
  // Configuration
  
  int<lower = 0, upper = 1> estimate_ate;
}

transformed data {
}

parameters {
  vector[num_all_treatment_coef] hyper_beta;
  
  vector<lower = 0>[num_all_treatment_coef] strata_beta_tau;
  cholesky_factor_corr[num_all_treatment_coef] strata_beta_L_corr_mat;
  matrix<lower = 0>[num_all_treatment_coef, num_strata] cluster_beta_tau;
  cholesky_factor_corr[num_all_treatment_coef] cluster_beta_L_corr_mat[num_strata];
  
  real<lower = 2, upper = 7> stratum_student_df;
  real<lower = 2, upper = 7> cluster_student_df;

  matrix[num_all_treatment_coef, num_strata] strata_beta_raw;
  matrix[num_all_treatment_coef, num_clusters] cluster_beta_raw;
  // // matrix[num_all_treatment_coef, num_strata] QR_strata_beta;
  // 
  // real<lower = 0> cluster_effect_tau;
  // vector[num_clusters] cluster_effect;
}

transformed parameters {
  // matrix[num_all_treatment_coef, num_strata] strata_beta_raw;
  // matrix[num_all_treatment_coef, num_strata] strata_beta; // = R_inv_treatment_design_matrix * QR_strata_beta;
  matrix[num_all_treatment_coef, num_strata] strata_beta; 
  matrix[num_all_treatment_coef, num_clusters] cluster_beta;
  
  {
    int cluster_pos = 1;
    matrix[num_all_treatment_coef, num_all_treatment_coef] strata_beta_L_vcov = diag_pre_multiply(strata_beta_tau, strata_beta_L_corr_mat);
    strata_beta =
      rep_matrix(hyper_beta, num_strata) + strata_beta_L_vcov * strata_beta_raw;

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      matrix[num_all_treatment_coef, num_all_treatment_coef] cluster_beta_L_vcov = diag_pre_multiply(cluster_beta_tau[, stratum_index], 
                                                                                                     cluster_beta_L_corr_mat[stratum_index]);

      cluster_beta[, cluster_pos:cluster_end] = rep_matrix(strata_beta[, stratum_index], curr_num_clusters)
         + cluster_beta_L_vcov * cluster_beta_raw[, cluster_pos:cluster_end];
         // + diag_matrix(cluster_beta_tau[, stratum_index]) * cluster_beta_raw[, cluster_pos:cluster_end];

      cluster_pos = cluster_end + 1;
    }

    // strata_beta_raw = strata_beta_L_vcov \ (strata_beta - rep_matrix(hyper_beta, num_strata));
  }
}

model {
  hyper_beta[1] ~ normal(0, hyper_intercept_sigma);
  hyper_beta[2:num_all_treatment_coef] ~ normal(0, hyper_coef_sigma);
  
  to_vector(strata_beta_raw) ~ student_t(stratum_student_df, 0, 1);
  // to_vector(strata_beta_raw) ~ normal(0, 1);
  to_vector(cluster_beta_raw) ~ student_t(cluster_student_df, 0, 1);
  // to_vector(cluster_beta_raw) ~ normal(0, 1);

  strata_beta_tau ~ normal(0, scale_sigma);
  to_vector(cluster_beta_tau) ~ normal(0, cluster_scale_sigma);

  strata_beta_L_corr_mat ~ lkj_corr_cholesky(lkj_df);

  // cluster_effect_tau ~ normal(0, cluster_intercept_scale_sigma);
  // cluster_effect ~ normal(0, cluster_effect_tau);
  
  {
    vector[num_obs] latent_var = rep_vector(0, num_obs);
   
    int obs_pos = 1; 
    int cluster_pos = 1;

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      cluster_beta_L_corr_mat[stratum_index] ~ lkj_corr_cholesky(lkj_df);

      for (cluster_pos_index in cluster_pos:cluster_end) {
        int curr_cluster_id = strata_cluster_ids[cluster_pos_index];
        int curr_cluster_num_obs = cluster_sizes[curr_cluster_id];
        int obs_end = obs_pos + curr_cluster_num_obs - 1;
          
        latent_var[obs_pos:obs_end] =
          // treatment_design_matrix[obs_pos:obs_end] * hyper_beta;
          // treatment_design_matrix[obs_pos:obs_end] * strata_beta[, stratum_index];
          treatment_design_matrix[obs_pos:obs_end] * cluster_beta[, curr_cluster_id];
          // treatment_design_matrix[obs_pos:obs_end] * strata_beta[, stratum_index] + rep_vector(cluster_effect[curr_cluster_id], curr_cluster_num_obs);
          // Q_treatment_design_matrix[obs_pos:obs_end] * QR_strata_beta[, stratum_index] + rep_vector(cluster_effect[curr_cluster_id], curr_cluster_num_obs);
                           
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
    cluster_latent_var_map = (treatment_map_design_matrix[ate_treatments] * cluster_beta)';
      // (treatment_map_design_matrix[ate_treatments] * rep_matrix(hyper_beta, num_clusters))';
      // (treatment_map_design_matrix[ate_treatments] * cluster_beta)';
      // (treatment_map_design_matrix[ate_treatments] * strata_beta[, cluster_stratum_ids])';
      // + rep_matrix(cluster_effect, num_ate_treatments);

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
