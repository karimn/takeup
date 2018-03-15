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
  
  matrix treatment_deworming_day_rng(matrix missing_latent_var, int use_logit) {
    int curr_treatment_size = rows(missing_latent_var); 
    int num_deworming_days = cols(missing_latent_var);
    
    matrix[curr_treatment_size, num_deworming_days + 1] deworming_days_mask = rep_matrix(0, curr_treatment_size, num_deworming_days + 1); 
    
    for (obs_ids_index in 1:curr_treatment_size) {
      vector[num_deworming_days + 1] deworming_day_prob = rep_vector(1, num_deworming_days + 1);
      vector[num_deworming_days + 1] cumul_prod_alphas = rep_vector(1, num_deworming_days + 1);
    
      real prev_not_deworm_prob = 1; 
      
      for (day_index in 1:num_deworming_days) {
        cumul_prod_alphas[day_index] = prev_not_deworm_prob;
        
        if (use_logit) {
          deworming_day_prob[day_index] = inv_logit(missing_latent_var[obs_ids_index, day_index]);
        } else {
          deworming_day_prob[day_index] = 1 - gumbel_cdf(- missing_latent_var[obs_ids_index, day_index], 0, 1); 
        }
        
        prev_not_deworm_prob = prev_not_deworm_prob * (1 - deworming_day_prob[day_index]);
      }
      
      cumul_prod_alphas[num_deworming_days + 1] = prev_not_deworm_prob;
      
      deworming_day_prob = deworming_day_prob .* cumul_prod_alphas; 
      
      deworming_days_mask[obs_ids_index, categorical_rng(deworming_day_prob)] = 1;
    }
  
    return deworming_days_mask; 
  }
  
  matrix treatment_cell_deworming_day_prop_rng(int[,] treatment_ids, 
                                               int[] missing_obs_ids, 
                                               // int[] missing_cluster_id,
                                               int[] missing_level_id, // Either stratum or cluster
                                               int[] missing_treatment_sizes,
                                               int[] observed_treatment_sizes,
                                               // matrix census_covar_latent_var,
                                               matrix[] latent_var_map, // matrix[num_clusters, num_deworming_days][num_ate_treatments]
                                               int[] observed_dewormed_any,
                                               int use_logit) {
    int num_treatment_ids = size(treatment_ids);
    int num_deworming_days = cols(latent_var_map[1]);
    int missing_treatment_pos = 1;
    int observed_treatment_pos = 1;
    matrix[num_treatment_ids, num_deworming_days + 1] cell_deworming_day = rep_matrix(0, num_treatment_ids, num_deworming_days + 1);
    
    matrix[num_deworming_days + 1, num_deworming_days + 1] days_diag = diag_matrix(rep_vector(1, num_deworming_days + 1));
    
    for (treatment_ids_index in 1:num_treatment_ids) {
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_ids_index];
      int missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      int curr_observed_treatment_size = observed_treatment_sizes[treatment_ids_index];
      
      matrix[curr_missing_treatment_size, num_deworming_days + 1] missing_deworming_days_mask = 
        treatment_deworming_day_rng(latent_var_map[treatment_ids_index, missing_level_id[missing_treatment_pos:missing_treatment_end]], use_logit);
                                    
      if (curr_observed_treatment_size > 0) {
        int observed_treatment_end = observed_treatment_pos + curr_observed_treatment_size - 1;
        
        matrix[curr_observed_treatment_size, num_deworming_days + 1] observed_deworming_days_mask = 
          days_diag[observed_dewormed_any[observed_treatment_pos:observed_treatment_end]];
      
        cell_deworming_day[treatment_ids_index, 1:(num_deworming_days + 1)] = (diagonal(crossprod(missing_deworming_days_mask)) 
          + diagonal(crossprod(observed_deworming_days_mask)))' / (curr_missing_treatment_size + curr_observed_treatment_size);
          
        observed_treatment_pos = observed_treatment_end + 1;
      } else {
        cell_deworming_day[treatment_ids_index, 1:(num_deworming_days + 1)] = 
          diagonal(crossprod(missing_deworming_days_mask))' / curr_missing_treatment_size;
      }
        
      missing_treatment_pos = missing_treatment_end + 1;
    }
    
    return(cell_deworming_day);
  }
}

data {
  // Configuration
  
  int<lower = 0, upper = 1> dynamic_model;
 
  // 0: complementary log-log (for the dynamic model)
  // 1: logistic
  int<lower = 0, upper = 1> model_link_type; 
  
  int<lower = 0, upper = 1> estimate_ate;
  
  // This allows us to specify how many levels the model has:
  // 1: Only the toplevel parameters
  // 2: Strata level (counties)
  // 3: Cluster level (villages)
  int<lower = 1, upper = 3> model_levels;
  
  int<lower = 0, upper = 1> use_cluster_re; // If model_levels < 2 this flag determines if a cluster RE is used.
  int<lower = 0, upper = 1> use_cluster_identity_corr; // Uncorrelated parameters at the cluster level (must have model_levels > 2)
  
  // Sizes
  
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
  
  int<lower = 0, upper = num_all_treatment_coef> num_private_value_dist_col;
  int<lower = 1, upper = num_all_treatment_coef> private_value_dist_col[num_private_value_dist_col];
  int<lower = 1, upper = num_all_treatment_coef> not_private_value_dist_col[num_all_treatment_coef - num_private_value_dist_col - 1];
  
  int<lower = 1, upper = num_all_treatments> num_cluster_level_treatments;
  int<lower = 1> within_cluster_treatment_sizes[num_cluster_level_treatments];
  int<lower = 0, upper = num_cluster_level_treatments> unique_within_cluster_treatment_sizes[2]; // Control and treated clusters have different sizes
  matrix<lower = 0, upper = 1>[sum(within_cluster_treatment_sizes), num_all_treatment_coef] within_cluster_treatment_map;
  // Left inverse of within_cluster_treatment_map
  matrix[num_all_treatment_coef, sum(within_cluster_treatment_sizes)] within_cluster_treatment_map_ginv; 
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs]; // Observations' stratum IDs
  int<lower = 1, upper = num_strata> cluster_stratum_ids[num_clusters]; // Clusters' stratum IDs
  int<lower = 1, upper = num_obs> cluster_obs_ids[num_obs]; // Observation IDs ordered by cluster
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; // Observations in strata 
  int<lower = 1, upper = num_clusters> strata_num_clusters[num_strata]; // Clusters in strata
  int<lower = 1, upper = num_obs> cluster_sizes[num_clusters]; // Observations in clusters
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
  
  // Deworming outcomes
  
  int num_deworming_days;
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
  int<lower = 1, upper = num_deworming_days + 1> dewormed_day_any[num_obs];
  
  // Dynamics
  
  int<lower = 1> param_poly_order;
  
  int<lower = 0, upper = 1> relevant_latent_var_map[num_deworming_days + 1, num_deworming_days];
  
  int<lower = 1, upper = num_obs * num_deworming_days> num_relevant_obs_days;
  int<lower = 1> num_param_dyn_coef;
  
  int<lower = 1, upper = num_clusters> cluster_id_long[num_relevant_obs_days];
  
  int<lower = 0, upper = 1> relevant_dewormed_any_daily[num_relevant_obs_days];
  int<lower = 1, upper = num_deworming_days> dewormed_day_long[num_relevant_obs_days];
  
  matrix[num_relevant_obs_days, num_param_dyn_coef] param_dyn_treatment_design_matrix_long;
  matrix[num_relevant_obs_days, num_param_dyn_coef] Q_param_dyn_treatment_design_matrix_long;
  matrix[num_param_dyn_coef, num_param_dyn_coef] R_inv_param_dyn_treatment_design_matrix_long;
  
  matrix[num_relevant_obs_days, num_all_treatment_coef] treatment_design_matrix_long;
  matrix[num_relevant_obs_days, num_all_treatment_coef] Q_treatment_design_matrix_long;
  matrix[num_all_treatment_coef, num_all_treatment_coef] R_inv_treatment_design_matrix_long;
  
  matrix[num_deworming_days, num_param_dyn_coef] param_dyn_treatment_map[num_all_treatments];
  
  // Covariates
  
  int<lower = 1> num_census_covar_coef;
  matrix[num_obs, num_census_covar_coef] census_covar_dm;
 
  // Counterfactuals and ATE
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_treatments;
  int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments, dynamic_model + 1]; 
  
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
  
  int observed_dewormed_day[num_observed_obs_ids];
  int observed_dewormed_any[num_observed_obs_ids];
  
  // Constants for hyperpriors 
  
  real<lower = 0> scale_sigma;
  real<lower = 0> cluster_scale_sigma;
  
  real<lower = 0> hyper_coef_sigma;
  real<lower = 0> hyper_intercept_sigma;
  
  real<lower = 0> private_value_dist_sigma;
  
  real<lower = 0> lkj_df;
  
}

transformed data {
  // Long (dynamic) versions of strata/cluster sizes 
  int<lower = 0> relevant_daily_strata_sizes[num_strata];
  int<lower = 0> relevant_daily_cluster_sizes[num_clusters]; // By cluster ID
  
  
  int<lower = 0> strata_num_all_treatment_coef = model_levels > 1 ? num_all_treatment_coef : 0;
  int<lower = 0> cluster_num_all_treatment_coef = 0;
  
  int<lower = 0, upper = num_all_treatments> num_within_cluster_rows = 0;
  int<lower = 0, upper = num_all_treatments> num_cluster_tau = 0;
  int<lower = 0, upper = num_all_treatments> num_true_cluster_tau = 0;
  int<lower = 0, upper = num_cluster_level_treatments> num_control_cluster_treatments = 0;
  int<lower = 0, upper = num_cluster_level_treatments> num_treated_cluster_treatments = 0;
  int<lower = 0, upper = num_all_treatments> num_within_cluster_control_treatments = 0;
  int<lower = 0, upper = num_all_treatments> num_within_cluster_treated_treatments = 0;
  
  int<lower = 0> strata_num_deworming_days = dynamic_model && model_levels > 1 ? num_deworming_days : 0;
  int<lower = 0> cluster_num_deworming_days = dynamic_model && model_levels > 2 ? num_deworming_days : 0;
  int<lower = 0> strata_num_param_dyn_coef = dynamic_model && model_levels > 1 ? num_param_dyn_coef : 0;
  int<lower = 0> cluster_num_param_dyn_coef = dynamic_model && model_levels > 1 ? num_param_dyn_coef : 0;
 
  if (model_levels > 2) { 
    num_treated_cluster_treatments = num_cluster_level_treatments;
    num_within_cluster_rows = sum(within_cluster_treatment_sizes);
    num_cluster_tau = num_within_cluster_rows;
    num_true_cluster_tau = num_all_treatment_coef;
    
    cluster_num_all_treatment_coef = num_all_treatment_coef;
    
    if (!use_cluster_identity_corr) {
      num_within_cluster_control_treatments = unique_within_cluster_treatment_sizes[1];
      num_within_cluster_treated_treatments = unique_within_cluster_treatment_sizes[2];
    }
    
    for (cluster_level_treat_index in 1:num_cluster_level_treatments) {
      if (within_cluster_treatment_sizes[cluster_level_treat_index] == unique_within_cluster_treatment_sizes[1]) {
        num_control_cluster_treatments += 1;
        num_treated_cluster_treatments -= 1;
      }
    }
  } else if (use_cluster_re) {
    num_cluster_tau = 1; 
    cluster_num_all_treatment_coef = 1;
  }
  
  { // Dynamic model data
    int stratum_pos = 1;
    int cluster_pos = 1;
    int cluster_obs_pos = 1;
   
    for (stratum_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[stratum_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;

      for (cluster_pos_index in cluster_pos:cluster_end) {
        int curr_cluster_num_obs = cluster_sizes[cluster_pos_index];
        
        int cluster_obs_end = cluster_obs_pos + curr_cluster_num_obs - 1;
        
        relevant_daily_cluster_sizes[cluster_pos_index] = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any[cluster_obs_pos:cluster_obs_end]]));
        
        cluster_obs_pos = cluster_obs_end + 1;
      }

      relevant_daily_strata_sizes[stratum_index] = sum(to_array_1d(relevant_latent_var_map[dewormed_day_any[stratum_pos:stratum_end]]));

      stratum_pos = stratum_end + 1;
      cluster_pos = cluster_end + 1;
    }

    if (sum(relevant_daily_strata_sizes) != num_relevant_obs_days) {
      reject("Unexpected relevant strata sizes ", sum(relevant_daily_strata_sizes), " ", num_relevant_obs_days);
    }
  }
}

parameters {
  vector[num_all_treatment_coef] hyper_beta;
  vector[dynamic_model ? num_deworming_days - 1 : 0] hyper_baseline_dyn_effect; // Dynamic
  vector[dynamic_model ? num_param_dyn_coef : 0] hyper_treat_beta_dyn_effect; // Dynamic
 
  // Student t degrees of freedom parameters with uniform priors
  real<lower = 2, upper = 8> stratum_student_df;
  real<lower = 2, upper = 8> stratum_baseline_dyn_student_df; // Dynamic
  real<lower = 2, upper = 8> stratum_dyn_student_df; // Dynamic
  real<lower = 2, upper = 8> cluster_student_df;
  
  matrix[strata_num_all_treatment_coef, num_strata] strata_beta_raw;
  vector<lower = 0>[strata_num_all_treatment_coef] strata_beta_tau;
  cholesky_factor_corr[strata_num_all_treatment_coef] strata_beta_L_corr_mat;
  
  vector<lower = 0>[max(strata_num_deworming_days - 1, 0)] strata_baseline_dyn_effect_tau; // Dynamic
  cholesky_factor_corr[max(strata_num_deworming_days - 1, 0)] strata_baseline_dyn_effect_L_corr_mat; // Dynamic
  
  vector<lower = 0>[strata_num_param_dyn_coef] strata_beta_dyn_effect_tau; // Dynamic
  matrix[max(strata_num_deworming_days - 1, 0), num_strata] strata_baseline_dyn_effect_raw; // Dynamic
  matrix[strata_num_param_dyn_coef, num_strata] QR_strata_beta_dyn_effect; // Dynamic
  
  matrix[cluster_num_all_treatment_coef, num_clusters] cluster_beta;
  matrix<lower = 0>[num_cluster_tau, num_strata] cluster_beta_tau;
  cholesky_factor_corr[num_within_cluster_control_treatments] control_cluster_L_corr_mat[num_strata, num_control_cluster_treatments];
  cholesky_factor_corr[num_within_cluster_treated_treatments] treated_cluster_L_corr_mat[num_strata, num_treated_cluster_treatments];
}

transformed parameters {
  vector[dynamic_model ? num_deworming_days : 0] hyper_full_baseline_dyn_effect;  
  matrix[dynamic_model ? num_deworming_days : 0, num_strata] strata_full_baseline_dyn_effect = 
    rep_matrix(0, dynamic_model ? num_deworming_days : 0, num_strata);  
  matrix[dynamic_model ? num_deworming_days : 0, num_clusters] cluster_full_baseline_dyn_effect;  
  
  // These are always specified even if the model is not multilevel; they'll just have the same values one level up.
  matrix[num_all_treatment_coef, num_strata] strata_beta = rep_matrix(0, num_all_treatment_coef, num_strata); 
  matrix[num_all_treatment_coef, num_clusters] effective_cluster_beta = rep_matrix(0, num_all_treatment_coef, num_clusters);
  matrix[num_within_cluster_rows, num_clusters] within_cluster_beta_raw;
  matrix<lower = 0>[num_true_cluster_tau, num_strata] true_cluster_beta_tau;
  
  matrix[dynamic_model ? num_param_dyn_coef : 0, num_strata] strata_beta_dyn_effect; 
  matrix[strata_num_param_dyn_coef, num_strata] strata_beta_dyn_effect_raw;  
  
  cholesky_factor_cov[strata_num_all_treatment_coef] strata_beta_L_vcov = diag_pre_multiply(strata_beta_tau, strata_beta_L_corr_mat);
  
  {
    int cluster_pos = 1;
    matrix[dynamic_model ? num_deworming_days - 1 : 0, dynamic_model ? num_deworming_days - 1 : 0] strata_baseline_dyn_effect_L_vcov;  
   
    if (dynamic_model) { 
      hyper_full_baseline_dyn_effect = append_row(0, hyper_baseline_dyn_effect);
    }
    
    if (model_levels > 1) {
      strata_beta = rep_matrix(hyper_beta, num_strata) + strata_beta_L_vcov * strata_beta_raw;
      
      if (dynamic_model) { 
        strata_beta_dyn_effect = R_inv_param_dyn_treatment_design_matrix_long * QR_strata_beta_dyn_effect;
        strata_beta_dyn_effect_raw = diag_matrix(strata_beta_dyn_effect_tau) \ (strata_beta_dyn_effect - rep_matrix(hyper_treat_beta_dyn_effect, num_strata));
        
        strata_baseline_dyn_effect_L_vcov = diag_pre_multiply(strata_baseline_dyn_effect_tau, strata_baseline_dyn_effect_L_corr_mat);
        
        strata_full_baseline_dyn_effect = rep_matrix(hyper_full_baseline_dyn_effect, num_strata)  
          + append_row(rep_row_vector(0, num_strata), strata_baseline_dyn_effect_L_vcov * strata_baseline_dyn_effect_raw);
      }
    } else {
      strata_beta = rep_matrix(hyper_beta, num_strata);
      
      if (dynamic_model) { 
        strata_full_baseline_dyn_effect = rep_matrix(hyper_full_baseline_dyn_effect, num_strata);  
        strata_beta_dyn_effect = rep_matrix(hyper_treat_beta_dyn_effect, num_strata);
      }
    }

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      if (model_levels > 2) {
        matrix[num_within_cluster_rows, num_within_cluster_rows] within_cluster_beta_L_corr_mat = rep_matrix(0, num_within_cluster_rows, num_within_cluster_rows);
        matrix[num_within_cluster_rows, num_within_cluster_rows] within_cluster_beta_L_vcov_mat;  
        matrix[num_all_treatment_coef, num_all_treatment_coef] true_cluster_beta_vcov_mat;  
        vector[num_within_cluster_rows] stratum_cluster_level_mu = within_cluster_treatment_map * strata_beta[, stratum_index];
        
        matrix[num_within_cluster_rows, curr_num_clusters] within_cluster_beta = within_cluster_treatment_map * cluster_beta[, cluster_pos:cluster_end];
        
        if (use_cluster_identity_corr) {
          within_cluster_beta_L_corr_mat = diag_matrix(rep_vector(1, num_within_cluster_rows));
        } else {
          int control_pos = 1;
          int treated_pos = 1;
          int corr_pos = 1;
          
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
        }
      
        within_cluster_beta_L_vcov_mat = diag_pre_multiply(cluster_beta_tau[, stratum_index], within_cluster_beta_L_corr_mat);
        true_cluster_beta_vcov_mat = tcrossprod(within_cluster_treatment_map_ginv * within_cluster_beta_L_vcov_mat);
        true_cluster_beta_tau[, stratum_index] = diagonal(true_cluster_beta_vcov_mat);
        
        within_cluster_beta_raw[, cluster_pos:cluster_end] = 
          within_cluster_beta_L_vcov_mat \ (within_cluster_beta - rep_matrix(stratum_cluster_level_mu, curr_num_clusters));
      } else if (use_cluster_re) {
        effective_cluster_beta[1, cluster_pos:cluster_end] = 
          strata_beta[1, stratum_index] + cluster_beta_tau[1, stratum_index] * cluster_beta[1, cluster_pos:cluster_end];
        
        effective_cluster_beta[2:num_all_treatment_coef, cluster_pos:cluster_end] = 
          rep_matrix(strata_beta[2:num_all_treatment_coef, stratum_index], curr_num_clusters);
      }  else {
        effective_cluster_beta[, cluster_pos:cluster_end] = rep_matrix(strata_beta[, stratum_index], curr_num_clusters);
      }
      
      cluster_full_baseline_dyn_effect[, cluster_pos:cluster_end] = rep_matrix(strata_full_baseline_dyn_effect[, stratum_index], curr_num_clusters);  
      
      cluster_pos = cluster_end + 1;
    }
    
    if (model_levels > 2 && !use_cluster_re) {
      effective_cluster_beta = cluster_beta;
    }  
  }
}

model {
  hyper_beta[1] ~ normal(0, hyper_intercept_sigma);
  hyper_beta[not_private_value_dist_col] ~ normal(0, hyper_coef_sigma);
  hyper_beta[private_value_dist_col] ~ normal(0, private_value_dist_sigma);
  
  if (dynamic_model) {
    hyper_baseline_dyn_effect ~ normal(0, hyper_coef_sigma);
    to_vector(hyper_treat_beta_dyn_effect) ~ normal(0, hyper_coef_sigma);
  }
 
  if (model_levels > 1) { 
    if (dynamic_model) {
      to_vector(strata_baseline_dyn_effect_raw) ~ student_t(stratum_baseline_dyn_student_df, 0, 1);
      to_vector(strata_beta_dyn_effect_raw) ~ student_t(stratum_dyn_student_df, 0, 1);
      strata_baseline_dyn_effect_tau ~ normal(0, scale_sigma);
      strata_beta_dyn_effect_tau ~ normal(0, scale_sigma);
      strata_baseline_dyn_effect_L_corr_mat ~ lkj_corr_cholesky(lkj_df);
    }
    
    to_vector(strata_beta_raw) ~ student_t(stratum_student_df, 0, 1);
    strata_beta_tau ~ normal(0, scale_sigma);
    strata_beta_L_corr_mat ~ lkj_corr_cholesky(lkj_df);
 
    if (model_levels > 2 || use_cluster_re) { 
      to_vector(within_cluster_beta_raw) ~ student_t(cluster_student_df, 0, 1);
      to_vector(true_cluster_beta_tau) ~ normal(0, cluster_scale_sigma);
    } else if (use_cluster_re) {
      to_vector(cluster_beta) ~ student_t(cluster_student_df, 0, 1);
      to_vector(cluster_beta_tau) ~ normal(0, cluster_scale_sigma);
    }
  }
  
  {
    vector[dynamic_model ? num_relevant_obs_days : num_obs] latent_var = rep_vector(0, dynamic_model ? num_relevant_obs_days : num_obs);
   
    int obs_pos = 1; 
    int cluster_pos = 1;
    int relevant_daily_cluster_pos = 1;

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      int control_pos = 1;
      int treated_pos = 1;
      
      if (model_levels > 2 && !use_cluster_identity_corr) {
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
        
        if (dynamic_model) { 
          int curr_relevant_daily_cluster_size = relevant_daily_cluster_sizes[cluster_pos_index];
          int relevant_daily_cluster_end = relevant_daily_cluster_pos + curr_relevant_daily_cluster_size - 1;
          
          if (model_levels > 1) {
            latent_var[relevant_daily_cluster_pos:relevant_daily_cluster_end] =
              treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * effective_cluster_beta[, cluster_pos_index]
              + strata_full_baseline_dyn_effect[dewormed_day_long[relevant_daily_cluster_pos:relevant_daily_cluster_end], stratum_index]
            + Q_param_dyn_treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * QR_strata_beta_dyn_effect[, stratum_index];
          } else {
            latent_var[relevant_daily_cluster_pos:relevant_daily_cluster_end] =
              treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * effective_cluster_beta[, cluster_pos_index] 
              + strata_full_baseline_dyn_effect[dewormed_day_long[relevant_daily_cluster_pos:relevant_daily_cluster_end], stratum_index]
              + param_dyn_treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * hyper_treat_beta_dyn_effect;
          }
                           
          relevant_daily_cluster_pos = relevant_daily_cluster_end + 1;
        } else {
          int obs_end = obs_pos + curr_cluster_num_obs - 1;
       
          latent_var[obs_pos:obs_end] =
            treatment_design_matrix[obs_pos:obs_end] * effective_cluster_beta[, cluster_pos_index];
            
          obs_pos = obs_end + 1;
        }
      }

      cluster_pos = cluster_end + 1;
    }
    
    if (dynamic_model) {
      if (model_link_type == 0) { 
        relevant_dewormed_any_daily ~ bernoulli(inv_cloglog(latent_var));
      } else if (model_link_type == 1) {
        relevant_dewormed_any_daily ~ bernoulli_logit(latent_var);
      }
    } else {
      if (model_link_type == 0) { 
        dewormed_any ~ bernoulli(inv_cloglog(latent_var));
      } else if (model_link_type == 1) {
        dewormed_any ~ bernoulli_logit(latent_var);
      }
    }
  }
}

generated quantities { 
  // Finite sample estimands
  matrix<lower = 0>[num_ate_treatments, dynamic_model ? num_deworming_days + 1 : 0] est_deworming_days;  
  vector<lower = 0, upper = 1>[num_ate_treatments] est_takeup = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] est_takeup_ate = rep_vector(0, num_ate_pairs);
  
  // vector<lower = 0, upper = 1>[num_deworming_days] sp_est_daily_takeup[(estimate_ate && save_sp_estimates) ? num_ate_treatments : 0];
  // matrix<lower = 0, upper = 1>[num_strata, num_deworming_days] stratum_sp_est_daily_takeup[(estimate_ate && save_sp_estimates) ? num_ate_treatments : 0];
  // matrix<lower = 0, upper = 1>[num_clusters, num_deworming_days] cluster_sp_est_daily_takeup[(estimate_ate && save_sp_estimates) ? num_ate_treatments : 0];
 
  // Superpopulation estimands 
  vector<lower = 0, upper = 1>[num_ate_treatments] sp_est_takeup = rep_vector(0, num_ate_treatments);
  matrix<lower = 0, upper = 1>[num_strata, num_ate_treatments] stratum_sp_est_takeup = rep_matrix(0, num_strata, num_ate_treatments);
  matrix<lower = 0, upper = 1>[num_clusters, num_ate_treatments] cluster_sp_est_takeup = rep_matrix(0, num_clusters, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] sp_est_takeup_ate = rep_vector(0, num_ate_pairs);
  matrix<lower = -1, upper = 1>[num_strata, num_ate_pairs] stratum_sp_est_takeup_ate = rep_matrix(0, num_strata, num_ate_pairs);
  matrix<lower = -1, upper = 1>[num_clusters, num_ate_pairs] cluster_sp_est_takeup_ate = rep_matrix(0, num_clusters, num_ate_pairs);

  corr_matrix[num_ate_treatments] strata_treatment_corr = rep_matrix(0, num_ate_treatments, num_ate_treatments);
  
  vector[num_ate_treatments] hyper_latent_var_map; 
  matrix[num_strata, num_ate_treatments] stratum_latent_var_map; 
  
  if (dynamic_model) {
    est_deworming_days = rep_matrix(0, num_ate_treatments, num_deworming_days + 1);
  }

  if (estimate_ate) {
    int stratum_pos = 1;
    matrix[num_clusters, num_ate_treatments] cluster_latent_var_map; 
    matrix[num_all_treatment_coef, num_all_treatment_coef] strata_beta_vcov = tcrossprod(strata_beta_L_vcov); 
    matrix[num_ate_treatments, num_ate_treatments] strata_treatment_sandwich = (diag_matrix(diagonal(strata_beta_vcov)) \ treatment_map_design_matrix[ate_treatments[, 1]])';
    
    strata_treatment_corr = quad_form(strata_beta_vcov, strata_treatment_sandwich);
    
    hyper_latent_var_map = treatment_map_design_matrix[ate_treatments[, 1]] * hyper_beta;
    stratum_latent_var_map = (treatment_map_design_matrix[ate_treatments[, 1]] * strata_beta)';
    cluster_latent_var_map = (treatment_map_design_matrix[ate_treatments[, 1]] * effective_cluster_beta)';
    
    if (dynamic_model) {
      // vector[num_deworming_days] hyper_daily_latent_var_map[num_ate_treatments];
      // matrix[num_strata, num_deworming_days] stratum_daily_latent_var_map[num_ate_treatments];
      matrix[num_clusters, num_deworming_days] cluster_daily_latent_var_map[num_ate_treatments];
      
      for (ate_treat_index in 1:num_ate_treatments) {
        // if (save_sp_estimates) {
        //   hyper_latent_var_map[ate_treat_index] = rep_vector(hyper_latent_var_map_day1[ate_treat_index], num_deworming_days)
        //     + hyper_full_baseline_dyn_effect
        //     + param_dyn_treatment_map[ate_treatments[ate_treat_index, 2]] * hyper_treat_beta_dyn_effect;
        // 
        //   stratum_latent_var_map[ate_treat_index] =
        //     rep_matrix(stratum_latent_var_map_day1[, ate_treat_index], num_deworming_days)
        //     + strata_full_baseline_dyn_effect'
        //     + (param_dyn_treatment_map[ate_treatments[ate_treat_index, 2]] * strata_beta_dyn_effect)';
        // }
        
        cluster_daily_latent_var_map[ate_treat_index] = 
          rep_matrix(cluster_latent_var_map[, ate_treat_index], num_deworming_days) 
          + strata_full_baseline_dyn_effect[, cluster_stratum_ids]' 
          + (param_dyn_treatment_map[ate_treatments[ate_treat_index, 2]] * strata_beta_dyn_effect[, cluster_stratum_ids])'; 
        
        // if (save_sp_estimates) {  
        //   if (use_logit) {
        //     sp_est_daily_takeup[ate_treat_index] = inv_logit(hyper_latent_var_map[ate_treat_index]);
        //     stratum_sp_est_daily_takeup[ate_treat_index] = inv_logit(stratum_latent_var_map[ate_treat_index]);
        //     cluster_sp_est_daily_takeup[ate_treat_index] = inv_logit(cluster_latent_var_map[ate_treat_index]);
        //   } else {
        //     for (day_index in 1:num_deworming_days) {
        //       sp_est_daily_takeup[ate_treat_index, day_index] = 1 - gumbel_cdf(- hyper_latent_var_map[ate_treat_index, day_index], 0, 1);
        // 
        //       for (stratum_index in 1:num_strata) {
        //         stratum_sp_est_daily_takeup[ate_treat_index, stratum_index, day_index] =
        //           1 - gumbel_cdf(- stratum_latent_var_map[ate_treat_index, stratum_index, day_index], 0, 1);
        //       }
        // 
        //       for (cluster_index in 1:num_clusters) {
        //         cluster_sp_est_daily_takeup[ate_treat_index, cluster_index, day_index] =
        //           1 - gumbel_cdf(- cluster_latent_var_map[ate_treat_index, cluster_index, day_index], 0, 1);
        //       }
        //     }
        //   }
        // }
      }
      
      est_deworming_days =
        treatment_cell_deworming_day_prop_rng(ate_treatments,
                                              missing_obs_ids,
                                              missing_treatment_cluster_id,
                                              missing_treatment_sizes,
                                              observed_treatment_sizes,
                                              cluster_daily_latent_var_map,
                                              observed_dewormed_day,
                                              model_link_type == 1);
  
      est_takeup = 1 - est_deworming_days[, num_deworming_days + 1];
    } else {
      est_takeup =
        treatment_cell_deworming_prop_rng(num_ate_treatments,
                                          missing_obs_ids,
                                          missing_treatment_cluster_id,
                                          missing_treatment_sizes,
                                          observed_treatment_sizes,
                                          cluster_latent_var_map,
                                          observed_dewormed_any);
                                          
      sp_est_takeup = inv_logit(hyper_latent_var_map);  
      sp_est_takeup_ate = sp_est_takeup[ate_pairs[, 1]] - sp_est_takeup[ate_pairs[, 2]];
      
      stratum_sp_est_takeup = inv_logit(stratum_latent_var_map); 
      stratum_sp_est_takeup_ate = stratum_sp_est_takeup[, ate_pairs[, 1]] - stratum_sp_est_takeup[, ate_pairs[, 2]];
      
      cluster_sp_est_takeup = inv_logit(cluster_latent_var_map); 
      cluster_sp_est_takeup_ate = cluster_sp_est_takeup[, ate_pairs[, 1]] - cluster_sp_est_takeup[, ate_pairs[, 2]];
    }

    est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }
}
