functions {
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
  int<lower = 0> num_obs;
  int<lower = 1> num_all_treatments; // # of treatment cells
  int<lower = 1> num_all_treatment_coef; // # of linear model coefficients minus intercept
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  // Below references to "maps" refer to a finite set of treatment cells that has a one-to-many relationship to actual observations
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix; // Design matrix generated from treatment map
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix;
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
  
  int num_deworming_days;
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
  int<lower = 1, upper = num_deworming_days + 1> dewormed_day_any[num_obs];
 
  int num_dewormed; 
  int strata_dewormed_sizes[num_strata];
  int dewormed_ids[num_dewormed];
  int<lower = 1> stratum_dewormed_index[num_dewormed];
  
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
  matrix[num_relevant_obs_days, num_census_covar_coef] census_covar_dm_long;
 
  // Counterfactuals and ATE
  
  int<lower = 0, upper = num_all_treatments * num_all_treatments> num_ate_treatments;
  int<lower = 1, upper = num_all_treatments> ate_treatments[num_ate_treatments, 2]; 
  
  
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
  
  // Constants for hyperpriors 
  
  real<lower = 0> scale_sigma;
  real<lower = 0> cluster_scale_sigma;
  
  real<lower = 0> hyper_coef_sigma;
  real<lower = 0> hyper_intercept_sigma;
  
  real<lower = 0> lkj_df;
  
  // Configuration
  
  int<lower = 0, upper = 1> estimate_ate;
  int<lower = 0, upper = 1> use_logit; 
  
  // This allows us to specify how many levels the model has:
  // 1: Only the toplevel parameters
  // 2: Strata level (counties)
  // 3: Cluster level (villages)
  int<lower = 1, upper = 3> model_levels;
}

transformed data {
  // Long (dynamic) versions of strata/cluster sizes 
  
  int<lower = 0> relevant_daily_strata_sizes[num_strata];
  int<lower = 0> relevant_daily_cluster_sizes[num_clusters]; // By cluster ID
  
  int<lower = 0> strata_num_all_treatment_coef = model_levels > 1 ? num_all_treatment_coef : 0;
  int<lower = 0> cluster_num_all_treatment_coef = model_levels > 2 ? num_all_treatment_coef : 0;
  int<lower = 0> strata_num_deworming_days = model_levels > 1 ? num_deworming_days : 0;
  int<lower = 0> cluster_num_deworming_days = model_levels > 2 ? num_deworming_days : 0;
  int<lower = 0> strata_num_param_dyn_coef = model_levels > 1 ? num_param_dyn_coef : 0;
  int<lower = 0> cluster_num_param_dyn_coef = model_levels > 1 ? num_param_dyn_coef : 0;
  
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
  
  {
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
  vector[num_all_treatment_coef] hyper_beta_day1;
  vector[num_deworming_days - 1] hyper_baseline_dyn_effect;
  vector[num_param_dyn_coef] hyper_treat_beta_dyn_effect;
  
  real<lower = 2, upper = 7> stratum_day1_student_df;
  real<lower = 2, upper = 7> stratum_baseline_dyn_student_df;
  real<lower = 2, upper = 7> stratum_dyn_student_df;
  
  real<lower = 2, upper = 7> cluster_day1_student_df;
  
  vector<lower = 0>[strata_num_all_treatment_coef] strata_beta_day1_tau;
  cholesky_factor_corr[strata_num_all_treatment_coef] strata_beta_day1_L_corr_mat;
  
  vector<lower = 0>[max(strata_num_deworming_days - 1, 0)] strata_baseline_dyn_effect_tau;
  cholesky_factor_corr[max(strata_num_deworming_days - 1, 0)] strata_baseline_dyn_effect_L_corr_mat;
  
  vector<lower = 0>[strata_num_param_dyn_coef] strata_beta_dyn_effect_tau;
  
  matrix[strata_num_all_treatment_coef, num_strata] strata_beta_day1_raw;
  matrix[max(strata_num_deworming_days - 1, 0), num_strata] strata_baseline_dyn_effect_raw;
  matrix[strata_num_param_dyn_coef, num_strata] QR_strata_beta_dyn_effect;
  
  matrix[cluster_num_all_treatment_coef, num_clusters] cluster_beta_day1;
  matrix<lower = 0>[num_within_cluster_rows, num_strata] cluster_beta_day1_tau;
  cholesky_factor_corr[unique_within_cluster_treatment_sizes[1]] control_cluster_day1_L_corr_mat[num_strata, num_control_cluster_treatments];
  cholesky_factor_corr[unique_within_cluster_treatment_sizes[2]] treated_cluster_day1_L_corr_mat[num_strata, num_treated_cluster_treatments];
}

transformed parameters {
  vector[num_deworming_days] hyper_full_baseline_dyn_effect;  
  matrix[num_deworming_days, num_strata] strata_full_baseline_dyn_effect;  
  matrix[num_deworming_days, num_clusters] cluster_full_baseline_dyn_effect;  
  
  matrix[num_all_treatment_coef, num_clusters] effective_cluster_beta_day1;
  
  matrix[num_param_dyn_coef, num_strata] strata_beta_dyn_effect; 
  matrix[strata_num_param_dyn_coef, num_strata] strata_beta_dyn_effect_raw;  
  
  matrix[num_all_treatment_coef, num_strata] strata_beta_day1;
  
  matrix[num_within_cluster_rows, num_clusters] within_cluster_beta_day1_raw; 

  {
    int cluster_pos = 1;
    
    matrix[num_all_treatment_coef, num_all_treatment_coef] strata_beta_day1_L_vcov;
    matrix[num_deworming_days - 1, num_deworming_days - 1] strata_baseline_dyn_effect_L_vcov;  
      
    hyper_full_baseline_dyn_effect = append_row(0, hyper_baseline_dyn_effect);
      
    if (model_levels > 1) {
      strata_beta_dyn_effect = R_inv_param_dyn_treatment_design_matrix_long * QR_strata_beta_dyn_effect;
      strata_beta_dyn_effect_raw = diag_matrix(strata_beta_dyn_effect_tau) \ (strata_beta_dyn_effect - rep_matrix(hyper_treat_beta_dyn_effect, num_strata));
      
      strata_beta_day1_L_vcov = diag_pre_multiply(strata_beta_day1_tau, strata_beta_day1_L_corr_mat);
      strata_baseline_dyn_effect_L_vcov = diag_pre_multiply(strata_baseline_dyn_effect_tau, strata_baseline_dyn_effect_L_corr_mat);
      
      strata_beta_day1 = rep_matrix(hyper_beta_day1, num_strata) + (strata_beta_day1_L_vcov * strata_beta_day1_raw);
      
      strata_full_baseline_dyn_effect = rep_matrix(hyper_full_baseline_dyn_effect, num_strata)  
        + append_row(rep_row_vector(0, num_strata), strata_baseline_dyn_effect_L_vcov * strata_baseline_dyn_effect_raw);
        
        if (model_levels > 2) {
          effective_cluster_beta_day1 = cluster_beta_day1;
        }
    } else {
      strata_full_baseline_dyn_effect = rep_matrix(hyper_full_baseline_dyn_effect, num_strata);  
      strata_beta_day1 = rep_matrix(hyper_beta_day1, num_strata);
      strata_beta_dyn_effect = rep_matrix(hyper_treat_beta_dyn_effect, num_strata);
    }
                 
    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      if (model_levels > 2) {
        matrix[num_within_cluster_rows, curr_num_clusters] within_cluster_beta_day1 = within_cluster_treatment_map * cluster_beta_day1[, cluster_pos:cluster_end];
        
        matrix[num_within_cluster_rows, num_within_cluster_rows] within_cluster_beta_day1_L_corr_mat = rep_matrix(0, num_within_cluster_rows, num_within_cluster_rows);
        matrix[num_within_cluster_rows, num_within_cluster_rows] within_cluster_beta_day1_L_vcov_mat;  
        vector[num_within_cluster_rows] stratum_cluster_level_mu = within_cluster_treatment_map * strata_beta_day1[, stratum_index];
        
        int control_pos = 1;
        int treated_pos = 1;
        int corr_pos = 1;
        
        for (cluster_level_treat_index in 1:num_cluster_level_treatments) {
          int corr_end; 
          
          if (within_cluster_treatment_sizes[cluster_level_treat_index] == unique_within_cluster_treatment_sizes[1]) {
            corr_end = corr_pos + unique_within_cluster_treatment_sizes[1] - 1;
            within_cluster_beta_day1_L_corr_mat[corr_pos:corr_end, corr_pos:corr_end] = control_cluster_day1_L_corr_mat[stratum_index, control_pos];
            control_pos += 1; 
          } else {
            corr_end = corr_pos + unique_within_cluster_treatment_sizes[2] - 1;
            within_cluster_beta_day1_L_corr_mat[corr_pos:corr_end, corr_pos:corr_end] = treated_cluster_day1_L_corr_mat[stratum_index, treated_pos];
            treated_pos += 1; 
          }
          
          corr_pos = corr_end + 1;
        }
        
        within_cluster_beta_day1_L_vcov_mat = diag_pre_multiply(cluster_beta_day1_tau[, stratum_index], within_cluster_beta_day1_L_corr_mat);
        
        within_cluster_beta_day1_raw[, cluster_pos:cluster_end] = 
          within_cluster_beta_day1_L_vcov_mat \ (within_cluster_beta_day1 - rep_matrix(stratum_cluster_level_mu, curr_num_clusters));
      } else {
        effective_cluster_beta_day1[, cluster_pos:cluster_end] = rep_matrix(strata_beta_day1[, stratum_index], curr_num_clusters);
      } 
      
      cluster_full_baseline_dyn_effect[, cluster_pos:cluster_end] = rep_matrix(strata_full_baseline_dyn_effect[, stratum_index], curr_num_clusters);  
      cluster_pos = cluster_end + 1;
    }
  }
}

model {
  hyper_beta_day1[1] ~ normal(0, hyper_intercept_sigma);
  
  hyper_beta_day1[2:num_all_treatment_coef] ~ normal(0, hyper_coef_sigma);
  
  hyper_baseline_dyn_effect ~ normal(0, hyper_coef_sigma);
  
  to_vector(hyper_treat_beta_dyn_effect) ~ normal(0, hyper_coef_sigma);
 
  if (model_levels > 1) { 
    to_vector(strata_beta_day1_raw) ~ student_t(stratum_day1_student_df, 0, 1);
    to_vector(strata_baseline_dyn_effect_raw) ~ student_t(stratum_baseline_dyn_student_df, 0, 1);
    to_vector(strata_beta_dyn_effect_raw) ~ student_t(stratum_dyn_student_df, 0, 1);
  
    strata_beta_day1_tau ~ normal(0, scale_sigma);
    strata_baseline_dyn_effect_tau ~ normal(0, scale_sigma);
    strata_beta_dyn_effect_tau ~ normal(0, scale_sigma);
  
    strata_beta_day1_L_corr_mat ~ lkj_corr_cholesky(lkj_df);
    strata_baseline_dyn_effect_L_corr_mat ~ lkj_corr_cholesky(lkj_df);
  }
  
  if (model_levels > 2) {
    to_vector(within_cluster_beta_day1_raw) ~ student_t(cluster_day1_student_df, 0, 1);
    to_vector(cluster_beta_day1_tau) ~ normal(0, cluster_scale_sigma);
  }
  
  {
    vector[num_relevant_obs_days] latent_var = rep_vector(0, num_relevant_obs_days);
    
    int cluster_pos = 1;
    int relevant_daily_cluster_pos = 1;

    for (stratum_index in 1:num_strata) {
      int curr_num_clusters = strata_num_clusters[stratum_index];
      int cluster_end = cluster_pos + curr_num_clusters - 1;
      
      if (model_levels > 2) {
        int control_pos = 1;
        int treated_pos = 1;
        
        for (cluster_level_treat_index in 1:num_cluster_level_treatments) {
          if (within_cluster_treatment_sizes[cluster_level_treat_index] == unique_within_cluster_treatment_sizes[1]) {
            control_cluster_day1_L_corr_mat[stratum_index, control_pos] ~ lkj_corr_cholesky(lkj_df);
            control_pos += 1; 
          } else {
            treated_cluster_day1_L_corr_mat[stratum_index, treated_pos] ~ lkj_corr_cholesky(lkj_df);
            treated_pos += 1; 
          }
        }
      }

      for (cluster_pos_index in cluster_pos:cluster_end) {
        int curr_cluster_num_obs = cluster_sizes[cluster_pos_index];
        
        int curr_relevant_daily_cluster_size = relevant_daily_cluster_sizes[cluster_pos_index];
        int relevant_daily_cluster_end = relevant_daily_cluster_pos + curr_relevant_daily_cluster_size - 1;
        
        if (model_levels > 1) {
          latent_var[relevant_daily_cluster_pos:relevant_daily_cluster_end] =
            treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * effective_cluster_beta_day1[, cluster_pos_index] 
            + strata_full_baseline_dyn_effect[dewormed_day_long[relevant_daily_cluster_pos:relevant_daily_cluster_end], stratum_index]
            + Q_param_dyn_treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * QR_strata_beta_dyn_effect[, stratum_index];
        } else {
          latent_var[relevant_daily_cluster_pos:relevant_daily_cluster_end] =
            treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * effective_cluster_beta_day1[, cluster_pos_index] 
            + strata_full_baseline_dyn_effect[dewormed_day_long[relevant_daily_cluster_pos:relevant_daily_cluster_end], stratum_index]
            + param_dyn_treatment_design_matrix_long[relevant_daily_cluster_pos:relevant_daily_cluster_end] * hyper_treat_beta_dyn_effect;
        }
                           
        relevant_daily_cluster_pos = relevant_daily_cluster_end + 1;
      }

      cluster_pos = cluster_end + 1;
    }
    
    if (use_logit) { 
      relevant_dewormed_any_daily ~ bernoulli_logit(latent_var);
    } else {
      relevant_dewormed_any_daily ~ bernoulli(inv_cloglog(latent_var));
    }
  }
}

generated quantities { 
  matrix<lower = 0>[num_ate_treatments, num_deworming_days + 1] est_deworming_days = rep_matrix(0, num_ate_treatments, num_deworming_days + 1);
  vector<lower = 0, upper = 1>[num_ate_treatments] est_takeup = rep_vector(0, num_ate_treatments);
  vector<lower = -1, upper = 1>[num_ate_pairs] est_takeup_ate = rep_vector(0, num_ate_pairs);
  
  
  if (estimate_ate) {
    int stratum_pos = 1;
    
    // vector<lower = 0, upper = 1>[num_deworming_days] sp_est_daily_takeup[num_ate_treatments];
    // matrix<lower = 0, upper = 1>[num_strata, num_deworming_days] stratum_sp_est_daily_takeup[num_ate_treatments];
    // matrix<lower = 0, upper = 1>[num_clusters, num_deworming_days] cluster_sp_est_daily_takeup[num_ate_treatments];
    
    // vector[num_ate_treatments] hyper_latent_var_map_day1;
    // matrix[num_strata, num_ate_treatments] stratum_latent_var_map_day1;
    matrix[num_clusters, num_ate_treatments] cluster_latent_var_map_day1;
    
    // vector[num_deworming_days] hyper_latent_var_map[num_ate_treatments];
    // matrix[num_strata, num_deworming_days] stratum_latent_var_map[num_ate_treatments];
    matrix[num_clusters, num_deworming_days] cluster_latent_var_map[num_ate_treatments];
    
    // hyper_latent_var_map_day1 = treatment_map_design_matrix[ate_treatments[, 1]] * hyper_beta_day1;
    // stratum_latent_var_map_day1 = (treatment_map_design_matrix[ate_treatments[, 1]] * strata_beta_day1)';
    cluster_latent_var_map_day1 = (treatment_map_design_matrix[ate_treatments[, 1]] * effective_cluster_beta_day1)';

    for (ate_treat_index in 1:num_ate_treatments) {
      // hyper_latent_var_map[ate_treat_index] = rep_vector(hyper_latent_var_map_day1[ate_treat_index], num_deworming_days) 
      //   + hyper_full_baseline_dyn_effect
      //   + param_dyn_treatment_map[ate_treatments[ate_treat_index, 2]] * hyper_treat_beta_dyn_effect; 
      //   
      // stratum_latent_var_map[ate_treat_index] = 
      //   rep_matrix(stratum_latent_var_map_day1[, ate_treat_index], num_deworming_days) 
      //   + strata_full_baseline_dyn_effect' 
      //   + (param_dyn_treatment_map[ate_treatments[ate_treat_index, 2]] * strata_beta_dyn_effect)'; 
      
      cluster_latent_var_map[ate_treat_index] = 
        rep_matrix(cluster_latent_var_map_day1[, ate_treat_index], num_deworming_days) 
        + strata_full_baseline_dyn_effect[, cluster_stratum_ids]' 
        + (param_dyn_treatment_map[ate_treatments[ate_treat_index, 2]] * strata_beta_dyn_effect[, cluster_stratum_ids])'; 
        
      // if (use_logit) {
      //   sp_est_daily_takeup[ate_treat_index] = inv_logit(hyper_latent_var_map[ate_treat_index]); 
      //   stratum_sp_est_daily_takeup[ate_treat_index] = inv_logit(stratum_latent_var_map[ate_treat_index]); 
      //   cluster_sp_est_daily_takeup[ate_treat_index] = inv_logit(cluster_latent_var_map[ate_treat_index]); 
      // } else {
      //   for (day_index in 1:num_deworming_days) {
      //     sp_est_daily_takeup[ate_treat_index, day_index] = 1 - gumbel_cdf(- hyper_latent_var_map[ate_treat_index, day_index], 0, 1); 
      //     
      //     for (stratum_index in 1:num_strata) {
      //       stratum_sp_est_daily_takeup[ate_treat_index, stratum_index, day_index] = 
      //         1 - gumbel_cdf(- stratum_latent_var_map[ate_treat_index, stratum_index, day_index], 0, 1); 
      //     }
      //     
      //     for (cluster_index in 1:num_clusters) {
      //       cluster_sp_est_daily_takeup[ate_treat_index, cluster_index, day_index] = 
      //         1 - gumbel_cdf(- cluster_latent_var_map[ate_treat_index, cluster_index, day_index], 0, 1); 
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
                                            cluster_latent_var_map,
                                            observed_dewormed_day,
                                            use_logit);

    est_takeup = 1 - est_deworming_days[, 13];
    est_takeup_ate = est_takeup[ate_pairs[, 1]] - est_takeup[ate_pairs[, 2]];
  }
}
