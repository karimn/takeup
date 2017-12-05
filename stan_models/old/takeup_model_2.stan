functions {
  matrix treatment_cell_takeup_rng(int[] treatment_ids, 
                                   int[] missing_obs_ids, 
                                   int[] missing_stratum_id,
                                   int[] missing_cluster_id,
                                   vector missing_bracelet_util_diff,
                                   int[] private_value_calendar_coef,
                                   int[] private_value_bracelet_coef,
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
        
        (stratum_treatment_coef[missing_stratum_id[missing_treatment_pos:missing_treatment_end], private_value_calendar_coef] *
           treatment_map_dm[treatment_ids[treatment_ids_index], private_value_bracelet_coef]') .* 
          (1 + missing_bracelet_util_diff[missing_treatment_pos:missing_treatment_end]) +
          
        missing_name_matched[missing_treatment_pos:missing_treatment_end] .* stratum_name_matching_effect[missing_stratum_id[missing_treatment_pos:missing_treatment_end]] +
        (diag_treatment_map_dm[treatment_ids[treatment_ids_index]] * stratum_treatment_name_matching_interact[, missing_stratum_id[missing_treatment_pos:missing_treatment_end]])' .*
          missing_name_matched[missing_treatment_pos:missing_treatment_end];
      
      vector[curr_missing_treatment_size] all_takeup;
      real all_observed_takeup_total = sum(observed_dewormed_any[observed_treatment_pos:observed_treatment_end]); 
      real monitored_observed_takeup_total = 
          sum(observed_dewormed_any[observed_treatment_pos:observed_treatment_end] .* (1 - observed_name_matched[observed_treatment_pos:observed_treatment_end]));
          
      real all_cell_num_obs = curr_missing_treatment_size + curr_observed_treatment_size;
      real monitored_cell_num_obs = all_cell_num_obs - 
        sum(missing_name_matched[missing_treatment_pos:missing_treatment_end]) - sum(observed_name_matched[observed_treatment_pos:observed_treatment_end]);
      
      for (missing_obs_ids_index in 1:curr_missing_treatment_size) {
        all_takeup[missing_obs_ids_index] = bernoulli_logit_rng(missing_latent_utility[missing_obs_ids_index]);  
      }
     
      cell_takeup_prop[treatment_ids_index, 1] = (sum(all_takeup) + all_observed_takeup_total) / all_cell_num_obs;
      cell_takeup_prop[treatment_ids_index, 2] = 
        (sum(all_takeup .* (1 - missing_name_matched[missing_treatment_pos:missing_treatment_end])) + monitored_observed_takeup_total) / monitored_cell_num_obs;
        
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
  int not_private_value_bracelet_coef[num_all_treatment_coef - num_private_value_bracelet_coef];
  
  int<lower = 0> num_bracelet_treated;
  int<lower = 1, upper = num_obs> bracelet_treated_id[num_bracelet_treated];
  int<lower = 0, upper = num_obs> strata_bracelet_sizes[num_strata];
  
  // Below references to "maps" refer to a finite set of treatment cells that has a one-to-many relationship to actual observations
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix; // Design matrix generated from treatment map
  int<lower = 1, upper = num_all_treatments> obs_treatment[num_obs]; // ID of observed treatment (from treatment map)
  
  vector<lower = 0, upper = 1>[num_obs] name_matched;
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs]; // Observation indices ordered by treatment ID (from treatment map)
 
  // Counterfactual information for finite sample analysis
  
  int<lower = 1, upper = num_obs> stratum_covar_id[num_obs]; // Observation indices ordered by stratum and endline covariate missingness
  
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
  
  int phone_owner_indices[num_obs];
 
  // Deworming outcomes
  
  int num_deworming_days;
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
  int<lower = 1, upper = num_deworming_days + 1> dewormed_day_any[num_obs];
 
  int num_dewormed; 
  int strata_dewormed_sizes[num_strata];
  int dewormed_ids[num_dewormed];
  
  matrix[num_deworming_days, num_deworming_days] hazard_day_map;
  matrix[num_deworming_days + 1, num_deworming_days] hazard_day_triangle_map;
  
  // WTP data
  
  int<lower = 0> num_wtp_obs;
  
  int<lower = 1, upper = num_obs> wtp_strata_sizes[num_strata]; 

  vector<lower = -1, upper = 1>[num_wtp_obs] gift_choice;
  vector<lower = -1, upper = 1>[num_wtp_obs] wtp_response;
  vector<lower = 0>[num_wtp_obs] wtp_offer;
  
  real<lower = 0> tau_mu_wtp_diff;
  real<lower = 0> mu_wtp_df_student_t;
  
  real<lower = 0> tau_sigma_wtp_diff;
  real<lower = 0> sigma_wtp_df_student_t;
  
  real<lower = 0> sigma_ksh_util_gamma;
  
  real<lower = 0, upper = 10> wtp_utility_df; // TODO put hyperprior on this parameter
}

transformed data {
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
  
  int num_not_private_value_bracelet_coef = num_all_treatment_coef - num_private_value_bracelet_coef;
  
  // Constants for hyperpriors 
  real scale_df = 3;
  real scale_sigma = 1;
  
  real coef_df = 7;
  real coef_sigma = 10;
  
  vector<lower = -0.5, upper = 0.5>[2] hyper_ksh_util_gamma_raw = [ 0, 0 ]'; 

  diag_treatment_map_dm[1] = rep_row_vector(0, num_all_treatment_coef);
  diag_treatment_map_dm[2:num_all_treatments] = diag_matrix(rep_vector(1, num_all_treatment_coef));
  
  diag_treatment_dm = diag_treatment_map_dm[obs_treatment];
}

parameters {
  // Modelled parameters
  
  real<lower = 0, upper = 1> hyper_baseline_day1_takeup; // Uniform[0, 1]
  
  row_vector<lower = 0>[num_strata] stratum_hazard_effect;
  vector<lower = 0>[num_deworming_days] day_hazard_effect;
  vector<lower = 0>[num_clusters] cluster_hazard_effect;
  
  // Gamma dist parameters with uniform hyperprior
  real<lower = 0, upper = 5> stratum_hazard_effect_alphabeta; 
  real<lower = 0, upper = 5> day_hazard_effect_alphabeta; 
  real<lower = 0, upper = 5> cluster_hazard_effect_alphabeta; 
  
  vector[num_not_private_value_bracelet_coef] hyper_beta_raw; // No tau for hyper parameter; coef_sigma is the SD
  vector[num_not_private_value_bracelet_coef] stratum_beta_raw[num_strata];
  vector<lower = 0>[num_not_private_value_bracelet_coef] stratum_tau_treatment;
  
  // vector[num_strata] stratum_name_matching_effect_raw;
  // real<lower = 0> tau_stratum_name_matching;
  // 
  // matrix[num_all_treatment_coef, num_strata] stratum_treatment_name_matching_interact_raw;
  // real<lower = 0> tau_stratum_treatment_name_matching_interact;
  // 
  // vector[num_census_covar_coef] hyper_census_covar_coef_raw;
  
  // WTP parameters

  real hyper_mu_wtp_diff_raw;
  vector[num_strata] mu_wtp_diff_raw;
  real<lower = 0> sigma_wtp_diff;

  // vector<lower = -0.5, upper = 0.5>[2] hyper_ksh_util_gamma_raw;

  vector[num_obs] latent_bracelet_val_diff_raw;
}

transformed parameters {
  real<lower = 0> hyper_baseline_day1_hazard = - log(hyper_baseline_day1_takeup);
  
  matrix<lower = 0>[num_deworming_days, num_strata] baseline_hazards = 
    hyper_baseline_day1_hazard * rep_matrix(stratum_hazard_effect, num_deworming_days) .* rep_matrix(day_hazard_effect, num_strata);
    
  vector[num_not_private_value_bracelet_coef] hyper_beta = hyper_beta_raw * coef_sigma;
  // vector[num_census_covar_coef] hyper_census_covar_coef = hyper_census_covar_coef_raw * coef_sigma;
  
  matrix[num_strata, num_all_treatment_coef] stratum_beta_mat;
  
  // vector[num_strata] stratum_name_matching_effect = stratum_name_matching_effect_raw * tau_stratum_name_matching;
  // matrix[num_all_treatment_coef, num_strata] stratum_treatment_name_matching_interact =
  //   stratum_treatment_name_matching_interact_raw * tau_stratum_treatment_name_matching_interact;
    
  vector[num_obs] latent_bracelet_util_diff;
    
  vector[num_obs] latent_utility = rep_vector(0, num_obs);
  // 
    
  // WTP parameters
  real hyper_mu_wtp_diff = hyper_mu_wtp_diff_raw * tau_mu_wtp_diff;
  vector[num_strata] mu_wtp_diff = (hyper_mu_wtp_diff_raw * tau_mu_wtp_diff) + mu_wtp_diff_raw * tau_mu_wtp_diff;
  
  vector<upper = 0>[num_obs] hetero_kappa = rep_vector(0, num_obs); 
  vector[num_strata] stratum_lp;
  
  {
    int stratum_pos = 1;
    int bracelet_val_stratum_pos = 1;
    int dewormed_stratum_pos = 1;
    

    for (strata_index in 1:num_strata) {
      vector[num_all_treatment_coef] stratum_beta = rep_vector(0, num_all_treatment_coef);
      vector[num_deworming_days + 1] stratum_triangle_sum_lambda = hazard_day_triangle_map * baseline_hazards[, strata_index];
      vector[num_deworming_days] stratum_dewormed_day_lambda = hazard_day_map * baseline_hazards[, strata_index];
    
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      
      int curr_bracelet_val_stratum_size = strata_bracelet_sizes[strata_index];
      int bracelet_val_stratum_end = bracelet_val_stratum_pos + curr_bracelet_val_stratum_size - 1;
      
      int curr_dewormed_stratum_size = strata_dewormed_sizes[strata_index];
      int dewormed_stratum_end = dewormed_stratum_pos + curr_dewormed_stratum_size - 1;
      
      latent_bracelet_util_diff[stratum_pos:stratum_end] =
        (mu_wtp_diff[strata_index] + (latent_bracelet_val_diff_raw[stratum_pos:stratum_end] * sigma_wtp_diff)) .* hyper_ksh_util_gamma_raw[phone_owner_indices[stratum_pos:stratum_end]];
      
      stratum_beta[not_private_value_bracelet_coef] = hyper_beta + stratum_beta_raw[strata_index] .* stratum_tau_treatment;
              
      latent_utility[stratum_pos:stratum_end] =
          // census_covar_dm[stratum_pos:stratum_end] * hyper_census_covar_coef +
          treatment_design_matrix[stratum_pos:stratum_end] * stratum_beta + 
          (treatment_design_matrix[stratum_pos:stratum_end, private_value_bracelet_coef] * stratum_beta[private_value_calendar_coef]) .*
             (1 + latent_bracelet_util_diff[stratum_pos:stratum_end]);
          // name_matched[stratum_pos:stratum_end] * stratum_name_matching_effect[strata_index] +
          // (diag_treatment_dm[stratum_pos:stratum_end] * stratum_treatment_name_matching_interact[, strata_index]) .* name_matched[stratum_pos:stratum_end];
         
      hetero_kappa[stratum_pos:stratum_end] = (- exp(latent_utility[stratum_pos:stratum_end])) .* cluster_hazard_effect[cluster_id[stratum_pos:stratum_end]];
    
      stratum_lp[strata_index] = 
        sum(hetero_kappa[stratum_pos:stratum_end] .* stratum_triangle_sum_lambda[dewormed_day_any[stratum_pos:stratum_end]]) + 
        sum(log1m_exp(hetero_kappa[dewormed_ids[dewormed_stratum_pos:dewormed_stratum_end]] .*
                        stratum_dewormed_day_lambda[dewormed_day_any[dewormed_ids[dewormed_stratum_pos:dewormed_stratum_end]]]));
          
      stratum_beta_mat[strata_index] = stratum_beta';
      
      stratum_pos = stratum_end + 1;
      bracelet_val_stratum_pos = bracelet_val_stratum_end + 1;
      dewormed_stratum_pos = dewormed_stratum_end + 1;
    }
  }
}

model {
  // WTP sampling
  
  hyper_mu_wtp_diff_raw ~ student_t(mu_wtp_df_student_t, 0, 1);
  mu_wtp_diff_raw ~ student_t(mu_wtp_df_student_t, 0, 1);
  sigma_wtp_diff ~ student_t(sigma_wtp_df_student_t, 0, tau_sigma_wtp_diff);

  {
    int wtp_stratum_pos = 1;

    for (strata_index in 1:num_strata) {
      int curr_wtp_stratum_size = wtp_strata_sizes[strata_index];
      int wtp_stratum_end = wtp_stratum_pos + curr_wtp_stratum_size - 1;

      for (wtp_obs_index in wtp_stratum_pos:wtp_stratum_end) {
        if (wtp_response[wtp_obs_index] == -1) {
          if (gift_choice[wtp_obs_index] == -1) {
            target += student_t_lcdf(- wtp_offer[wtp_obs_index] | wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff);
          } else {
            target += student_t_lccdf(wtp_offer[wtp_obs_index] | wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff);
          }
        } else {
          target += log(gift_choice[wtp_obs_index] * (student_t_cdf(gift_choice[wtp_obs_index] * wtp_offer[wtp_obs_index], wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff) -
                                                      student_t_cdf(0, wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff)));
        }
      }

      wtp_stratum_pos = wtp_stratum_end + 1;
    }
  }

  latent_bracelet_val_diff_raw ~ student_t(wtp_utility_df, 0, 1);
  
  // Take-up sampling
  
  stratum_hazard_effect ~ gamma(stratum_hazard_effect_alphabeta, stratum_hazard_effect_alphabeta);
  day_hazard_effect ~ gamma(day_hazard_effect_alphabeta, day_hazard_effect_alphabeta);
  cluster_hazard_effect ~ gamma(cluster_hazard_effect_alphabeta, cluster_hazard_effect_alphabeta);
  
  hyper_beta_raw ~ student_t(coef_df, 0, 1); 
  stratum_beta_raw ~ multi_student_t(coef_df, 
                                     rep_vector(0.0, num_not_private_value_bracelet_coef), 
                                     diag_matrix(rep_vector(1, num_not_private_value_bracelet_coef)));
  stratum_tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  
  // stratum_name_matching_effect_raw ~ student_t(coef_df, 0, 1);
  // tau_stratum_name_matching ~ student_t(scale_df, 0, scale_sigma);
  // to_vector(stratum_treatment_name_matching_interact_raw) ~ student_t(coef_df, 0, 1);
  // tau_stratum_treatment_name_matching_interact ~ student_t(scale_df, 0, scale_sigma);

  // hyper_census_covar_coef_raw ~ student_t(coef_df, 0, 1);

  target += sum(stratum_lp);
}

generated quantities {
  /*
  matrix<lower = 0, upper = 1>[num_non_phone_owner_treatments, 2] non_phone_takeup; 
  matrix<lower = 0, upper = 1>[num_phone_owner_treatments, 2] phone_takeup; 
  vector<lower = -1, upper = 1>[num_non_phone_owner_ate_pairs] non_phone_takeup_ate;
  vector<lower = -1, upper = 1>[num_phone_owner_ate_pairs] phone_takeup_ate;
  vector[num_non_phone_owner_ate_pairs] non_phone_takeup_ate_percent;
  vector[num_phone_owner_ate_pairs] phone_takeup_ate_percent;
  
  {
    vector[num_missing_non_phone_owner_obs_ids] missing_non_phone_bracelet_util_diff = latent_bracelet_util_diff[missing_non_phone_owner_obs_ids];
    vector[num_missing_phone_owner_obs_ids] missing_phone_bracelet_util_diff = latent_bracelet_util_diff[missing_phone_owner_obs_ids];
    
    non_phone_takeup = treatment_cell_takeup_rng(non_phone_owner_treatments,
                                             missing_non_phone_owner_obs_ids,
                                             non_phone_missing_treatment_stratum_id,
                                             non_phone_missing_treatment_cluster_id,
                                             missing_non_phone_bracelet_util_diff,
                                             private_value_calendar_coef,
                                             private_value_bracelet_coef,
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
                                             stratum_treatment_name_matching_interact,
                                             stratum_beta_mat,
                                             to_vector(observed_non_phone_dewormed_any));
                                             
    phone_takeup = treatment_cell_takeup_rng(phone_owner_treatments,
                                             missing_phone_owner_obs_ids,
                                             phone_missing_treatment_stratum_id,
                                             phone_missing_treatment_cluster_id,
                                             missing_phone_bracelet_util_diff,
                                             private_value_calendar_coef,
                                             private_value_bracelet_coef,
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
                                             stratum_treatment_name_matching_interact,
                                             stratum_beta_mat,
                                             to_vector(observed_phone_dewormed_any));
  }
                                           
  non_phone_takeup_ate = non_phone_takeup[non_phone_owner_ate_pairs[, 1], 1] - non_phone_takeup[non_phone_owner_ate_pairs[, 2], 1];
  phone_takeup_ate = phone_takeup[phone_owner_ate_pairs[, 1], 1] - phone_takeup[phone_owner_ate_pairs[, 2], 1];
  non_phone_takeup_ate_percent = non_phone_takeup_ate ./ non_phone_takeup[non_phone_owner_ate_pairs[, 2], 2];
  phone_takeup_ate_percent = phone_takeup_ate ./ phone_takeup[phone_owner_ate_pairs[, 2], 2];*/
}
