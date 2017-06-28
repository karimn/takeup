functions {
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
  
  int<lower = 0> num_bracelet_treated;
  int<lower = 1, upper = num_obs> bracelet_treated_id[num_bracelet_treated];
  int<lower = 0, upper = num_obs> strata_bracelet_sizes[num_strata];
  
  // int<lower = 0> num_calendar_treated;
  // int<lower = 1, upper = num_obs> calendar_treated_id[num_calendar_treated];
  // int<lower = 0, upper = num_obs> strata_calendar_sizes[num_strata];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1, upper = num_obs> strata_sizes[num_strata]; 
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
  
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
  
  real<lower = 0, upper = 10> wtp_utility_df; // TODO put hyperprior on this parameter
  
  // Constants for hyperpriors 
  real scale_df;
  real scale_sigma;
  
  real coef_df;
  real coef_sigma;
  real bracelet_val_coef_sigma;
}

transformed data {
  matrix[num_obs, num_all_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[obs_treatment];
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
  vector<lower = 0>[num_all_treatment_coef] tau_treatment;
  
  // WTP parameters
  
  real hyper_mu_wtp_diff_raw;
  vector[num_strata] mu_wtp_diff_raw;
  real<lower = 0> sigma_wtp_diff;
  
  vector[num_bracelet_treated] latent_bracelet_val_diff_raw;
  // real mu_calendar_util_raw;
  real ksh_util_gamma_raw;
  // real<lower = 0> sigma_calendar_util;
  // vector[num_calendar_treated] latent_calendar_util_raw;

  // real<upper = 0> bracelet_val_coef_raw;
  // vector[num_strata] bracelet_val_coef_raw;
}

transformed parameters {
  real hyper_intercept = inv_logit(hyper_baseline_takeup);
  
  // real<upper = 0> bracelet_val_coef = bracelet_val_coef_raw * bracelet_val_coef_sigma; // No hyperparam; not sure best was to pick a negative only param
  // real<upper = 0> bracelet_val_coef = bracelet_val_coef_raw * bracelet_val_coef_sigma; // No hyperparam; not sure best was to pick a negative only param
  // vector[num_strata] bracelet_val_coef = hyper_bracelet_val_coef + bracelet_val_coef_raw * coef_sigma;
  real ksh_util_gamma = ksh_util_gamma_raw * coef_sigma;
  
  vector[num_all_treatment_coef] hyper_beta = hyper_beta_raw * coef_sigma;
  
  vector[num_strata] stratum_intercept = hyper_intercept + stratum_intercept_raw * tau_stratum_intercept;
  vector[num_clusters] cluster_effects = cluster_effects_raw * tau_cluster_effect;
  
  vector[num_all_treatment_coef] beta[num_strata];
  vector[num_obs] latent_utility = rep_vector(0, num_obs);
  
  // WTP parameters
  real hyper_mu_wtp_diff = hyper_mu_wtp_diff_raw * tau_mu_wtp_diff;
  vector[num_strata] mu_wtp_diff = hyper_mu_wtp_diff + mu_wtp_diff_raw * tau_mu_wtp_diff;
  
  // real mu_calendar_util = mu_calendar_util_raw * coef_sigma;
  // vector[num_calendar_treated] latent_calendar_util;
  
  vector[num_bracelet_treated] latent_bracelet_val_diff;
  
  {
    int stratum_pos = 1;
    int bracelet_val_stratum_pos = 1;
    int calendar_util_stratum_pos = 1;
    vector[num_obs] full_bracelet_util_diff = rep_vector(0, num_obs);
    // vector[num_obs] full_calendar_util = rep_vector(0, num_obs);
    

    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = stratum_pos + curr_stratum_size - 1;
      int curr_bracelet_val_stratum_size = strata_bracelet_sizes[strata_index];
      int bracelet_val_stratum_end = bracelet_val_stratum_pos + curr_bracelet_val_stratum_size - 1;
      // int curr_calendar_util_stratum_size = strata_calendar_sizes[strata_index];
      // int calendar_util_stratum_end = calendar_util_stratum_pos + curr_calendar_util_stratum_size - 1;
      
      // latent_calendar_util[calendar_util_stratum_pos:calendar_util_stratum_end] =
      //   mu_calendar_util + latent_calendar_util_raw[calendar_util_stratum_pos:calendar_util_stratum_end] * sigma_calendar_util;
      //   
      // for (i in calendar_util_stratum_pos:calendar_util_stratum_end) {
      //   full_calendar_util[calendar_treated_id[i]] = step(latent_calendar_util[i]);
      // }
        
      // full_bracelet_val_diff[bracelet_treated_id[bracelet_val_stratum_pos:bracelet_val_stratum_end]] =
      //   latent_bracelet_val_diff[bracelet_val_stratum_pos:bracelet_val_stratum_end];
      
      latent_bracelet_val_diff[bracelet_val_stratum_pos:bracelet_val_stratum_end] = 
        mu_wtp_diff[strata_index] + latent_bracelet_val_diff_raw[bracelet_val_stratum_pos:bracelet_val_stratum_end] * sigma_wtp_diff;
        
      full_bracelet_util_diff[bracelet_treated_id[bracelet_val_stratum_pos:bracelet_val_stratum_end]] = 
         - (latent_bracelet_val_diff[bracelet_val_stratum_pos:bracelet_val_stratum_end] * ksh_util_gamma);

      beta[strata_index] = (hyper_beta + beta_raw[strata_index] .* tau_treatment) .* [ 1, ksh_util_gamma ]';

      latent_utility[stratum_pos:stratum_end] =
          stratum_intercept[strata_index] +
          cluster_effects[cluster_id[stratum_pos:stratum_end]] +
          treatment_design_matrix[stratum_pos:stratum_end] * beta[strata_index] +
          full_bracelet_util_diff[stratum_pos:stratum_end];

      stratum_pos = stratum_end + 1;
      bracelet_val_stratum_pos = bracelet_val_stratum_end + 1;
    }
  }
}

model {
  stratum_intercept_raw ~ student_t(coef_df, 0, 1); 
  tau_stratum_intercept ~ student_t(scale_df, 0, scale_sigma);
  
  cluster_effects_raw ~ student_t(coef_df, 0, 1); 
  tau_cluster_effect ~ student_t(scale_df, 0, scale_sigma);
  
  hyper_beta_raw ~ student_t(coef_df, 0, 1);
  beta_raw ~ multi_student_t(coef_df, rep_vector(0.0, num_all_treatment_coef), diag_matrix(rep_vector(1, num_all_treatment_coef)));
  tau_treatment ~ student_t(scale_df, 0, scale_sigma);
  
  // hyper_bracelet_val_coef_raw ~ student_t(coef_df, 0, 1);
  // bracelet_val_coef_raw ~ student_t(coef_df, 0, 1);
  
  // WTP sampling
  
  hyper_mu_wtp_diff_raw ~ student_t(mu_wtp_df_student_t, 0, 1);
  mu_wtp_diff_raw ~ student_t(mu_wtp_df_student_t, 0, 1);
  sigma_wtp_diff ~ student_t(sigma_wtp_df_student_t, 0, tau_sigma_wtp_diff); 
  
  ksh_util_gamma_raw ~ student_t(coef_df, 0, 1);
  
  latent_bracelet_val_diff_raw ~ student_t(wtp_utility_df, 0, 1);
  
  {
    int wtp_stratum_pos = 1;
    // int bracelet_val_stratum_pos = 1;

    for (strata_index in 1:num_strata) {
      int curr_wtp_stratum_size = wtp_strata_sizes[strata_index];
      int wtp_stratum_end = wtp_stratum_pos + curr_wtp_stratum_size - 1;
      // int curr_bracelet_val_stratum_size = strata_bracelet_sizes[strata_index];
      // int bracelet_val_stratum_end = bracelet_val_stratum_pos + curr_bracelet_val_stratum_size - 1;

      for (i in wtp_stratum_pos:wtp_stratum_end) {
        if (wtp_response[i] == -1) {
          if (gift_choice[i] == -1) {
            target += student_t_lcdf(- wtp_offer[i] | wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff);
          } else {
            target += student_t_lccdf(wtp_offer[i] | wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff);
          }
        } else {
          target += log(gift_choice[i] * (student_t_cdf(gift_choice[i] * wtp_offer[i], wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff)
                    - student_t_cdf(0, wtp_utility_df, mu_wtp_diff[strata_index], sigma_wtp_diff)));
        }
      }


      wtp_stratum_pos = wtp_stratum_end + 1;
      // bracelet_val_stratum_pos = bracelet_val_stratum_end + 1;
    }
  }
  
  // mu_calendar_util_raw ~ student_t(coef_df, 0, 1);
  // sigma_calendar_util ~ uniform(0, gamma_util * tau_sigma_wtp_diff);
  // latent_calendar_util_raw ~ student_t(coef_df, 0, 1);
  
  dewormed_any ~ bernoulli_logit(latent_utility);
}

generated quantities {
}
