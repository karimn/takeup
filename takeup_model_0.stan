functions {
  real[] takeup_proportion_rng(int[] eval_treatment_ids, // the treatment IDs to evaluate proportions for  
                               matrix treatment_map_design_matrix,
                               matrix name_match_interact_map_design_matrix,
                               int[] dewormed_any,
                               vector[] beta,
                               vector stratum_intercept,
                               vector cluster_effects,
                               int[] treatment_id,
                               int[] missing_treatment_id,
                               int[] treatment_sizes,
                               int[] missing_treatment_sizes,
                               int[] missing_stratum_id_individ,
                               int[] missing_cluster_id) { 
    int num_eval_treatment_prop = size(eval_treatment_ids);
    int total_num_coef = cols(treatment_map_design_matrix) + cols(name_match_interact_map_design_matrix);
    real takeup_proportion[num_eval_treatment_prop]; 
    int num_strata = num_elements(stratum_intercept);
    real missing_link_model[num_strata, num_eval_treatment_prop];
    
    for (strata_index in 1:num_strata) {
      matrix[num_eval_treatment_prop, total_num_coef] all_treat_map_design_matrix = 
        append_col(treatment_map_design_matrix[eval_treatment_ids], name_match_interact_map_design_matrix[eval_treatment_ids]); 
        
      missing_link_model[strata_index, ] = to_array_1d(all_treat_map_design_matrix * beta[strata_index] + stratum_intercept[strata_index]);
    }
    
    for (treatment_index_index in 1:num_eval_treatment_prop) {
      int treatment_index = eval_treatment_ids[treatment_index_index]; 
      int treatment_pos;  
      int missing_treatment_pos; 
      int curr_treatment_size = treatment_sizes[treatment_index];
      int curr_missing_treatment_size = missing_treatment_sizes[treatment_index];
      int treatment_end;
      int missing_treatment_end;
      int missing_stratum_id[curr_missing_treatment_size];
      int treatment_missing_cluster_id[curr_missing_treatment_size];
      real missing_treatment_link[curr_missing_treatment_size];
      real missing_dewormed[curr_missing_treatment_size];
      int curr_missing_id[curr_missing_treatment_size];
      
      if (treatment_index > 1) {
        treatment_pos = sum(treatment_sizes[1:(treatment_index - 1)]) + 1;
        missing_treatment_pos = sum(missing_treatment_sizes[1:(treatment_index - 1)]) + 1;
      } else {
        treatment_pos = 1;
        missing_treatment_pos = 1;
      }
      
      treatment_end = treatment_pos + curr_treatment_size - 1;
      missing_treatment_end = missing_treatment_pos + curr_missing_treatment_size - 1;
      curr_missing_id = missing_treatment_id[missing_treatment_pos:missing_treatment_end];
      missing_stratum_id = missing_stratum_id_individ[curr_missing_id];
      treatment_missing_cluster_id = missing_cluster_id[curr_missing_id];
      missing_treatment_link = missing_link_model[missing_stratum_id, treatment_index_index];
      
      for(i in 1:curr_missing_treatment_size) {
        missing_dewormed[i] = bernoulli_rng(Phi(missing_treatment_link[i] + cluster_effects[treatment_missing_cluster_id[i]]));
      }
      
      takeup_proportion[treatment_index_index] = 
        (sum(missing_dewormed) + sum(dewormed_any[treatment_id[treatment_pos:treatment_end]])) /
        (curr_missing_treatment_size + curr_treatment_size);
    }
    
    return(takeup_proportion);
  }
}

data {
  int<lower = 0> num_obs;
  int<lower = 1> num_missing;
  int<lower = 1> num_all_treatments;
  int<lower = 1> num_all_treatment_coef;
  int<lower = 1> num_name_match_interact_coef;
  int<lower = 1> num_eval_treatment_prop;
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  matrix[num_all_treatments, num_all_treatment_coef] treatment_map_design_matrix;
  matrix[num_all_treatments, num_name_match_interact_coef] name_match_interact_map_design_matrix;
  int<lower = 1, upper = num_all_treatments> obs_treatment[num_obs];
  int<lower = 1, upper = num_all_treatments> eval_treatment_prop_id[num_eval_treatment_prop];
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_obs];
  int<lower = 1> strata_sizes[num_strata]; 
  
  int<lower = 1, upper = num_obs> treatment_id[num_obs];
  int<lower = 1, upper = num_missing> missing_treatment_id[num_missing];
  int<lower = 1> treatment_sizes[num_all_treatments];
  int<lower = 1> missing_treatment_sizes[num_all_treatments];
 
  // Binary deworming outcome 
  
  int<lower = 0, upper = 1> dewormed_any[num_obs];
}

transformed data {
  matrix[num_obs, num_all_treatment_coef + num_name_match_interact_coef] all_treatment_design_matrix = 
    append_col(treatment_map_design_matrix[obs_treatment], name_match_interact_map_design_matrix[obs_treatment]);

  // int<lower = 1, upper = num_strata> stratum_id_individ[num_obs] = stratum_id[cluster_id];
  int<lower = 1, upper = num_strata> missing_stratum_id[num_missing] = stratum_id[missing_treatment_id];
  int<lower = 1, upper = num_clusters> missing_cluster_id[num_missing] = cluster_id[missing_treatment_id];
  
  // Unmodeled parameters
  
  vector[num_all_treatment_coef + num_name_match_interact_coef] mu = rep_vector(0, num_all_treatment_coef + num_name_match_interact_coef);
  vector[num_all_treatment_coef] tau_treatment = rep_vector(1, num_all_treatment_coef);
  vector[num_name_match_interact_coef] tau_name_match_interact = rep_vector(sqrt(0.25), num_name_match_interact_coef);
  vector[num_all_treatment_coef + num_name_match_interact_coef] tau = append_row(tau_treatment, tau_name_match_interact);
  cov_matrix[num_all_treatment_coef + num_name_match_interact_coef] Sigma = diag_matrix(tau);
}

parameters {
  real<lower=-3, upper=3> mu_strata; // uniformly distributed prior

  vector[num_strata] stratum_intercept;
  vector[num_clusters] cluster_effects;
  vector[num_all_treatment_coef + num_name_match_interact_coef] beta[num_strata];
  
  // Scale hyperparameters for betas 
  
  // corr_matrix[num_treatments] Omega_beta;
  // vector<lower = 0>[num_all_treatment_coef] tau_beta;
  // real<lower = 0> tau_stratum_effect;
  // real<lower = 0> tau_cluster_effect;
}

transformed parameters {
  // cov_matrix[num_all_treatment_coef] Sigma_beta = diag_matrix(tau_beta);
}

model {
  //Omega_beta ~ lkj_corr(50);
  //tau_beta ~ gamma(2, 1/10);
  
  // tau_beta ~ normal(0, 100); // cauchy(0, 20); // This is a very weakly informative prior
  // tau_stratum_effect ~ normal(0, 10); // cauchy(0, 2.5);
  // tau_cluster_effect ~ normal(0, 10); // cauchy(0, 2.5);
  
  stratum_intercept ~ normal(mu_strata, 1); # tau_stratum_effect); 
  cluster_effects ~ normal(0, 1); # tau_cluster_effect);
  
  beta ~ multi_normal(mu, Sigma); // For now assuming no correlation between effects
 
  {
    int strata_pos = 1;
    vector[num_obs] link_model;  
  
    // Looping over strata in the data; Stan doesn't support ragged arrays, so can't split the data into arrays.
    // Otherwise, will have to loop by row, which will slow things down dramatically.
    for (strata_index in 1:num_strata) {
      int curr_stratum_size = strata_sizes[strata_index];
      int stratum_end = strata_pos + curr_stratum_size - 1;
      
      link_model[strata_pos:stratum_end] = all_treatment_design_matrix[strata_pos:stratum_end] * beta[strata_index];
        
      strata_pos = strata_pos + curr_stratum_size;
    }
    
    link_model = link_model + stratum_intercept[stratum_id] + cluster_effects[cluster_id]; 
  
    dewormed_any ~ bernoulli(Phi(link_model)); // Probit
    //dewormed_any ~ bernoulli_logit(link_model);
  } 
}

generated quantities {
  real<lower = 0, upper = 1> takeup_proportion[num_eval_treatment_prop] =
    takeup_proportion_rng(
      eval_treatment_prop_id, treatment_map_design_matrix, name_match_interact_map_design_matrix, dewormed_any, beta, stratum_intercept, cluster_effects,
      treatment_id, missing_treatment_id, treatment_sizes, missing_treatment_sizes, missing_stratum_id, missing_cluster_id
    );
}
