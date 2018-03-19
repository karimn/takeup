data {
  int<lower = 0> num_obs;
  int<lower = 0, upper = num_obs> num_with_links_obs;
  int<lower = 0, upper = num_with_links_obs> num_thinks_other_knows;
  int<lower = 1, upper = num_obs> num_treatments;
  int<lower = 1> num_treatment_coef;
  
  int<lower = 1, upper = num_treatments> control_treatment_map[num_treatments];
  
  int<lower = 1, upper = num_treatments> treatment_ids[num_obs];
  int<lower = 1, upper = num_treatments> control_treatment_ids[num_obs];
  
  int<lower = 0, upper = num_obs> treatment_obs_sizes[num_treatments];
  int<lower = 0, upper = num_obs> treatment_missing_sizes[num_treatments];
  
  int<lower = 1, upper = num_obs> treatment_obs_ids[sum(treatment_obs_sizes)];
  int<lower = 1, upper = num_obs> treatment_missing_ids[sum(treatment_missing_sizes)];
  
  int<lower = 0, upper = num_obs> know_treatment_obs_sizes[num_treatments];
  int<lower = 0, upper = num_obs> know_treatment_missing_sizes[num_treatments];
  
  int<lower = 1, upper = num_obs> know_treatment_obs_ids[sum(know_treatment_obs_sizes)];
  int<lower = 1, upper = num_obs> know_treatment_missing_ids[sum(know_treatment_missing_sizes)];
  
  int<lower = 0, upper = num_obs> other_know_treatment_obs_sizes[num_treatments];
  int<lower = 0, upper = num_obs> other_know_treatment_missing_sizes[num_treatments];
  
  int<lower = 1, upper = num_obs> other_know_treatment_obs_ids[sum(other_know_treatment_obs_sizes)];
  int<lower = 1, upper = num_obs> other_know_treatment_missing_ids[sum(other_know_treatment_missing_sizes)];
  
  int<lower = 1, upper = num_obs> with_links_ids[num_with_links_obs];
  int<lower = 1, upper = num_obs> thinks_other_knows_ids[num_thinks_other_knows];
  
  int<lower = 0, upper = 10> obs_know_person[num_obs];
  int<lower = 0, upper = 10> thinks_other_knows[num_obs];
  int<lower = 0, upper = 10> thinks_other_knows_yes[num_obs];
  
  matrix[num_treatments, num_treatment_coef] treatment_map_design_matrix;
  
  int<lower = 1> village_size[num_obs];
}

transformed data {
  int obs_per_person = 10;
  int num_second_order_dm_col = num_treatment_coef * 4 + 3;
  
  matrix[num_obs, num_treatment_coef] treatment_design_matrix = treatment_map_design_matrix[treatment_ids];
  matrix[num_obs, num_treatment_coef] control_treatment_design_matrix = treatment_map_design_matrix[control_treatment_ids];
 
  vector<lower = 1>[num_obs] vec_village_size = to_vector(village_size); 
  vector[num_obs] stdz_village_size = (vec_village_size - mean(vec_village_size)) / sd(vec_village_size);
  
  int with_links_obs_index[num_obs] = rep_array(-1, num_obs);
  int thinks_other_knows_obs_index[num_obs] = rep_array(-1, num_obs);
  
  for (with_links_index in 1:num_with_links_obs) {
    with_links_obs_index[with_links_ids[with_links_index]] = with_links_index;
  }
  
  for (thinks_other_knows_index in 1:num_thinks_other_knows) {
    thinks_other_knows_obs_index[thinks_other_knows_ids[thinks_other_knows_index]] = thinks_other_knows_index;
  }
}

parameters {
  vector[num_treatment_coef] beta_know_person;
  vector[num_obs] intercept_know_person;
  real<lower = 0> intercept_sd_know_person;
  
  vector[num_second_order_dm_col] beta_thinks_other_knows;
  vector[num_with_links_obs] intercept_thinks_other_knows;
  real<lower = 0> intercept_sd_thinks_other_knows;
  
  vector[num_second_order_dm_col] beta_thinks_other_knows_yes;
  vector[num_thinks_other_knows] intercept_thinks_other_knows_yes;
  real<lower = 0> intercept_sd_thinks_other_knows_yes;
}

transformed parameters {
  vector[num_obs] latent_var_know_person = intercept_know_person + treatment_design_matrix * beta_know_person;
  vector[num_obs] social_degree = vec_village_size .* inv_logit(latent_var_know_person);
  real mean_degree = mean(social_degree);
  real sd_degree = sd(social_degree);
  vector[num_obs] stdz_degree = (social_degree - mean_degree) / sd_degree;
}

model {
  intercept_sd_know_person ~ normal(0, 5);
  intercept_know_person ~ normal(0, intercept_sd_know_person);
  beta_know_person ~ normal(0, 1);
  
  intercept_sd_thinks_other_knows ~ normal(0, 5);
  intercept_thinks_other_knows ~ normal(0, intercept_sd_thinks_other_knows);
  beta_thinks_other_knows ~ normal(0, 1);
  
  intercept_sd_thinks_other_knows_yes ~ normal(0, 5);
  intercept_thinks_other_knows_yes ~ normal(0, intercept_sd_thinks_other_knows_yes);
  beta_thinks_other_knows_yes ~ normal(0, 1);
  
  {
    matrix[num_obs, num_second_order_dm_col] treatment_social_design_matrix;
    vector[num_with_links_obs] latent_var_thinks_other_knows; 
    vector[num_thinks_other_knows] latent_var_thinks_other_knows_yes; 
    
    obs_know_person ~ binomial_logit(obs_per_person, latent_var_know_person);
    
    treatment_social_design_matrix[, 1:num_treatment_coef] = treatment_design_matrix;
    treatment_social_design_matrix[, (num_treatment_coef + 1):(2 * num_treatment_coef)] = 
      treatment_design_matrix .* rep_matrix(stdz_degree, num_treatment_coef);
    treatment_social_design_matrix[, (2 * num_treatment_coef + 1):(3 * num_treatment_coef)] = 
      treatment_design_matrix .* rep_matrix(stdz_village_size, num_treatment_coef);
      
    treatment_social_design_matrix[, (4 * num_treatment_coef) + 1] = stdz_degree;
    treatment_social_design_matrix[, (4 * num_treatment_coef) + 2] = stdz_village_size;
    treatment_social_design_matrix[, (4 * num_treatment_coef) + 3] = stdz_degree .* stdz_village_size;
    
    treatment_social_design_matrix[, (3 * num_treatment_coef + 1):(4 * num_treatment_coef)] = 
      treatment_design_matrix .* rep_matrix(stdz_degree .* stdz_village_size, num_treatment_coef);
    
    latent_var_thinks_other_knows = intercept_thinks_other_knows + treatment_social_design_matrix[with_links_ids] * beta_thinks_other_knows;
    latent_var_thinks_other_knows_yes = intercept_thinks_other_knows_yes + treatment_social_design_matrix[thinks_other_knows_ids] * beta_thinks_other_knows_yes;
    
    thinks_other_knows[with_links_ids] ~ binomial_logit(obs_know_person[with_links_ids], latent_var_thinks_other_knows); 
    thinks_other_knows_yes[thinks_other_knows_ids] ~ binomial_logit(thinks_other_knows[thinks_other_knows_ids], latent_var_thinks_other_knows_yes); 
  }
}

generated quantities { 
  vector<lower = 0, upper = 1>[num_treatments] mean_know_person;
  vector<lower = 0, upper = 1>[num_treatments] mean_thinks_other_know;
  vector<lower = 0, upper = 1>[num_treatments] mean_thinks_other_know_yes;
  
  vector<lower = -1, upper = 1>[num_treatments] ate_know_person;
  vector<lower = -1, upper = 1>[num_treatments] ate_thinks_other_know;
  vector<lower = -1, upper = 1>[num_treatments] ate_thinks_other_know_yes;
  
  {
    int treatment_obs_pos = 1;
    int treatment_missing_pos = 1;
    int know_treatment_obs_pos = 1;
    int know_treatment_missing_pos = 1;
    int other_know_treatment_obs_pos = 1;
    int other_know_treatment_missing_pos = 1;
    
    for (treatment_index in 1:num_treatments) {
      int curr_treatment_obs_size = treatment_obs_sizes[treatment_index];
      int curr_treatment_missing_size = treatment_missing_sizes[treatment_index];
      int treatment_obs_end = treatment_obs_pos + curr_treatment_obs_size - 1; 
      int treatment_missing_end = treatment_missing_pos + curr_treatment_missing_size - 1; 
      
      int know_curr_treatment_obs_size = know_treatment_obs_sizes[treatment_index];
      int know_curr_treatment_missing_size = know_treatment_missing_sizes[treatment_index];
      int know_treatment_obs_end = know_treatment_obs_pos + know_curr_treatment_obs_size - 1; 
      int know_treatment_missing_end = know_treatment_missing_pos + know_curr_treatment_missing_size - 1; 
      
      int other_know_curr_treatment_obs_size = other_know_treatment_obs_sizes[treatment_index];
      int other_know_curr_treatment_missing_size = other_know_treatment_missing_sizes[treatment_index];
      int other_know_treatment_obs_end = other_know_treatment_obs_pos + other_know_curr_treatment_obs_size - 1; 
      int other_know_treatment_missing_end = other_know_treatment_missing_pos + other_know_curr_treatment_missing_size - 1; 
      
      real mean_obs_know_person = mean(to_vector(obs_know_person[treatment_obs_ids[treatment_obs_pos:treatment_obs_end]]) ./ 10);
      vector[curr_treatment_missing_size] missing_know_person;
      
      matrix[num_obs, num_second_order_dm_col] treatment_social_design_matrix; 
      
      real mean_obs_thinks_other_know = mean(to_vector(thinks_other_knows[know_treatment_obs_ids[know_treatment_obs_pos:know_treatment_obs_end]]) 
                                             ./ to_vector(obs_know_person[know_treatment_obs_ids[know_treatment_obs_pos:know_treatment_obs_end]]));
      vector[know_curr_treatment_missing_size] missing_thinks_other_know;
      
      real mean_obs_thinks_other_know_yes = mean(to_vector(thinks_other_knows_yes[other_know_treatment_obs_ids[other_know_treatment_obs_pos:other_know_treatment_obs_end]]) 
                                             ./ to_vector(thinks_other_knows[other_know_treatment_obs_ids[other_know_treatment_obs_pos:other_know_treatment_obs_end]]));
      vector[other_know_curr_treatment_missing_size] missing_thinks_other_know_yes;
     
      matrix[num_obs, num_treatment_coef] curr_dm = rep_matrix(treatment_design_matrix[treatment_index], num_obs);
      
      treatment_social_design_matrix[, 1:num_treatment_coef] = curr_dm;
      treatment_social_design_matrix[, (num_treatment_coef + 1):(2 * num_treatment_coef)] = 
        curr_dm .* rep_matrix(stdz_degree, num_treatment_coef);
      treatment_social_design_matrix[, (2 * num_treatment_coef + 1):(3 * num_treatment_coef)] = 
        curr_dm .* rep_matrix(stdz_village_size, num_treatment_coef);
        
      treatment_social_design_matrix[, (4 * num_treatment_coef) + 1] = stdz_degree;
      treatment_social_design_matrix[, (4 * num_treatment_coef) + 2] = stdz_village_size;
      treatment_social_design_matrix[, (4 * num_treatment_coef) + 3] = stdz_degree .* stdz_village_size;
      
      treatment_social_design_matrix[, (3 * num_treatment_coef + 1):(4 * num_treatment_coef)] = 
        curr_dm .* rep_matrix(stdz_degree .* stdz_village_size, num_treatment_coef);
      
      for (missing_index in 1:curr_treatment_missing_size) {
        int curr_missing_id = treatment_missing_ids[treatment_missing_pos + missing_index - 1];
       
        missing_know_person[missing_index] = 
          binomial_rng(obs_per_person, 
                       inv_logit(intercept_know_person[curr_missing_id] + treatment_map_design_matrix[treatment_index] * beta_know_person));
      }
      
      for (know_missing_index in 1:know_curr_treatment_missing_size) {
        int know_curr_missing_id = know_treatment_missing_ids[know_treatment_missing_pos + know_missing_index - 1];
       
        missing_thinks_other_know[know_missing_index] = 
          binomial_rng(obs_know_person[know_curr_missing_id], 
                       inv_logit(intercept_thinks_other_knows[with_links_obs_index[know_curr_missing_id]] 
                       + treatment_social_design_matrix[know_curr_missing_id] * beta_thinks_other_knows));
      }
      
      for (other_know_missing_index in 1:other_know_curr_treatment_missing_size) {
        int other_know_curr_missing_id = other_know_treatment_missing_ids[other_know_treatment_missing_pos + other_know_missing_index - 1];
       
        missing_thinks_other_know_yes[other_know_missing_index] = 
          binomial_rng(thinks_other_knows[other_know_curr_missing_id], 
                       inv_logit(intercept_thinks_other_knows_yes[thinks_other_knows_obs_index[other_know_curr_missing_id]] 
                       + treatment_social_design_matrix[other_know_curr_missing_id] * beta_thinks_other_knows_yes));
      }
      
      mean_know_person[treatment_index] = 
        (curr_treatment_obs_size * mean_obs_know_person + curr_treatment_missing_size * mean(missing_know_person ./ 10)) 
        / (curr_treatment_obs_size + curr_treatment_missing_size);
        
      mean_thinks_other_know[treatment_index] = 
        (know_curr_treatment_obs_size * mean_obs_thinks_other_know 
         + (know_curr_treatment_missing_size 
            * mean(missing_thinks_other_know ./ to_vector(obs_know_person[know_treatment_missing_ids[know_treatment_missing_pos:know_treatment_missing_end]])))) 
        / (know_curr_treatment_obs_size + know_curr_treatment_missing_size);
        
      mean_thinks_other_know_yes[treatment_index] = 
        (other_know_curr_treatment_obs_size * mean_obs_thinks_other_know_yes
         + (other_know_curr_treatment_missing_size 
            * mean(missing_thinks_other_know_yes ./ to_vector(thinks_other_knows[other_know_treatment_missing_ids[other_know_treatment_missing_pos:other_know_treatment_missing_end]]))))
        / (other_know_curr_treatment_obs_size + other_know_curr_treatment_missing_size);
      
      treatment_obs_pos = treatment_obs_end + 1;
      treatment_missing_pos = treatment_missing_end + 1;
      know_treatment_obs_pos = know_treatment_obs_end + 1;
      know_treatment_missing_pos = know_treatment_missing_end + 1;
      other_know_treatment_obs_pos = other_know_treatment_obs_end + 1;
      other_know_treatment_missing_pos = other_know_treatment_missing_end + 1;
    }
  }
  
  ate_know_person = mean_know_person - mean_know_person[control_treatment_map];
  ate_thinks_other_know = mean_thinks_other_know - mean_thinks_other_know[control_treatment_map];
  ate_thinks_other_know_yes = mean_thinks_other_know_yes - mean_thinks_other_know_yes[control_treatment_map];
}
