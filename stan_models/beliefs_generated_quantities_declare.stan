matrix<lower = 0, upper = 1>[num_beliefs_obs, num_dist_group_treatments] obs_prob_1ord; // = inv_logit(centered_obs_beta_1ord * beliefs_treatment_map_design_matrix');
vector<lower = -1, upper = 1>[num_dist_group_treatments] prob_1ord;
vector<lower = -1, upper = 1>[num_beliefs_ate_pairs] ate_1ord;

matrix<lower = 0, upper = 1>[num_beliefs_obs, num_dist_group_treatments] obs_prob_2ord; // = inv_logit(centered_obs_beta_2ord * beliefs_treatment_map_design_matrix');
vector<lower = -1, upper = 1>[num_dist_group_treatments] prob_2ord;
vector<lower = -1, upper = 1>[num_beliefs_ate_pairs] ate_2ord;
  
