row_vector[num_dist_group_treatments] hyper_beta_1ord;
 
matrix[beliefs_use_stratum_level ? num_counties : 0, num_dist_group_treatments] stratum_beta_1ord_raw;
matrix[beliefs_use_cluster_level ? num_clusters : 0, num_dist_group_treatments] cluster_beta_1ord_raw;
matrix[beliefs_use_obs_level ? num_beliefs_obs : 0, num_dist_group_treatments] obs_beta_1ord_raw;
 
row_vector<lower = 0>[beliefs_use_stratum_level ? num_dist_group_treatments : 0] stratum_beta_1ord_sd;
row_vector<lower = 0>[beliefs_use_cluster_level ? num_dist_group_treatments : 0] cluster_beta_1ord_sd;
row_vector<lower = 0>[beliefs_use_obs_level ? num_dist_group_treatments : 0] obs_beta_1ord_sd;
 
row_vector[num_dist_group_treatments] hyper_beta_2ord;
 
matrix[beliefs_use_stratum_level? num_counties : 0, num_dist_group_treatments] stratum_beta_2ord_raw;
matrix[beliefs_use_cluster_level ? num_clusters : 0, num_dist_group_treatments] cluster_beta_2ord_raw;
matrix[beliefs_use_obs_level ? num_beliefs_obs : 0, num_dist_group_treatments] obs_beta_2ord_raw;
 
row_vector<lower = 0>[beliefs_use_stratum_level ? num_dist_group_treatments : 0] stratum_beta_2ord_sd;
row_vector<lower = 0>[beliefs_use_cluster_level ? num_dist_group_treatments: 0] cluster_beta_2ord_sd;
row_vector<lower = 0>[beliefs_use_obs_level ? num_dist_group_treatments : 0] obs_beta_2ord_sd;
 
vector[num_beliefs_obs] obs_beta_common_raw;
real<lower = 0> obs_beta_common_sd;