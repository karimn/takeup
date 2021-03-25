{
  matrix[num_counties, num_dist_group_treatments] stratum_beta_1ord;
  matrix[num_clusters, num_dist_group_treatments] cluster_beta_1ord;
  matrix[num_beliefs_obs, num_dist_group_treatments] obs_beta_1ord;
  
  matrix[num_counties, num_dist_group_treatments] stratum_beta_2ord;
  matrix[num_clusters, num_dist_group_treatments] cluster_beta_2ord;
  matrix[num_beliefs_obs, num_dist_group_treatments] obs_beta_2ord;
  
  vector[num_beliefs_obs] obs_beta_common = obs_beta_common_raw * obs_beta_common_sd;
  
  stratum_beta_1ord = beliefs_use_stratum_level ? stratum_beta_1ord_raw .* rep_matrix(stratum_beta_1ord_sd, num_counties) : rep_matrix(0, num_counties, num_dist_group_treatments);
  cluster_beta_1ord = beliefs_use_cluster_level ? cluster_beta_1ord_raw .* rep_matrix(cluster_beta_1ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_dist_group_treatments);
  obs_beta_1ord = beliefs_use_obs_level ? obs_beta_1ord_raw .* rep_matrix(obs_beta_1ord_sd, num_beliefs_obs) : rep_matrix(0, num_beliefs_obs, num_dist_group_treatments);
  
  stratum_beta_2ord = beliefs_use_stratum_level ? stratum_beta_2ord_raw .* rep_matrix(stratum_beta_2ord_sd, num_counties) : rep_matrix(0, num_counties, num_dist_group_treatments);
  cluster_beta_2ord = beliefs_use_cluster_level ? cluster_beta_2ord_raw .* rep_matrix(cluster_beta_2ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_dist_group_treatments);
  obs_beta_2ord = beliefs_use_obs_level ? obs_beta_2ord_raw .* rep_matrix(obs_beta_2ord_sd, num_beliefs_obs) : rep_matrix(0, num_beliefs_obs, num_dist_group_treatments);
  
  centered_obs_beta_1ord = 
    rep_matrix(hyper_beta_1ord, num_beliefs_obs) + stratum_beta_1ord[beliefs_stratum_index] + cluster_beta_1ord[beliefs_cluster_index] + obs_beta_1ord + rep_matrix(obs_beta_common, num_dist_group_treatments);
  
  centered_obs_beta_2ord = 
    rep_matrix(hyper_beta_2ord, num_beliefs_obs) + stratum_beta_2ord[beliefs_stratum_index] + cluster_beta_2ord[beliefs_cluster_index] + obs_beta_2ord + rep_matrix(obs_beta_common, num_dist_group_treatments);
}