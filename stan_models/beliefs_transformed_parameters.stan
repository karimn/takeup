matrix[num_beliefs_obs, num_treatments] centered_obs_beta_1ord; 
matrix[num_beliefs_obs, num_treatments] centered_obs_beta_2ord;  

matrix[num_clusters, num_treatments] centered_cluster_beta_1ord; 
matrix[num_clusters, num_treatments] centered_cluster_beta_2ord;  

matrix[num_beliefs_obs, num_treatments] centered_obs_dist_beta_1ord; 
matrix[num_beliefs_obs, num_treatments] centered_obs_dist_beta_2ord;  

matrix[num_clusters, num_treatments] centered_cluster_dist_beta_1ord; 
matrix[num_clusters, num_treatments] centered_cluster_dist_beta_2ord;  

{
  matrix[num_counties, num_treatments] stratum_beta_1ord;
  matrix[num_clusters, num_treatments] cluster_beta_1ord;
  matrix[num_beliefs_obs, num_treatments] obs_beta_1ord;
  
  matrix[num_counties, num_treatments] stratum_beta_2ord;
  matrix[num_clusters, num_treatments] cluster_beta_2ord;
  matrix[num_beliefs_obs, num_treatments] obs_beta_2ord;
  
  matrix[num_counties, num_treatments] stratum_dist_beta_1ord;
  matrix[num_clusters, num_treatments] cluster_dist_beta_1ord;
  matrix[num_beliefs_obs, num_treatments] obs_dist_beta_1ord;
  
  matrix[num_counties, num_treatments] stratum_dist_beta_2ord;
  matrix[num_clusters, num_treatments] cluster_dist_beta_2ord;
  matrix[num_beliefs_obs, num_treatments] obs_dist_beta_2ord;
  
  vector[num_beliefs_obs] obs_beta_common = beliefs_use_indiv_intercept ?  obs_beta_common_raw * obs_beta_common_sd[1] : rep_vector(0, num_beliefs_obs);
  
  stratum_beta_1ord = beliefs_use_stratum_level ? stratum_beta_1ord_raw .* rep_matrix(stratum_beta_1ord_sd, num_counties) : rep_matrix(0, num_counties, num_treatments);
  cluster_beta_1ord = beliefs_use_cluster_level ? cluster_beta_1ord_raw .* rep_matrix(cluster_beta_1ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_treatments);
  obs_beta_1ord = beliefs_use_obs_level ? obs_beta_1ord_raw .* rep_matrix(obs_beta_1ord_sd, num_beliefs_obs) : rep_matrix(0, num_beliefs_obs, num_treatments);
  
  if (beliefs_use_dist) {
    stratum_dist_beta_1ord = beliefs_use_stratum_level ? stratum_dist_beta_1ord_raw .* rep_matrix(stratum_dist_beta_1ord_sd, num_counties) : rep_matrix(0, num_counties, num_treatments);
    cluster_dist_beta_1ord = beliefs_use_cluster_level ? cluster_dist_beta_1ord_raw .* rep_matrix(cluster_dist_beta_1ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_treatments);
    obs_dist_beta_1ord = beliefs_use_obs_level ? obs_dist_beta_1ord_raw .* rep_matrix(obs_dist_beta_1ord_sd, num_beliefs_obs) : rep_matrix(0, num_beliefs_obs, num_treatments);
  } 
  
  stratum_beta_2ord = beliefs_use_stratum_level ? stratum_beta_2ord_raw .* rep_matrix(stratum_beta_2ord_sd, num_counties) : rep_matrix(0, num_counties, num_treatments);
  cluster_beta_2ord = beliefs_use_cluster_level ? cluster_beta_2ord_raw .* rep_matrix(cluster_beta_2ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_treatments);
  obs_beta_2ord = beliefs_use_obs_level ? obs_beta_2ord_raw .* rep_matrix(obs_beta_2ord_sd, num_beliefs_obs) : rep_matrix(0, num_beliefs_obs, num_treatments);
  
  if (beliefs_use_dist) {
    stratum_dist_beta_2ord = beliefs_use_stratum_level ? stratum_dist_beta_2ord_raw .* rep_matrix(stratum_dist_beta_2ord_sd, num_counties) : rep_matrix(0, num_counties, num_treatments);
    cluster_dist_beta_2ord = beliefs_use_cluster_level ? cluster_dist_beta_2ord_raw .* rep_matrix(cluster_dist_beta_2ord_sd, num_clusters) : rep_matrix(0, num_clusters, num_treatments);
    obs_dist_beta_2ord = beliefs_use_obs_level ? obs_dist_beta_2ord_raw .* rep_matrix(obs_dist_beta_2ord_sd, num_beliefs_obs) : rep_matrix(0, num_beliefs_obs, num_treatments);
  } 
  
  centered_cluster_beta_1ord = rep_matrix(hyper_beta_1ord, num_clusters) + stratum_beta_1ord[cluster_county_id] + cluster_beta_1ord;
  centered_obs_beta_1ord = centered_cluster_beta_1ord[beliefs_cluster_index] + obs_beta_1ord + rep_matrix(obs_beta_common, num_treatments); 
 
  if (beliefs_use_dist) { 
    centered_cluster_dist_beta_1ord = rep_matrix(hyper_dist_beta_1ord, num_clusters) + stratum_dist_beta_1ord[cluster_county_id] + cluster_dist_beta_1ord;
    centered_obs_dist_beta_1ord = centered_cluster_dist_beta_1ord[beliefs_cluster_index] + obs_dist_beta_1ord;
  } else {
    centered_cluster_dist_beta_1ord = rep_matrix(0, num_clusters, num_treatments); 
    centered_obs_dist_beta_1ord = rep_matrix(0, num_beliefs_obs, num_treatments);
  }
  
  centered_cluster_beta_2ord = rep_matrix(hyper_beta_2ord, num_clusters) + stratum_beta_2ord[cluster_county_id] + cluster_beta_2ord;
  centered_obs_beta_2ord = centered_cluster_beta_2ord[beliefs_cluster_index] + obs_beta_2ord + rep_matrix(obs_beta_common, num_treatments); 
  
  if (beliefs_use_dist) { 
    centered_cluster_dist_beta_2ord = rep_matrix(hyper_dist_beta_2ord, num_clusters) + stratum_dist_beta_2ord[cluster_county_id] + cluster_dist_beta_2ord;
    centered_obs_dist_beta_2ord = centered_cluster_dist_beta_2ord[beliefs_cluster_index] + obs_dist_beta_2ord;
  } else {
    centered_cluster_dist_beta_2ord = rep_matrix(0, num_clusters, num_treatments); 
    centered_obs_dist_beta_2ord = rep_matrix(0, num_beliefs_obs, num_treatments);
  }
}