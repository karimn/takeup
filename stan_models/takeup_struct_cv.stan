vector[cross_validate ? (use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs) : 0] log_lik;
vector[cross_validate ? (use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs) : 0] log_lik_heldout;

if (cross_validate) {
  log_lik = rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs); 
  log_lik_heldout = rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs); 
  
  if (use_binomial) {
    for (cluster_index_index in 1:num_included_clusters) {
      int cluster_index = included_clusters[cluster_index_index];
      
      log_lik[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[cluster_index]);
    }
    
    for (cluster_index_index in 1:num_excluded_clusters) {
      int cluster_index = excluded_clusters[cluster_index_index];
      
      log_lik_heldout[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[cluster_index]);
    }
  } else {
    vector[num_clusters] temp_log_lik = rep_vector(0, num_clusters);
    
    for (obs_index_index in 1:num_included_monitored_obs) {
      int obs_index = included_monitored_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[cluster_index:cluster_index]);
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_excluded_monitored_obs) {
      int obs_index = excluded_monitored_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[cluster_index:cluster_index]);
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik_heldout[obs_index_index] = curr_log_lik;
      }
    }
    
    if (cluster_log_lik) {
      log_lik = temp_log_lik[included_clusters];
      log_lik_heldout = temp_log_lik[excluded_clusters];
    }
  }
}   
  