for (dist_group_index in 1:num_discrete_dist) { 
  for (dist_group_mix_index in 1:num_dist_group_mix) {
    if (lognormal_dist_model) {
      group_dist_mean[dist_group_index, dist_group_mix_index] ~ normal(hyper_dist_mean_mean, hyper_dist_mean_sd);
    } else {
      group_dist_mean[dist_group_index, dist_group_mix_index] ~ normal(hyper_dist_mean_mean, hyper_dist_mean_sd) T[0, ];
    }
  }
  
  group_dist_sd[dist_group_index] ~ normal(0, hyper_dist_sd_sd);
}

// for (cluster_index in 1:num_clusters) {
//   int dist_group_pos = 1;
//   int curr_assigned_dist_group = cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 2];
//   
//   if (fit_dist_model_to_data) {
//     cluster_standard_dist[cluster_index] ~ lognormal(group_dist_mean[curr_assigned_dist_group, 1], group_dist_sd[curr_assigned_dist_group, 1]);  
//   }
//   
//   for (dist_group_index in 1:num_discrete_dist) {
//     if (dist_group_index != curr_assigned_dist_group) {
//       missing_cluster_standard_dist[cluster_index, dist_group_pos] ~ lognormal(group_dist_mean[dist_group_index, 1], group_dist_sd[dist_group_index, 1]);
//       
//       dist_group_pos += 1;
//     }
//   }
// }

for (cluster_index in 1:num_clusters) {
  // int dist_group_pos = 1;
  int curr_assigned_dist_group = cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 2];
  
  if (fit_dist_model_to_data) {
    vector[num_dist_group_mix] group_dist_mix_lp;
    
    for (group_dist_mix_index in 1:num_dist_group_mix) {
      if (lognormal_dist_model) {
        // real temp_dist = (log(cluster_standard_dist[cluster_index]) - group_dist_mean[curr_assigned_dist_group, group_dist_mix_index]) / group_dist_sd[curr_assigned_dist_group, group_dist_mix_index]; 
        group_dist_mix_lp[group_dist_mix_index] =
          // normal_lpdf(temp_dist, 0, 1) 
          // - log(group_dist_sd[curr_assigned_dist_group, group_dist_mix_index] * cluster_standard_dist[cluster_index]) 
          lognormal_lpdf(cluster_standard_dist[cluster_index] | group_dist_mean[curr_assigned_dist_group, group_dist_mix_index], group_dist_sd[curr_assigned_dist_group, group_dist_mix_index])
          + log(group_dist_mix[curr_assigned_dist_group, group_dist_mix_index]); 
      } else {
        group_dist_mix_lp[group_dist_mix_index] = 
          normal_lpdf(cluster_standard_dist[cluster_index] | group_dist_mean[curr_assigned_dist_group, group_dist_mix_index], group_dist_sd[curr_assigned_dist_group, group_dist_mix_index])  
          + normal_lccdf(0 | group_dist_mean[curr_assigned_dist_group, group_dist_mix_index], group_dist_sd[curr_assigned_dist_group, group_dist_mix_index]) +
          + log(group_dist_mix[curr_assigned_dist_group, group_dist_mix_index]); 
      }
    }
    
    target += log_sum_exp(group_dist_mix_lp);
  }
  
  // for (dist_group_index in 1:num_discrete_dist) {
  //   if (dist_group_index != curr_assigned_dist_group) {
  //     vector[num_dist_group_mix] missing_group_dist_mix_lp;
  //     
  //     for (group_dist_mix_index in 1:num_dist_group_mix) {
  //       if (lognormal_dist_model) {
  //         real temp_dist = (missing_cluster_standard_dist[cluster_index, dist_group_pos] - group_dist_mean[dist_group_index, group_dist_mix_index]) / group_dist_sd[dist_group_index, group_dist_mix_index]; 
  //         
  //         missing_group_dist_mix_lp[group_dist_mix_index] =
  //           normal_lpdf(temp_dist | 0, 1) 
  //           - log(group_dist_sd[dist_group_index, group_dist_mix_index])
  //           - log(missing_cluster_standard_dist[cluster_index, dist_group_pos])
  //           // lognormal_lpdf(missing_cluster_standard_dist[cluster_index, dist_group_pos] | group_dist_mean[dist_group_index, group_dist_mix_index], group_dist_sd[dist_group_index, group_dist_mix_index])
  //             + log(group_dist_mix[dist_group_index, group_dist_mix_index]); 
  //       } else {
  //         missing_group_dist_mix_lp[group_dist_mix_index] =
  //           normal_lpdf(missing_cluster_standard_dist[cluster_index, dist_group_pos] | group_dist_mean[dist_group_index, group_dist_mix_index], group_dist_sd[dist_group_index, group_dist_mix_index]);
  //           
  //         if (missing_cluster_standard_dist[cluster_index, dist_group_pos] < 0) {
  //           missing_group_dist_mix_lp[group_dist_mix_index] += negative_infinity();
  //         } else {
  //           missing_group_dist_mix_lp[group_dist_mix_index] += normal_lccdf(0 | group_dist_mean[dist_group_index, group_dist_mix_index], group_dist_sd[dist_group_index, group_dist_mix_index])
  //             + log(group_dist_mix[dist_group_index, group_dist_mix_index]); 
  //         }
  //       }
  //     }
  //     
  //     target += log_sum_exp(missing_group_dist_mix_lp);
  //     
  //     dist_group_pos += 1;
  //   }
  // }
}
