int<lower = 1> num_obs; // Actual observations
int<lower = 1> num_clusters;
int<lower = 1> num_counties;
int<lower = 1> num_treatments;
int<lower = 1> num_discrete_dist;

int<lower = 1, upper = num_clusters> obs_cluster_id[num_obs];
int<lower = 1, upper = num_counties> obs_county_id[num_obs];
int<lower = 1, upper = num_counties> cluster_county_id[num_clusters];

int<lower = 1> cluster_treatment_map[num_treatments * num_discrete_dist, 2];
int<lower = 1, upper = num_treatments * num_discrete_dist> cluster_assigned_dist_group_treatment[num_clusters]; 

