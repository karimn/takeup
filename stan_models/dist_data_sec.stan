int<lower = 0, upper = 1> fit_dist_model_to_data;
int<lower = 0, upper = 1> lognormal_dist_model;
int<lower = 1> num_dist_group_mix;

vector[num_clusters] cluster_standard_dist; // Standardized distance to treatment

// Hyperparameters

real hyper_dist_mean_mean;
real<lower = 0> hyper_dist_mean_sd;
real<lower = 0> hyper_dist_sd_sd;