data {
  int<lower = 0> num_obs;
  int<lower = 1> num_treatments;
  int<lower = 1> num_clusters;
  int<lower = 1> num_strata;
  
  int<lower = 1, upper = num_clusters> cluster_id[num_obs];
  int<lower = 1, upper = num_strata> stratum_id[num_clusters];
 
  matrix[num_obs, num_treatments] treatment_design_matrix;
 
  // Binary deworming outcome 
  int<lower = 0, upper = 1> dewormed_any[num_obs];
}

parameters {
  vector[num_treatments] beta_cluster[num_clusters];
  vector[num_treatments] beta_stratum[num_strata];
  
  vector[num_treatments] mu;
  
  corr_matrix[num_treatments] Omega_strata;
  vector<lower = 0>[num_treatments] tau_strata;
  
  corr_matrix[num_treatments] Omega_clusters;
  vector<lower = 0>[num_treatments] tau_clusters;
}

model {
  // Some hyperpriors with temp parameters
  
  Omega_strata ~ lkj_corr(2);
  tau_strata ~ cauchy(0, 2.5);
  
  Omega_clusters ~ lkj_corr(2);
  tau_clusters ~ cauchy(0, 2.5);
  
  mu ~ normal(0, 100);
  
  // Model treatment effects hierarchical over clusters and strata
  
  for (stratum_index in 1:num_strata) {
    beta_stratum[stratum_index] ~ multi_normal(mu, quad_form_diag(Omega_strata, tau_strata));
  }
  
  for (cluster_index in 1:num_clusters) {
    beta_cluster[cluster_index] ~ multi_normal(beta_stratum[stratum_id[cluster_index]], quad_form_diag(Omega_clusters, tau_clusters));
  }
  
  // Model the outcome (deworming) in response to experimental treatment
  
  {
    vector[num_obs] treatment_beta_cluster;
    
    for (i in 1:num_obs) {
      treatment_beta_cluster[i] = treatment_design_matrix[i] * beta_cluster[cluster_id[i]];
    }
    
    dewormed_any ~ bernoulli_logit(treatment_beta_cluster);
  }
}
