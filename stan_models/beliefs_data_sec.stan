int<lower = 0, upper = 1> beliefs_use_obs_level;
int<lower = 0, upper = 1> beliefs_use_cluster_level;
int<lower = 0, upper = 1> beliefs_use_stratum_level;

int know_table_A_sample_size;

int<lower = 0> num_beliefs_obs;
int<lower = 1, upper = num_obs> beliefs_obs_index[num_beliefs_obs];

int<lower = 0, upper = know_table_A_sample_size> num_recognized[num_beliefs_obs];
int<lower = 0, upper = know_table_A_sample_size> num_knows_1ord[num_beliefs_obs];
int<lower = 0, upper = know_table_A_sample_size> num_knows_2ord[num_beliefs_obs];

matrix[num_treatments * num_discrete_dist, num_treatments * num_discrete_dist] beliefs_treatment_map_design_matrix;

int<lower = 0> num_beliefs_ate_pairs; 
int<lower = 1, upper = num_treatments * num_discrete_dist> beliefs_ate_pairs[num_beliefs_ate_pairs, 2];