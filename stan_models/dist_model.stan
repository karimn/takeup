data {
#include base_data_sec.stan
#include dist_data_sec.stan
}

parameters {
#include dist_parameters_sec.stan
}

model {
#include dist_model_sec.stan 
}
