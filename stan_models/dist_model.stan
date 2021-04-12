functions {
#include /../multilvlr/util.stan
}

data {
#include base_data_sec.stan
#include dist_data_sec.stan
}

transformed data {
#include base_transformed_data_declare.stan 
}

parameters {
#include dist_parameters_sec.stan
}

model {
#include dist_model_sec.stan 
}

generated quantities {
#include dist_generated_quantities_declare.stan
#include dist_generated_quantities_define.stan
}
