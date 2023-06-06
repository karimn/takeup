array[num_B_treatments, num_mu_treatments] vector[num_optim_distances] optim_w;
array[num_B_treatments, num_mu_treatments] vector[num_optim_distances] visibility;
array[num_B_treatments, num_mu_treatments] vector[num_optim_distances] net_private_benefit;


for (B_treatment_index in 1:num_B_treatments) {
    for (mu_treatment_index in 1:num_mu_treatments) {


        net_private_benefit[B_treatment_index, mu_treatment_index] = rep_vector(structural_cluster_benefit[1, B_treatment_index], num_optim_distances) - param_dist_cost(
            optim_distances,
            rep_vector(cluster_linear_dist_cost[1, B_treatment_index], num_optim_distances), 
            rep_vector(cluster_quadratic_dist_cost[1, B_treatment_index], num_optim_distances)
        );



        visibility[B_treatment_index, mu_treatment_index] = calculate_mu_rep(
            rep_array(mu_treatment_index, num_optim_distances), // treatment ids
            optim_distances, // distance
            base_mu_rep, // baseline mu rep
            1, // baseline thing
            beliefs_treatment_map_design_matrix, // 
            rep_matrix(centered_cluster_beta_1ord[1, :], num_optim_distances),
            rep_matrix(centered_cluster_dist_beta_1ord[1, :], num_optim_distances),
            mu_rep_type
        );


    }
}
