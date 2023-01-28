array[num_B_treatments, num_mu_treatments] vector[num_optim_distances] optim_w;



for (B_treatment_index in 1:num_B_treatments) {
    for (mu_treatment_index in 1:num_mu_treatments) {
        vector[num_optim_distances]  visibility;
        vector[num_optim_distances] net_private_benefit;


        net_private_benefit = rep_vector(structural_cluster_benefit[1, B_treatment_index], num_optim_distances) - param_dist_cost(
            optim_distances,
            rep_vector(cluster_linear_dist_cost[1, B_treatment_index], num_optim_distances), 
            rep_vector(cluster_quadratic_dist_cost[1, B_treatment_index], num_optim_distances)
        );


        visibility = calculate_mu_rep(
            { mu_treatment_index }, // treatment ids
            optim_distances, // distance
            base_mu_rep, // baseline mu rep
            1, // baseline thing
            beliefs_treatment_map_design_matrix, // 
            rep_matrix(centered_cluster_beta_1ord[1, :], num_optim_distances),
            rep_matrix(centered_cluster_dist_beta_1ord[1, :], num_optim_distances),
            mu_rep_type
        );
        if (USE_MAP_IN_OPTIM) {
            optim_w[B_treatment_index, mu_treatment_index, :] = map_find_fixedpoint_solution(
                net_private_benefit,
                visibility,
                rep_vector(total_error_sd[B_treatment_index], num_optim_distances),
                rep_vector(u_sd[B_treatment_index], num_optim_distances),

                use_u_in_delta,
                alg_sol_rel_tol,
                alg_sol_f_tol,
                alg_sol_max_steps
            );
        } else {
            for (dist_index in 1:num_optim_distances)
            optim_w[B_treatment_index, mu_treatment_index, dist_index] = find_fixedpoint_solution(
                net_private_benefit[dist_index],
                visibility[dist_index],
                total_error_sd[B_treatment_index],
                u_sd[B_treatment_index],

                use_u_in_delta,
                alg_sol_rel_tol,
                alg_sol_f_tol,
                alg_sol_max_steps
            );

        }
    }
}

// vector calculate_mu_rep(array[] int treatment_ids, vector dist,
//                         real base_mu_rep, real mu_beliefs_effect,
//                         matrix design_matrix,
//                         matrix beta, matrix dist_beta, int mu_rep_type) {

//     obs_cluster_mu_rep = calculate_mu_rep(
//       cluster_incentive_treatment_id, cluster_standard_dist, 
//       base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord,
//       mu_rep_type);