for (treatment_index in 1:num_dist_group_treatments) {
  prob_1ord[treatment_index] = mean(obs_prob_1ord[, treatment_index]);
  prob_2ord[treatment_index] = mean(obs_prob_2ord[, treatment_index]);
}

ate_1ord = prob_1ord[beliefs_ate_pairs[, 1]] - prob_1ord[beliefs_ate_pairs[, 2]];
ate_2ord = prob_2ord[beliefs_ate_pairs[, 1]] - prob_2ord[beliefs_ate_pairs[, 2]];
  
