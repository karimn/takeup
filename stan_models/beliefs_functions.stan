vector calculate_beliefs_latent_predictor(matrix design_matrix, matrix beta, matrix dist_beta, vector dist) {
  if (rows(design_matrix) == 1) {
    return (beta * design_matrix[1]') + ((dist_beta * design_matrix[1]') .* dist);
  } else {
    return rows_dot_product(design_matrix, beta) + (rows_dot_product(design_matrix, dist_beta) .* dist);
  }
}