library(tidyverse)
source("dist_structural_util.R")
library(microbenchmark)
library(truncnorm)
library(testthat)

analytical_delta(3, 0.1) 
analytical_delta_bounded(3, 0.1, bounds = c(-1000, 1000))
analytical_delta_bounded(3, 0.1, bounds = c(-Inf, Inf))


ub = 3
lb = -3
bandwidth = 1
df = expand.grid(
  w = seq(from = lb - bandwidth, to = ub + bandwidth, 0.1 ),
  u_sd = seq(from = 0.1, to = 2,  by = 0.05)
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    delta = calculate_delta(w, sqrt(1 + u_sd^2), u_sd),
    analytical_delta = analytical_delta(w, u_sd), 
    delta_deriv = calculate_delta_deriv(w, sqrt(1 + u_sd^2), u_sd),
    analytical_delta_deriv = analytical_delta_deriv(w, u_sd),
    delta_bounded = calculate_delta_bounded(w, sqrt(1 + u_sd^2), u_sd, bounds = c(lb, ub)),
    analytical_delta_bounded = analytical_delta_bounded(w, u_sd, bounds = c(lb, ub))
  )


comp_df = df %>%
  pivot_longer(
    delta:analytical_delta_bounded
  ) %>%
  mutate(
    type = if_else(str_detect(name, "analytical_"), "analytical", "numerical"), 
    name = str_remove(name, "analytical_")
  ) %>%
  spread(type, value)

make_plot = FALSE
if (make_plot == TRUE) {
  p_comp_plot = comp_df %>%
    ggplot(aes(
      x = numerical, 
      y = analytical, 
      colour = name
    )) +
    geom_point() +
    facet_wrap(~name, ncol = 1) +
    theme_bw() +
    guides(colour = "none") +
    geom_abline() +
    labs(
      title = "Numerical vs Analytical Delta Calculations"
    )
  ggsave(
    plot = p_comp_plot,
    "temp-plots/numerical-comp-plot.png",
    width = 8,
    height = 6,
    dpi = 500
  )
}


summ_comp_df = comp_df %>%
  mutate(
    diff = analytical - numerical
  ) %>%
  group_by(name) %>%
  summarise(
    median_diff = median(diff), 
    mean_diff = mean(diff), 
    median_abs_diff = median(abs(diff)), 
    mean_abs_diff = mean(abs(diff))
  )

test_that("Analytical and numerical align", {
  mad = summ_comp_df$median_abs_diff
  mean_abs_diff = summ_comp_df$mean_abs_diff
  map(
    mad,
    expect_lte, 1e-4
  )
  # Deriv a bit different
  map(
    mean_abs_diff,
    expect_lte, 1e-1
  )
})


benchmark = FALSE
if (benchmark == TRUE) {
  benchmark_results = microbenchmark(
    "analytical_delta" = analytical_delta(seq(from = -3, to = 3, by = 0.1), 0.2), 
    "numerical_delta" = calculate_delta(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2), 
    "analytical_delta_deriv" = analytical_delta_deriv(seq(from = -3, to = 3, by = 0.1), 0.2), 
    "numerical_delta_deriv" = calculate_delta(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2),
    "analytical_delta_bounded_vec" = 
        analytical_delta_bounded(seq(from = -3, to = 3, by = 0.1),
          0.2, c(-3, 3)),
    "numerical_delta_bounded" = calculate_delta_bounded(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2, c(-3, 3))
  )

  benchmark_results
}



