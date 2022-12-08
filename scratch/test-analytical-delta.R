library(tidyverse)
source("dist_structural_util.R")
library(microbenchmark)
library(truncnorm)
library(testthat)
library(rstan)


stan_owen_t_code <-
  '
  functions {
    vector stan_owen_t(vector x, vector y) {
      return owens_t(x, y);
   }
  }
'
expose_stan_functions(stanc(model_code = stan_owen_t_code))


test_that("Delta Equal at limits", {
  ad = analytical_delta(3, 0.1)
  ad_below = analytical_delta_bounded(3, 0.1, bounds = c(-3, 3))
  ad_b = analytical_delta_bounded(3, 0.1, bounds = c(-10000, 10000))

  ad_b_inf = analytical_delta_bounded(3, 0.1, bounds = c(-Inf, Inf))


  expect_equal(
    ad, ad_b
  )

  expect_equal(
    ad,
    ad_b_inf
  )

  expect_lte(ad_below, 3)

})



ub = 3
lb = -3
bandwidth = 3
df = expand.grid(
  w = seq(from = lb - bandwidth, to = ub + bandwidth, length.out = 100 ),
  u_sd = seq(from = 0.1, to = 2,  by = 0.05)
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    delta = calculate_delta(w, sqrt(1 + u_sd^2), u_sd),
    analytical_delta = analytical_delta(w, u_sd), 
    delta_deriv = calculate_delta_deriv(w, sqrt(1 + u_sd^2), u_sd),
    analytical_delta_deriv = analytical_delta_deriv(w, u_sd, delta_w = analytical_delta),
    delta_bounded = calculate_delta_bounded(w, sqrt(1 + u_sd^2), u_sd, bounds = c(lb, ub)),
    analytical_delta_bounded = analytical_delta_bounded(w, u_sd, bounds = c(lb, ub)),
    analytical_delta_deriv_bounded = analytical_delta_deriv_bounded(w, u_sd, bounds = c(lb, ub), delta_w = analytical_delta_bounded),
    analytical_conv_Fw = analytical_conv_Fw(w, u_sd, bounds = c(lb, ub))
  )




comp_df = df %>%
  pivot_longer(
    delta:analytical_conv_Fw
  ) %>%
  mutate(
    type = if_else(str_detect(name, "analytical_"), "analytical", "numerical"), 
    name = str_remove(name, "analytical_")
  ) %>%
  spread(type, value)

make_plot = TRUE
if (make_plot == TRUE) {
  p_comp_plot = comp_df %>%
    filter(!(name %in% c("delta_deriv_bounded", "conv_Fw"))) %>%
    ggplot(aes(
      x = analytical, 
      y = numerical, 
      colour = name
    )) +
    geom_point() +
    facet_wrap(~name, scales = "free", ncol = 1) +
    theme_bw() +
    guides(colour = "none") +
    geom_abline(linetype = "longdash") +
    labs(
      title = "Numerical vs Analytical Delta Calculations", 
      subtitle = str_glue("Delta bounded by {lb}, {ub}")
    )
  p_comp_plot
  ggsave(
    plot = p_comp_plot,
    "temp-plots/numerical-comp-plot.png",
    width = 8,
    height = 6,
    dpi = 500
  )

  comp_df %>%
    gather(variable, value, analytical, numerical ) %>%
    filter(variable == "analytical") %>%
    filter(name == "delta_deriv_bounded") %>%
    ggplot(aes(
      x = w, 
      y = pmin(value, 3), 
      colour = variable, 
      group = u_sd
    )) +
    facet_wrap(~variable) + 
    geom_line()  +
    geom_vline(xintercept = c(lb, ub), linetype = "longdash") +
    guides(colour = "none") +
    theme_bw() +
    labs(
      title = "Derivative of Delta, Bounded"
    )
  ggsave("temp-plots/delta-deriv-bounded.png", width = 8,  height = 6, dpi = 500)

  comp_df %>%
    gather(variable, value, analytical, numerical ) %>%
    filter(name == "delta_bounded") %>%
    ggplot(aes(
      x = w, 
      y = pmin(value, 3), 
      colour = variable, 
      group = u_sd
    )) +
    facet_wrap(~variable) + 
    geom_line()  +
    theme_bw() +
    guides(colour = "none") +
    geom_vline(xintercept = c(lb, ub), linetype = "longdash") +
    labs(title = "Comparing Delta, Bounded")

ggsave("temp-plots/delta-bounded-numerical-comp.png",
  width = 8,
  height = 6,
  dpi = 500)

}


summ_comp_df = comp_df %>%
  filter(name != "delta_deriv_bounded") %>%
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
    "analytical_delta_bounded" = 
        analytical_delta_bounded(seq(from = -3, to = 3, by = 0.1),
          0.2, c(-3, 3)),
    "numerical_delta_bounded" = calculate_delta_bounded(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2, c(-3, 3)),
    "analytical_delta_deriv_bounded" = 
        analytical_delta_deriv_bounded(seq(from = -3, to = 3, by = 0.1),
          0.2, c(-3, 3))
  )

  benchmark_results
}



