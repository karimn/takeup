library(tidyverse)
source("dist_structural_util.R")
library(microbenchmark)
library(truncnorm)

  
ub = 3
lb = -3
bandwidth = 0
df = expand.grid(
  w = seq(from = lb - bandwidth, to = ub + bandwidth, 0.1 ),
  u_sd = seq(from = 0, to = 2,  by = 0.05)
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    # delta = calculate_delta(w, sqrt(1 + u_sd^2), u_sd),
    # ed_delta = analytical_delta(w, u_sd), 
    # delta_deriv = calculate_delta_deriv(w, sqrt(1 + u_sd^2), u_sd),
    # ed_delta_deriv = analytical_delta_deriv(w, u_sd),
    # Fw = calculate_Fw_bounded(w, sqrt(1 + u_sd^2), u_sd, bounds = c(lb, ub)),
    # M_plus = calculate_M_plus(w, sqrt(1 + u_sd^2), u_sd, bounds = c(lb, ub)),
    # M_minus = calculate_M_minus(w, sqrt(1 + u_sd^2), u_sd, bounds = c(lb, ub)),
     delta_bounded = calculate_delta_bounded(w, sqrt(1 + u_sd^2), u_sd, bounds = c(lb, ub)),
     ed_delta_bounded = analytical_delta_bounded(w, u_sd, bounds = c(lb, ub))
  )

df %>%
  filter(u_sd > 0) %>%
  ggplot(aes(
    x = ed_delta_bounded,
    y = delta_bounded
  )) +
  geom_point() +
  geom_abline()

  stop()

df %>%
    filter(u_sd > 0) %>%
    ggplot(aes(
        x = w, 
        y = delta_bounded, 
        group = u_sd
    )) +
    geom_line() +
    geom_line(aes(
        y = ed_delta_bounded
    ), 
    colour = "hotpink") +
    geom_hline( 
        yintercept = 3
    )

df %>%
  filter(u_sd > 0) %>%
  ggplot(aes(
    x = ed_delta_deriv,
    y = delta_deriv
  )) +
  geom_point() +
  geom_abline()

df %>%
  filter(u_sd > 0) %>%
  ggplot(aes(
    x = ed_delta,
    y = delta
  )) +
  geom_point() +
  geom_abline()

df %>%
  filter(u_sd > 0) %>%
  gather(variable, value, delta, ed_delta) %>%
  ggplot(aes(
    x = w, 
    y = value, 
    group = interaction(u_sd, variable), 
    colour = variable
  )) +
  geom_line()  +
  facet_wrap(~variable)

df %>%
  filter(u_sd > 0) %>%
  filter(abs(w) < 5) %>%
  gather(variable, value, delta, ed_delta) %>%
  ggplot(aes(
    x = w, 
    y = value, 
    group = interaction(u_sd, variable), 
    colour = variable
  )) +
  geom_line()  


df %>%
  # filter(abs(w) < 3 ) %>%
  group_by(u_sd) %>%
  filter(u_sd > 0) %>%
  # filter(!any(M_minus < -20)) %>%
  ggplot(aes(
    x = w, 
    y = delta_bounded, 
    group = u_sd
  )) +
  geom_line() +
  NULL  +
  labs(
    title = "Bounded Rep Returns"
  ) +
  geom_hline(yintercept = 3, linetype = "longdash") +
  theme_bw()

ggsave(
  "temp-plots/bounded-rep-returns.png",
  width = 8,
  height = 6,
  dpi = 500
)


df %>%
  # filter(abs(w) < 3 ) %>%
  group_by(u_sd) %>%
  filter(u_sd > 0) %>%
  # filter(!any(M_minus < -20)) %>%
  ggplot(aes(
    x = w, 
    y = M_plus - M_minus, 
    group = u_sd
  )) +
  geom_line() +
  NULL  +
  ylim(0, 10) +
  labs(title = "Attempt 2") +
  geom_hline(yintercept = 3, linetype = "longdash")
  # ylim(0, 3.5)



ed_df = tibble(
  w = seq(from = 0, to = 10, by = 0.1),

val = ed_func(
  seq(from = 0, to = 10, by = 0.1),
  -3,
  3,
  0.2
)
)

ed_df %>%
  ggplot(aes(
    x = w, 
    y = val
  )) +
  geom_point()

calculate_M_plus(seq(from = 0, to = 10, by = 0.1), sqrt(1 + 0.2), 0.2, c(-3, 3))
calculate_M_minus(seq(from = 0, to = -10, by = -0.1), sqrt(1 + 0.2), 0.2, c(-3, 3))

calculate_delta_bounded(seq(from = 0, to = -10, by = -0.1), sqrt(1 + 0.2), 0.2, c(-3, 3))

microbenchmark(
  "analytical_delta" = analytical_delta(seq(from = -3, to = 3, by = 0.1), 0.2), 
  "numerical_delta" = calculate_delta(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2), 
  "analytical_delta_deriv" = analytical_delta_deriv(seq(from = -3, to = 3, by = 0.1), 0.2), 
  "numerical_delta_deriv" = calculate_delta(seq(from = -3, to = 3, by = 0.1), sqrt(1 + 0.2^2), 0.2)
)