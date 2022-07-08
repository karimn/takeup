library(magrittr)
library(tidyverse)
library(EnvStats)

ggplot2::theme_set(theme_minimal())

calculate_delta <- function(v_cutoff, mean1, mean2, sd1, sd2, p.mix = 0.5) {
  F_v <- pnormMix(v_cutoff, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, p.mix = p.mix) 
  
  delta_part <- function(v) v * dnormMix(v, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, p.mix = p.mix) 
  
  e_honor <- integrate(delta_part, lower = v_cutoff, upper = Inf)$value / (1 - F_v)
  e_stigma <- integrate(delta_part, lower = -Inf, upper = v_cutoff)$value / F_v 

  e_honor - e_stigma  
} 

generate_v_cutoff_fixedpoint <- function(b, mu, ...) {
  function(v_cutoff) {
    v_cutoff + b + mu * calculate_delta(v_cutoff, ...)
  }
}

generate_equil_sim <- function(mean1, mean2, sd1, sd2, p.mix, nb_min = -2, nb_max = 2) {
  expand.grid(
    net_benefit = seq(nb_min, nb_max, 0.1),
    mu = seq(0.0, 1.0, 0.05)
  ) %>%
    rowwise() %>% 
    mutate(
      equilibrium_info = list(nleqslv(net_benefit, generate_v_cutoff_fixedpoint(net_benefit, mu, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, p.mix = p.mix))),
      v_cutoff = equilibrium_info$x,
      prob = 1 - pnormMix(v_cutoff, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, p.mix = p.mix),
      term_code = equilibrium_info$termcd,
      delta = calculate_delta(v_cutoff, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, p.mix = p.mix)
    ) %>% 
    ungroup()
}

generate_delta_sim <- function(mean1, mean2, sd1, sd2, p.mix) {
  expand.grid(
    v_cutoff = seq(-1.5, 1.5, 0.1)
  ) %>%
    rowwise() %>% 
    mutate(
      delta = calculate_delta(v_cutoff, mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, p.mix = p.mix)
    ) %>% 
    ungroup()
}

tibble(
  x = seq(-1.75, 1.75, 0.1),
  pdf_bi = dnormMix(x, mean1 = -0.75, sd1 = 0.4375, mean2 = 0.75, sd2 = 0.4375, p.mix = 0.5),
  cdf_bi = pnormMix(x, mean1 = -0.75, sd1 = 0.4375, mean2 = 0.75, sd2 = 0.4375, p.mix = 0.5),
  pdf_uni = dnorm(x),
  cdf_uni = pnorm(x),
) %>%
  pivot_longer(-x, names_pattern = r"{([pc]df)_(uni|bi)}", names_to = c(".value", "v_type")) %>% 
  {
  cowplot::plot_grid(
    ggplot(., aes(x, pdf, color = v_type)) +
      geom_line(show.legend = TRUE)
    
    # ggplot(., aes(x, cdf, color = v_type)) +
    #   geom_line()
  )
} 

equil_sim <- bind_rows(
  bi = generate_equil_sim(-0.75, 0.75, 0.4375, 0.4375, 0.5, nb_min = -3),
  uni = generate_equil_sim(0, 0, 1, 1, 0.0, nb_min = -2.5),
  
  .id = "v_type"
)

equil_sim %>% 
  ggplot(aes(net_benefit, v_cutoff)) +
  geom_point(, data = . %>% filter(term_code != 1)) +
  geom_line(aes(group = mu, color = mu), data = . %>% filter(term_code == 1)) +
  coord_cartesian(ylim = c(-3, 3)) +
  facet_wrap(vars(v_type)) +
  theme(legend.position = "bottom")

equil_sim %>% 
  ggplot(aes(net_benefit, prob)) +
  # geom_point(, data = . %>% filter(term_code != 1)) +
  geom_line(aes(group = mu, color = mu), data = . %>% filter(term_code == 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(vars(v_type)) +
  theme(legend.position = "bottom")

delta_sim <- bind_rows(
  bi = generate_delta_sim(-0.75, 0.75, 0.4375, 0.4375, 0.5),
  uni = generate_delta_sim(0, 0, 1, 1, 0.0),
  
  .id = "v_type"
)

delta_sim %>% 
  ggplot(aes(v_cutoff, delta)) +
  geom_line(aes(color = v_type))
