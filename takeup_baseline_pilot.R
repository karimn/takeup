# ---
# output: pdf_notebook
# ---

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)

# Load data ---------------------------------------------------------------

baseline.data <- read_csv("Baseline Survey.csv", col_types = list(date = col_date("%b %d, %Y"))) %>% 
  filter(date >= "2016-08-16")

# Knowledge ---------------------------------------------------------------

baseline.data %>% 
  select(who_worms) %>% 
  mutate(s_who_worms = str_split(who_worms, " ")) %>% 
  unnest %>% 
  mutate(s_who_worms = factor(s_who_worms, levels = c(1:4, 97:99), labels = c("Children", "Adult", "Sick", "Healthy", "Other", "Don't Know", "Prefer not say"))) %>% 
  ggplot +
  geom_bar(aes(s_who_worms))

# Beliefs -----------------------------------------------------------------

baseline.data %>% 
  select(dworm_rate, bet_rate) %>% 
  gather(key, rate) %>% 
  ggplot +
  geom_bar(aes(factor(rate), fill = factor(key)), position = "dodge") 

baseline.data %>% 
  ggplot +
  geom_bar(aes(factor(bet_rate))) +
  labs(x = "Deworm Rate (Bet)")

# Social norms ------------------------------------------------------------

plot.reward.fine.yes.no <- function(.data, question) {
  .data %>%
    select_(.dots = sprintf("%s_%s", c("reward", "fine"), question)) %>% 
    gather(key, val) %>% 
    ggplot() +
    geom_bar(aes(factor(val, levels = 0:1, labels = c("No", "Yes")), fill = key), position = "dodge") +
    theme(legend.position = "bottom") +
    xlab("")
}

plot.reward.fine.amount <- function(.data, question) {
  .data %>%
    select_(.dots = sprintf("%s_%s_amount", c("reward", "fine"), question)) %>% 
    gather(key, val) %>% 
    ggplot() +
    geom_histogram(aes(val, fill = key), binwidth = 100, position = "dodge") +
    theme(legend.position = "bottom") +
    xlab("")
}

baseline.data %>%
  plot.reward.fine.yes.no("well") 

baseline.data %>% 
  plot.reward.fine.amount("well")
  
baseline.data %>%
  plot.reward.fine.yes.no("greet") 

baseline.data %>% 
  plot.reward.fine.amount("greet")

baseline.data %>%
  plot.reward.fine.yes.no("vac") 

baseline.data %>% 
  plot.reward.fine.amount("vac")

baseline.data %>%
  plot.reward.fine.yes.no("worm") 

baseline.data %>% 
  plot.reward.fine.amount("worm")

baseline.data %>%
  plot.reward.fine.yes.no("latrine") 

baseline.data %>% 
  plot.reward.fine.amount("latrine")

baseline.data %>%
  select(starts_with("priority")) %>% 
  gather(key, val) %>% 
  ggplot +
  geom_bar(aes(factor(val, levels = 1:5, labels = c("Well", "Greet", "Vac", "Deworm", "Latrine")), fill = key), position = "dodge") +
  xlab("")
  