#' ---
#' title: "TakeUp Pilot Summary Analysis"
#' author: "Anne Karing (Berkeley) and Karim Naguib (Evidence Action)"
#' output: pdf_document
#' ---

#+ include=FALSE

library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(sp)
library(rgeos)
library(rgdal)

source("takeup_pilot_util.R")

load("pilot_forms.RData")

special.days <- bind_rows(data.frame(day = 4:5, cluster = 10, note = "High school"),
                          data.frame(day = 6, cluster = 6, note = "End of bracelets"),
                          data.frame(day = c(6, 3, 5), cluster = c(5, 7, 4), note = "Market day"),
                          data.frame(day = 7, cluster = 6, note = "Church"),
                          data.frame(day = 8, cluster = 10, note = "End of airtime"),
                          data.frame(day = 8, cluster = 6, note = "People waiting in households"))

calc.cumul.take.up <- function(.data, grouping = "cluster") { 
  group_by_(.data, .dots = grouping) %>% 
  arrange(day) %>% 
  mutate(cumul.take.up = cumsum(take.up),
         cumul.take.up.prop = cumsum(take.up.prop),
         current.last.day = day == max(day[!is.na(take.up)])) %>% 
  ungroup
}

pilot.vill.summ.2a <- form.2a.data %>% 
  filter(!d2d, cluster.targeted.village, !is.na(PoT)) %>% 
  left_join(village.data[, c("village.code", "min.dist.pot")], by = "village.code") %>% 
  mutate(dist.range = cut(min.dist.pot, breaks = c(seq(0, 2, by = 0.5), 3, 4, 5, 10), include.lowest = T)) %>%  
  group_by(dist.range, cluster) %>% 
  mutate(unique.vill = !duplicated(village.code),
         dist.range.pop = sum(ifelse(unique.vill, village.pop, 0))) %>% 
  group_by(deworming.day, incentive, mda.type, d2d, add = TRUE) %>% 
  summarize(take.up = n(), dist.range.pop = first(dist.range.pop)) %>% 
  ungroup %>% 
  mutate(take.up.prop = take.up / dist.range.pop,
         day = deworming.day) %>% 
  calc.cumul.take.up(c("cluster", "dist.range")) 

pilot.summ.2a <- form.2a.data %>% 
  mutate(high.school = Age <= 15) %>% 
  filter(cluster != 10 | !high.school) %>% 
  group_by(cluster, cluster.pop, deworming.day, incentive, mda.type, d2d) %>% 
  summarize(take.up = n()) %>% 
  ungroup %>% 
  mutate(take.up.prop = take.up / cluster.pop,
         day = deworming.day) %>% 
  calc.cumul.take.up %>% 
  left_join(special.days, by = c("cluster", "day"))

# pilot.summ.2c <- form.2c.data %>% 
#   rename(cluster = Cluster) %>% 
#   left_join(cluster.pop.data, by = "cluster") %>% 
#   group_by(cluster, clust)

pilot.summ <- read_csv("~/Downloads/Pilot Deworming Stats - Simplified (11).csv") %>% 
  select(-1) %>% 
  set_names(paste0("cluster.", 1:10)) %>% 
  mutate(day = 1:11) %>% 
  gather(cluster, take.up, cluster.1:cluster.10) %>%
  separate(cluster, into = c("cl", "cluster"), sep = "\\.") %>% 
  select(-cl) %>% 
  mutate_each(funs(as.numeric(.)), take.up, cluster) %>% 
  filter(day <= 11) %>% 
  left_join(clust.pop.data, by = c("cluster" = "Cluster")) %>% 
  mutate(take.up.prop = take.up/cluster.pop) %>% 
  calc.cumul.take.up %>% 
  left_join(special.days)

base.plot.pilot.aggregate <- function(.data, y.axis, y.axis.label, title, .group, .color) {
  day.breaks <- seq_len(max(.data$day) - 1)
  day.labels <- day.breaks %>% sprintf("%d\n%s", ., ((. - 1) + as.Date("2016-3-7")) %>% wday(label = TRUE))
  
  ggplot(.data, aes_string(x = "day", y = y.axis, group = .group, color = .color)) + 
    geom_line() +
    geom_point(aes(shape = mda.type), size = 3, stroke = 1) +
    scale_x_continuous("Day", breaks = day.breaks, labels = day.labels) +
    scale_y_continuous(y.axis.label) +
    scale_shape_discrete("Point-of-Treatment", solid = FALSE, labels = c("Village Central Location", "Clinic", "Door-to-door", "Market", "Market and School", "School")) +
    ggtitle(title) +
    theme(legend.position = "bottom")
}

plot.pilot.aggregate <- function(.data, y.axis = "cumul.take.up.prop", y.axis.label = "Take-up Proportion (of targeted villages)", title = "Cumulative Proportional Take-up") {
  .data %>% 
    base.plot.pilot.aggregate(y.axis, y.axis.label, title, .group = "factor(cluster)", .color = "factor(incentive)") + 
    geom_text_repel(aes(label = sprintf("Cluster %d", cluster)), nudge_x = .5, size = 3, segment.color = NA, data = . %>% filter(current.last.day)) +
    geom_label_repel(aes(label = note), color = "black", nudge_x = -0.2, size = 2.5, data = . %>% filter(!is.na(note))) +
    # geom_label_repel(aes(label = note), color = "black", nudge_x = -0.2, nudge_y = 0.2, size = 2.5, data = . %>% filter(!is.na(note))) +
    geom_point(color = "black", size = 1, data = . %>% filter(!is.na(note))) +
    scale_color_discrete("Incentive/Signal", labels = c("Airtime", "Bracelet", "Ink", "None")) 
}

plot.pilot.village.aggregate <- function(.data, y.axis = "cumul.take.up.prop", y.axis.label = "Take-up Proportion (of targeted villages)", title = "Cumulative Proportional Take-up") {
  .data %>% 
    base.plot.pilot.aggregate(y.axis, y.axis.label, title, .group = "factor(P)", .color = "min.dist.pot") + 
    scale_color_continuous("Distance") +
    facet_wrap(~ incentive)
}

#' # Take-up By Cluster

#+ echo=FALSE, fig.width=8, fig.height=7, warning=FALSE

pilot.summ.2a %>% plot.pilot.aggregate
pilot.summ.2a %>%
  plot.pilot.aggregate("cumul.take.up", "Absolute Take-up", "Cumulative Absolute Take-up")

#+ echo=FALSE
pilot.summ.2a %>% 
  select(cluster, day, cumul.take.up.prop) %>% 
  spread(day, cumul.take.up.prop) %>% 
  mutate_each(funs(ifelse(is.na(.), "", sprintf("%.2f", .))), 2:11) %>% 
  knitr::kable(col.names = c("Cluster", paste("Day", 1:10))) #, digits = 2)

#+ echo=FALSE, fig.width=8, fig.height=7, warning=FALSE, eval=FALSE

form.2a.data %>% 
  # left_join(clust.pop, by = c("cluster" = "Cluster")) %>% 
  ggplot(aes(Age, group = factor(cluster))) +
  # geom_density() +
  geom_density(aes(color = factor(incentive), linetype = mda.type))

form.2a.data %>%
  # left_join(clust.pop, by = c("cluster" = "Cluster")) %>% 
  ggplot(aes(Age, group = deworming.day)) +
  # geom_density() +
  geom_density(aes(fill = factor(incentive), color = factor(deworming.day), linetype = mda.type), alpha = 0.1) +
  facet_wrap(~ cluster)

#' # Take-up Daily Rates

#+ echo=FALSE, fig.width=8, fig.height=7, warning=FALSE

plot.daily.take.up <- function(.data, y.var="take.up.prop", y.label="Take-up Proportion") { 
    ggplot(.data, aes_string(x = "deworming.day", y = y.var, color = "factor(prop.type)")) + 
      geom_ribbon(aes_string(ymax = y.var, ymin = "0", fill = "factor(prop.type)"), alpha = 0.2) +
      scale_x_continuous("Day", breaks = 1:10) +
      ylab(y.label) +
      scale_color_discrete("") +
      scale_fill_discrete("") +
      facet_grid(d2d ~ incentive, margins = FALSE, labeller = labeller(d2d = . %>% ifelse(., "D2D", "Centralized"))) +
      theme(legend.position = "bottom") + 
      ggtitle("Daily Take-up")
}

form.2a.data %>% 
  group_by(d2d, incentive, deworming.day, targeted.village) %>% 
  summarize(take.up = n()) %>% 
  ungroup %>% 
  left_join(clust.pop.data %>% group_by(d2d, incentive) %>% summarize(incentive.pop=sum(cluster.pop)), by = c("incentive", "d2d")) %>% 
  mutate(take.up.prop = take.up / incentive.pop,
         prop.type = ifelse(targeted.village, "Targeted", "Not Targeted")) %>% 
  bind_rows(group_by(., d2d, incentive, deworming.day) %>% 
              summarize(take.up = sum(take.up), take.up.prop=sum(take.up.prop), prop.type = "Total")) %>% {
    plot.daily.take.up(.) %>% plot
    plot.daily.take.up(., "take.up", "Absolute Take-up") %>% plot
}  

#' # Distance Analysis

#+ echo=FALSE, cache=TRUE

vill.pot.dist <- form.2a.data %>% 
  filter(!d2d, !is.na(village.code), targeted.village) %>% {
    sp.village <- as.data.frame(select(., village.lat, village.lon))
    sp.pot <- as.data.frame(select(., pot.lat, pot.lon))
    
    coordinates(sp.village) <- sp.village[, c("village.lon", "village.lat")] 
    coordinates(sp.pot) <- sp.pot[, c("pot.lon", "pot.lat")] 
    
    proj4string(sp.village) <- CRS(wgs.84)
    proj4string(sp.pot) <- CRS(wgs.84)
    
    sp.village %<>% spTransform(CRS(kenya.proj4))  
    sp.pot %<>% spTransform(CRS(kenya.proj4))  
    
    mutate(., village.pot.dist = diag(gDistance(sp.village, sp.pot, byid = TRUE))/1000)
  } 

vill.pot.dist.mean <- vill.pot.dist %>% 
  group_by(cluster, deworming.day) %>% 
  mutate(mean.dist = mean(village.pot.dist, na.rm = TRUE))

#+ echo=FALSE, fig.width=8, fig.height=7, warning=FALSE, message=FALSE
pilot.vill.summ.2a %>% 
  ggplot(aes(x = day, y = cumul.take.up.prop, group = dist.range, color = dist.range)) +
  geom_point(aes(shape = incentive)) +
  geom_line() +
  labs(x = "Deworming Day", y = "Cumulative Take-up") +
  scale_color_discrete("Village-PoT Distance Range (km)") +
  scale_shape_discrete("Incentive") +
  facet_wrap(~ cluster, labeller = as_labeller(. %>% sprintf("Cluster %s", .)), scales = "free_y") +
  theme(legend.position = "bottom") +
  ggtitle("Proportional Take-up separated by village-PoT distance ranges")


#+ echo=FALSE, fig.height=10, fig.width=8
vill.pot.dist %>% 
  ggplot(aes(x = factor(deworming.day), y = village.pot.dist, color = factor(incentive))) +
  geom_violin(trim = FALSE) +
  geom_point(aes(y = mean.dist), data = vill.pot.dist.mean) +
  labs(x = "Deworming Day", y = "Village-PoT distance (km)") +
  scale_color_discrete("Incentive") +
  coord_cartesian(ylim = c(0, 5)) +
  facet_wrap(~ cluster, ncol = 2, labeller = as_labeller(. %>% sprintf("Cluster %s", .))) + 
  theme(legend.position = "bottom") +
  ggtitle("Violin density plot of distance between villages and PoT")

#+ echo=FALSE, fig.height=8
vill.pot.dist %>% 
  ggplot(aes(x = village.pot.dist)) + 
  geom_density(aes(fill = factor(incentive)), alpha = 0.2) + 
  coord_cartesian(xlim = c(0, 5)) +
  scale_fill_discrete("Incentive") +
  labs(x = "Distance") +
  facet_wrap(~ cluster, scales = "free_y", ncol = 2, labeller = as_labeller(. %>% sprintf("Cluster %s", .))) +
  theme(legend.position = "bottom") + 
  ggtitle("Distribution of Village-PoT distance over intervention period")

#+ echo=FALSE, fig.height=8, eval=FALSE
vill.pot.dist %>% 
  ggplot(aes(x = deworming.day, y = ..density.., color = village.pot.dist, group = village.pot.dist)) + 
  geom_freqpoly(binwidth = 1) +
  facet_wrap(~ cluster, scales = "free_y", ncol = 2) 
  