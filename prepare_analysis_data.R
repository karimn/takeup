library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(lubridate)

source("analysis_util.R")

# This data was prepared in takeup_field_notebook.Rmd
census.data <- read_rds(file.path("data", "takeup_census.rds"))
baseline.data <- read_rds(file.path("data", "takeup_baseline_data.rds"))
takeup.data <- read_rds(file.path("data", "takeup.rds"))
all.endline.data <- read_rds(file.path("data", "all_endline.rds"))
reconsent.data <- read_rds(file.path("data", "reconsent.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))
sms.content.data <- read_rds(file.path("data", "takeup_sms_treatment.rds"))

endline.data <- prepare.endline.data(all.endline.data)
consent.dewormed.reports <- prepare.consent.dewormed.data(all.endline.data, reconsent.data)
analysis.data <- prepare.analysis.data(census.data, takeup.data, endline.data, consent.dewormed.reports, cluster.strat.data)

old.baseline.data <- baseline.data
baseline.data %<>% prepare.baseline.data

cluster.takeup.data <- prepare.cluster.takeup.data(analysis.data)

social.info.data <- sms.content.data %>% 
  filter(!is.na(social.info)) %>% 
  left_join(select(census.data, cluster.id, village, assigned.treatment, sms.treatment, KEY.individ), "KEY.individ") %>% 
  filter(!is.na(cluster.id)) %>% 
  left_join(select(cluster.strat.data, cluster.id, dist.pot.group), "cluster.id") %>% 
  group_by(deworming.day, assigned.treatment, dist.pot.group) %>% 
  do({
    interq.info <- quantile(.$social.info/10, c(0.25, 0.75))
    
    summarize(., 
              cumul = mean(social.info)/10,
              qant.25 = interq.info[1],
              qant.75 = interq.info[2]) 
  }) %>% 
  ungroup %>% 
  rename(dewormed.day = deworming.day) 

save(all.endline.data, endline.data, consent.dewormed.reports, analysis.data, baseline.data, cluster.takeup.data, 
     census.data, reconsent.data, takeup.data, sms.content.data, social.info.data,
     file = file.path("data", "analysis.RData"))