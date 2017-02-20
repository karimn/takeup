# Setup -------------------------------------------------------------------

library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(lubridate)
library(sp)
library(rgeos)

library(econometr)

source("analysis_util.R")

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

# Load the data ----------------------------------------------------------
load(file.path("data", "takeup_village_pot_dist.RData"))

# This data was prepared in takeup_field_notebook.Rmd
census.data <- read_rds(file.path("data", "takeup_census.rds"))
baseline.data <- read_rds(file.path("data", "takeup_baseline_data.rds"))
takeup.data <- read_rds(file.path("data", "takeup.rds"))
all.endline.data <- read_rds(file.path("data", "all_endline.rds"))
reconsent.data <- read_rds(file.path("data", "reconsent.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))
sms.content.data <- read_rds(file.path("data", "takeup_sms_treatment.rds"))
pot.info <- read_rds(file.path("data", "pot_info.rds"))

# HH distance to PoT ------------------------------------------------------

census.data %<>% 
  left_join(pot.info %>% 
              filter(!is.na(wave)) %>% 
              transmute(cluster.id, 
                        pot.lon = lon.post.rct.update,
                        pot.lat = lat.post.rct.update),
            "cluster.id") %>% 
  group_by(cluster.id) %>% 
  do({
    if (is.na(first(.$cluster.id))) {
      individ.dist.to.pot <- NA
    } else {
      cluster.pot.location <- convert.to.sp(slice(., 1), ~ pot.lon + pot.lat, wgs.84) %>% 
        spTransform(kenya.proj4)
      
      individ.locations <- convert.to.sp(., ~ lon + lat, wgs.84) %>% 
        spTransform(kenya.proj4)
      
      individ.dist.to.pot <- c(gDistance(individ.locations, cluster.pot.location, byid = TRUE))
    }
  
    mutate(., dist.to.pot = individ.dist.to.pot)
  }) %>% 
  ungroup

# Village center distance to PoT ------------------------------------------

village.centers %<>% 
  select(cluster.id, lon, lat) %>% 
  left_join(pot.info %>% 
              filter(!is.na(wave)) %>% 
              transmute(cluster.id, 
                        pot.lon = lon.post.rct.update,
                        pot.lat = lat.post.rct.update),
            "cluster.id") %>% {
    cluster.pot.locations <- convert.to.sp(., ~ pot.lon + pot.lat, wgs.84) %>% 
      spTransform(kenya.proj4)
    
    cluster.center.locations <- convert.to.sp(., ~ lon + lat, wgs.84) %>% 
      spTransform(kenya.proj4)
    
    mutate(., dist.to.pot = diag(gDistance(cluster.center.locations, cluster.pot.locations, byid = TRUE)))
  } %>% 
  left_join(select(cluster.strat.data, cluster.id, assigned.treatment, dist.pot.group), "cluster.id")

# Other prep --------------------------------------------------------------

endline.data <- prepare.endline.data(all.endline.data, census.data, cluster.strat.data)
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
              num.obs = n(),
              qant.25 = interq.info[1],
              qant.75 = interq.info[2]) 
  }) %>%
  group_by(deworming.day, assigned.treatment) %>% 
  mutate(all.dist.cumul = weighted.mean(cumul, num.obs)) %>% 
  ungroup %>% 
  rename(dewormed.day = deworming.day) 

save(all.endline.data, endline.data, consent.dewormed.reports, analysis.data, baseline.data, cluster.takeup.data, 
     census.data, reconsent.data, takeup.data, sms.content.data, social.info.data, village.centers,
     file = file.path("data", "analysis.RData"))