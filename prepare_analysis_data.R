# Setup -------------------------------------------------------------------

library(magrittr)
library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(rgeos)

library(econometr)

source("analysis_util.R")

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

# Load the data ----------------------------------------------------------
load(file.path("data", "takeup_village_pot_dist.RData"))

# This data was prepared in takeup_field_notebook.Rmd
census.data <- read_rds(file.path("data", "takeup_census.rds")) %>% 
  rename(census.consent = consent) # Rename this to reduce chance of error

baseline.data <- read_rds(file.path("data", "takeup_baseline_data.rds"))
takeup.data <- read_rds(file.path("data", "takeup.rds"))
all.endline.data <- read_rds(file.path("data", "all_endline.rds"))
reconsent.data <- read_rds(file.path("data", "reconsent.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds")) # Study clusters metadata
sms.content.data <- read_rds(file.path("data", "takeup_sms_treatment.rds"))
pot.info <- read_rds(file.path("data", "pot_info.rds")) # Point of treatment information
wtp.data <- read_rds(file.path("data", "takeup_wtp.rds")) # Willingness-to-pay experiment data

# Geographic data for clusters in the study
pot_geo_info <- pot.info %>% 
  filter(!is.na(wave)) %>% 
  transmute(cluster.id, 
            pot.lon = lon.post.rct.update,
            pot.lat = lat.post.rct.update) 
  
# HH distance to PoT ------------------------------------------------------

# census.data %>% 
#   group_by(cluster.id) %>% 
#   group_modify(function(households, ...) {
#     hh_sf <- st_as_sf(households, coords = c("lon", "lat"), crs = wgs.84) %>% 
#       st_transform(kenya.proj4)
#     
#     # pot_sf <- st_as_sf(pot_geo_info, coords = c("pot.lon", "pot.lat"), crs = wgs.84)
#     
#     return(households)
#   }, pot_geo_info = pot_geo_info) %>% 
#   ungroup()
#   
#   
#   
#   nest_join(pot_geo_info, ., by = "cluster.id", name = "households") %>% 
#   # st_as_sf(coords = c("pot.lon", "pot.lat"), crs = wgs.84) %>% 
#   mutate(households = map(households, st_as_sf, coords = c("lon", "lat"), crs = wgs.84) %>% 
#            map2(geometry, ~ mutate(.x, 
#                                    dist.to.pot = st_distance(.x %>% st_transform(kenya.proj4),
#                                                              st_sfc(.y, crs = wgs.84) %>% st_transform(kenya.proj4)) %>% 
#                                      as.vector()))) 

# Calculate the distance between households and their assigned point-of-treatment 
census.data %<>% 
  left_join(pot_geo_info, by = "cluster.id") %>% 
  group_by(cluster.id) %>%
  group_modify(~ {
    if (is.na(first(.y$cluster.id))) {
      individ.dist.to.pot <- NA
    } else {
      cluster.pot.location <- convert.to.sp(slice(.x, 1), ~ pot.lon + pot.lat, wgs.84) %>% 
        spTransform(kenya.proj4)
      
      individ.locations <- convert.to.sp(.x, ~ lon + lat, wgs.84) %>% 
        spTransform(kenya.proj4)
      
      individ.dist.to.pot <- c(gDistance(individ.locations, cluster.pot.location, byid = TRUE))
    }
  
    mutate(.x, dist.to.pot = individ.dist.to.pot)
  }) %>% 
  ungroup()

# Village center distance to PoT ------------------------------------------

village.centers %<>% 
  select(cluster.id, lon, lat) %>% 
  left_join(pot_geo_info, by = "cluster.id") %>% {
    cluster.pot.locations <- convert.to.sp(., ~ pot.lon + pot.lat, wgs.84) %>% 
      spTransform(kenya.proj4)
    
    cluster.center.locations <- convert.to.sp(., ~ lon + lat, wgs.84) %>% 
      spTransform(kenya.proj4)
    
    mutate(., dist.to.pot = diag(gDistance(cluster.center.locations, cluster.pot.locations, byid = TRUE)))
  } %>% 
  left_join(select(cluster.strat.data, cluster.id, assigned.treatment, dist.pot.group), by = "cluster.id")

# Core data prep --------------------------------------------------------------

baseline.data %<>% 
  prepare.baseline.data(cluster.strat.data) 

endline.data <- prepare.endline.data(all.endline.data, census.data, cluster.strat.data)
consent.dewormed.reports <- prepare.consent.dewormed.data(all.endline.data, reconsent.data)
analysis.data <- prepare.analysis.data(census.data, takeup.data, endline.data, baseline.data, consent.dewormed.reports, cluster.strat.data, max.name.match.cost = 1)

# Outliers ----------------------------------------------------------------

cell.takeup.data <- prepare.cluster.takeup.data(analysis.data, monitored.only = FALSE, 
                                                add_group_by = c("phone_owner", "mon_status"))

outlier.cells <- cell.takeup.data %>% 
  filter(outlier) 

# WTP prep ----------------------------------------------------------------

wtp.data %<>% 
  mutate(first_choice = factor(first_choice, levels = 1:2, labels = c("bracelet", "calendar")),
         second_choice = if_else(first_choice == "bracelet", cal_plus_ksh, bra_plus_ksh) %>% 
           factor(levels = 1:2, labels = c("switch", "keep"))) %>% 
  select(-county) %>% # I want to use the more reliable one in cluster.strat.data 
  left_join(select(cluster.strat.data, cluster.id, assigned.treatment, dist.pot.group, county), by = "cluster.id")

# Social info report in SMS -----------------------------------------------

social.info.data <- sms.content.data %>% 
  filter(!is.na(social.info)) %>% 
  left_join(select(census.data, cluster.id, village, assigned.treatment, sms.treatment, KEY.individ), by = "KEY.individ") %>% 
  filter(!is.na(cluster.id)) %>% 
  left_join(select(cluster.strat.data, cluster.id, dist.pot.group), by = "cluster.id") %>% 
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

# Knowledge tables prep ---------------------------------------------------

# endline.know.table.data <- read_rds(file.path("data", "know_tables.rds"))

endline.know.table.data <- bind_rows(table.A = read_csv(file.path("data", "raw-data", "Endline Survey-survey-sec_D-tableA.csv")) %>% 
                                       transmute(num.recognized = recogniseA,
                                                 dewormed = dewormedA,
                                                 second.order = order_2ndA,
                                                 second.order.reason = order_2nd_reason,
                                                 relationship = relationshipA,
                                                 relationship.other = relationshipA_other,
                                                 times.seen = times_seenA,
                                                 visited = visitedA, 
                                                 know.other.index = instanceA,
                                                 PARENT_KEY, KEY),
                                     table.B = read_csv(file.path("data", "raw-data", "Endline Survey-survey-sec_D-tableB.csv")) %>% 
                                       transmute(num.recognized = recogniseB,
                                                 dewormed = dewormedB,
                                                 dewormed.know.only = dewormedBB,
                                                 know.other.index = instanceB + 0:9, 
                                                 know.other.index.2 = know.other.index + 1, 
                                                 PARENT_KEY, KEY),
                                     .id = "know.table.type") %>% 
  mutate(recognized = num.recognized > 0)

survey.know.list <- file.path("instruments", "SurveyCTO Forms", "Endline Survey", "Deployed Form Version", 
                              c("knowledge_list.csv", "kak_knowledge_list.csv")) %>% 
  map_df(read_csv, col_types = cols(knowledge_person_key = col_character(), survey.type = col_character()), .id = "wave") %>% 
  mutate(wave = as.integer(wave)) %>% 
  filter(wave == 1 | !is.na(survey.type))

endline.know.table.data %<>% 
  inner_join(select(endline.data, KEY.individ, wave, cluster.id, assigned.treatment, sms.treatment, dist.pot.group, KEY), 
             by = c("PARENT_KEY" = "KEY")) %>% 
  inner_join(select(survey.know.list, person_key, know.other.index, KEY.individ.other),
             by = c("KEY.individ" = "person_key", "know.other.index")) %>% 
  left_join(select(survey.know.list, person_key, know.other.index, KEY.individ.other), # Get the second person (table B)
            by = c("KEY.individ" = "person_key", "know.other.index.2" = "know.other.index"),
            suffix = c(".1", ".2")) %>% 
  # mutate_at(vars(recognized, dewormed, dewormed.know.only), 
  #           list(~ factor(., levels = c(0:1, 98), labels = c("no", "yes", "don't know")))) %>% 
  mutate(across(c(recognized, dewormed, dewormed.know.only), factor, levels = c(0:1, 98), labels = c("no", "yes", "don't know"))) %>% 
  # mutate_at(vars(second.order), 
  #           list(~ factor(., levels = c(1:2, 97:98), labels = c("yes", "no", "prefer not say", "don't know")))) %>% 
  mutate(across(second.order, factor, levels = c(1:2, 97:98), labels = c("yes", "no", "prefer not say", "don't know"))) %>% 
  # mutate_at(vars(relationship), 
  #           list(~ factor(., 
  #                       levels = c(1:5, 99), 
  #                       labels = c("hh member", "extended family", "friend", "neighbor", "church", "other")))) %>%
  mutate(across(relationship, factor, levels = c(1:5, 99), labels = c("hh member", "extended family", "friend", "neighbor", "church", "other"))) %>%
  mutate(relationship = if_else(fct_match(relationship, "other") & fct_match(str_to_lower(relationship.other), c("village member", "village mate", "same village", "village elder", "village mates", "villager")),
                                "village member", as.character(relationship)) %>% factor,
         dewormed = if_else(fct_match(know.table.type, "table.B") & is.na(dewormed), dewormed.know.only, dewormed)) %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, respondent.dewormed.any = dewormed.any),
            by = "KEY.individ") %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, actual.other.dewormed.any.1 = dewormed.any),
            by = c("KEY.individ.other.1" = "KEY.individ")) %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, actual.other.dewormed.any.2 = dewormed.any),
            by = c("KEY.individ.other.2" = "KEY.individ")) %>% 
  mutate(actual.other.dewormed.any.either = actual.other.dewormed.any.1 | (fct_match(know.table.type, "table.B") & actual.other.dewormed.any.2)) 
  # left_join(select(cluster.strat.data, cluster.id, assigned.treatment, dist.pot.group), "cluster.id") # Get it from endline data 

# Save data ---------------------------------------------------------------

save(all.endline.data, endline.data, consent.dewormed.reports, baseline.data, 
     analysis.data, 
     cell.takeup.data, outlier.cells,
     census.data, reconsent.data, takeup.data, sms.content.data, social.info.data, village.centers, wtp.data, endline.know.table.data,
     file = file.path("data", "analysis.RData"))

endline.data %>% 
  set_names(str_replace_all(names(.), "\\.|-", "_")) %>% 
  mutate_if(is.list, list(~ map_chr(., ~ str_c(.x, collapse = " ")))) %>%  
  haven::write_dta(file.path("data", "endline_data.dta"))

analysis.data %>% 
  set_names(str_replace_all(names(.), "\\.|-", "_")) %>% 
  haven::write_dta(file.path("data", "analysis_data.dta"))