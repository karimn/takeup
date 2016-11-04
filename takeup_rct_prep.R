library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(stringr)
library(validate)
library(haven)
library(purrr)
library(ggplot2)
library(ggmap)
library(ggalt)
library(ggrepel)
library(gridExtra)
library(scales)
library(broom)
library(parallel)
library(foreach)
library(sp)
library(rgeos)
library(maptools)
library(lme4)

config <- yaml::yaml.load_file("local_config.yaml")
doParallel::registerDoParallel(cores = config$cores) 
# Constants ---------------------------------------------------------------

kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

# Pilot Social knowledge --------------------------------------------------------

soc.know.data.1 <- read_csv("~/Data/TakeUp/pilot/social_deworming_knowledge_1.csv") %>% 
  set_names(str_replace_all(names(.), "\\s+", ".")) 
soc.know.data.2 <- read_csv("~/Data/TakeUp/pilot/social_deworming_knowledge_2.csv") %>% 
  set_names(str_replace_all(names(.), "\\s+", ".")) 

test.response <- function(resp.1, resp.2) {
  is.na(resp.1) | is.na(resp.2) | (resp.1 == resp.2)
}

soc.know.data <- full_join(soc.know.data.1, soc.know.data.2, by = c("Cluster", "Village.Code", "Table.Number", "ID")) %>% 
  set_names(str_replace(names(.), "\\.x", ".1") %>% str_replace("\\.y", ".2")) %>% 
  filter(test.response(Treated.1, Treated.2), 
         test.response(Day.Treated.1, Day.Treated.2), 
         test.response(Wear.Signal.1, Wear.Signal.2)) %>% 
  mutate_each(funs(factor(., 
                          levels = c("Y", "N", "D", "P", "M"), 
                          labels = c("Yes", "No", "Don't Know if Dewormed", "Don't Know Person", "Multiple Persons"))), starts_with("Treated")) %>% 
  mutate_each(funs(factor(., 
                          levels = c("D", "E", "L"), 
                          labels = c("Don't Know", "Early", "Late"))), starts_with("Day.Treated")) %>% 
  mutate_each(funs(factor(., 
                          levels = c("B", "K", "N", "D"), 
                          labels = c("Bracelet", "Ink", "None", "Don't Know"))), starts_with("Wear.Signal")) %>% 
  set_names(str_replace(names(.), "\\.1", "")) %>% {
    local({ 
      load("~/Code/dewormtheworld/pilot_instr/all.deworm.RData")
      
      left_join(., all.deworm.lists[, c("id", "Cluster", "village.code", "dewormed")], by = c("ID" = "id", "Village.Code" = "village.code", "Cluster"))
    })
  } %>% 
  left_join(clust.pop.data[, c("Cluster", "incentive")], by = "Cluster") %>% 
  mutate(resp.id = paste(Cluster, Village.Code, Table.Number, sep = "-"))

soc.know.data %>%
  group_by(incentive) %>% 
  mutate(incentive.n = n()) %>% 
  group_by(Treated, add = TRUE) %>% 
  summarize(prop = n()/first(incentive.n)) %>% 
  ggplot(aes(x = Treated, y = prop, fill = incentive)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  scale_fill_discrete("Incentive")
  # ggplot(aes(x = Treated, fill = incentive)) +
  # geom_bar(position = position_dodge())

yes.no.dk.treated.data <- soc.know.data %>%
  filter(Treated %in% c("Yes", "No", "Don't Know if Dewormed")) %>% 
  mutate(correct.know = (Treated == "Yes" & dewormed) | (Treated == "No" & !dewormed),
         know.status = ifelse(correct.know, 
                              ifelse(Treated == "Yes", "True Positive", "True Negative"),
                              ifelse(Treated == "Don't Know if Dewormed", 
                                     "Don't Know", 
                                     ifelse(Treated == "Yes" & !dewormed, "False Positive", "False Negative")))) %>% 
  mutate(signal = incentive %in% c("bracelet", "ink"),
         incentive2 = ifelse(signal, "signal", incentive))

yes.no.dk.treated.data %>% 
  filter(Treated %in% c("Yes", "No")) %>% 
  group_by(signal, incentive) %>% 
  mutate(incentive.n = n()) %>% 
  group_by(know.status, add = TRUE) %>% 
  summarize(prop = n()/first(incentive.n)) %>% 
  ggplot(aes(x = know.status, y = prop)) +
  geom_bar(aes(, fill = incentive, color = signal), position = position_dodge(), stat = "identity", alpha = 0.1) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1)) +
  xlab("Response") +
  ylab("Proportion of All Yes/No Responses") +
  scale_fill_discrete("Incentive", labels = c("Airtime", "Bracelet", "Ink", "None")) +
  scale_color_discrete("Social Signal", labels = c("No", "Yes")) +
  ggtitle("Social Deworming Knowledge") +
  theme(legend.position = "bottom")

yes.no.treated.data %>%
  block.bootstrap(10, "resp.id") %>%
  do(tidy(lm(correct.know ~ signal, data = .))) %>%
  # mutate(crse.std.error = bootstrap.vcov.clx()
  group_by(term) %>%
  summarize(mean.est = mean(estimate),
            bstrp.se = sd(estimate),
            bstrp.mde = bstrp.se * 2.8,
            bstrp.stat = mean.est/bstrp.se)

yes.no.treated.data %>% 
  nnet::multinom(know.status ~ signal, data = .) %>% 
  summary

soc.know.data %>%
  group_by(incentive) %>% 
  mutate(incentive.n = n()) %>% 
  group_by(Day.Treated, add = TRUE) %>% 
  summarize(prop = n()/first(incentive.n)) %>% 
  ggplot(aes(x = Day.Treated, y = prop, fill = incentive)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  scale_fill_discrete("Incentive")

know.signal.data <- soc.know.data %>%
  filter(incentive %in% c("bracelet", "ink"), #, "none", "airtime"),
         Treated %in% c("Yes", "No", "Don't Know if Dewormed")) %>% 
  mutate(correct.know = (Treated == "Yes" & dewormed) | (Treated == "No" & !dewormed),
         know.status = ifelse(correct.know, 
                              ifelse(Treated == "Yes", "True Positive", "True Negative"),
                              ifelse(Treated != "Don't Know if Dewormed", 
                                     ifelse(Treated == "Yes" & !dewormed, "False Positive", "False Negative"),
                                     as.character(Treated))),
    wear.approp.signal = ifelse(Wear.Signal %in% c("Bracelet", "Ink"), "Signaling", as.character(Wear.Signal)),
         wear.approp.signal = ifelse(wear.approp.signal == "Signaling", ifelse(str_to_lower(Wear.Signal) == incentive, "Correct Signal", "Incorrect Signal"), wear.approp.signal)) %>% 
  group_by(incentive) %>% 
  mutate(incentive.n = n()) %>% 
  group_by(wear.approp.signal, know.status, add = TRUE) %>% 
  summarize(prop = n()/first(incentive.n))

know.signal.data %>% 
  ggplot(aes(x = wear.approp.signal, y = prop, fill = know.status)) +
  geom_bar(position = position_dodge(), stat = "identity", size = 2) +
  scale_fill_brewer("Dewormed", palette="Spectral") +
  xlab("") +
  facet_wrap(~ incentive)

know.signal.data %>%
  ungroup %>% 
  mutate(wear.approp.signal = relevel(factor(wear.approp.signal), ref = "None")) %>% 
  nnet::multinom(wear.approp.signal ~ know.status, data = .) %>% 
  summary

# Population proportions for drugs ----------------------------------------

rct.counties <- c("Busia", "Siaya", "Kakamega")
busia.subcounties <- c("butula", "nambale", "teso south", "teso north") 
siaya.subcounties <- c("gem", "ugenya", "ugunja")

ke.constit.pop.data <- read_csv("~/Data/TakeUp/ke_constit_pop_2009.csv") %>%
  rename(subcounty = Constituency) %>% 
  filter(subcounty %in% str_to_upper(c("butula", "nambale", "gem", "ugenya"))) %>% 
  mutate(Total = as.numeric(Total)) %>% 
  select(subcounty, Total)

ke.teso.pop.data <- read_csv("~/Data/TakeUp/ke_pop_2009.csv") %>% 
  rename(area = `Area in Sq. Km.`) %>% 
  filter(grepl("TESO (SOUTH|NORTH)", District)) %>% 
  group_by(District) %>% 
  summarize_each(funs(sum), Male, Female, Total, Households, area) %>% 
  rename(subcounty = District) %>% 
  select(subcounty, Total)

ke.siaya.ugenya.pop.data <- read_csv("~/Data/TakeUp/ke_pop_2009.csv") %>% 
  filter(Division == "UGUNJA") %>% 
  summarize_each(funs(sum), Male, Female, Total, Households) %>% 
  mutate(subcounty = "UGUNJA") %>% 
  select(subcounty, Total)

rct.subcounty.pop.data <- bind_rows(ke.constit.pop.data, ke.teso.pop.data, ke.siaya.ugenya.pop.data) %>% 
  rename(pop = Total) %>% 
  mutate(pop.prop = pop/sum(pop),
         pills = ceiling(250000 * pop.prop))

# Geographic admin boundaries ---------------------------------------------

ke.lvl2.adm.data <- read_rds("~/Data/TakeUp/KEN_adm2.rds")
ke.lvl3.adm.data <- read_rds("~/Data/TakeUp/KEN_adm3.rds")

counties.adm.data <- ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ] #, "Vihiga"), ]
subcounties.adm.data <- counties.adm.data[!counties.adm.data$NAME_1 %in% c("Busia", "Siaya") | counties.adm.data$NAME_2 %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ]

ke.ward.lvl.adm.data <- ke.lvl3.adm.data 
ke.ward.lvl.adm.data@data %<>% 
  select(NAME_1, NAME_2, NAME_3) %>% 
  rename(county = NAME_1,
         constituency = NAME_2,
         ward = NAME_3)

rct.ward.lvl.adm.data <- ke.ward.lvl.adm.data %>% 
  magrittr::extract(.$county %in% rct.counties, ) %>% 
  magrittr::extract(!.$county %in% c("Busia", "Siaya") | .$constituency %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ) 

rct.ward.lvl.adm.data$ward.id <- rownames(rct.ward.lvl.adm.data@data)
  

kakamega.wards <- ke.ward.lvl.adm.data %>% 
  magrittr::extract(.$county == "Kakamega", ) 

# counties.adm.data[counties.adm.data$NAME_1 == "Busia" & counties.adm.data$NAME_2 == "Matayos", ] %>% plot(add=T, col=alpha("red", 0.2))
# counties.adm.data[counties.adm.data$NAME_1 == "Vihiga", ] %>% plot(add=T, col=alpha("red", 0.2))
# counties.adm.data[counties.adm.data$NAME_1 == "Siaya" & !counties.adm.data$NAME_2 %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ] %>% plot(add=T, col=alpha("green", 0.2))

# Health facilities geolocations ------------------------------------------

ke.health.facilities <- read_csv("~/Data/TakeUp/ke_health_facil.csv") %>% 
  mutate(lat = str_extract(Geolocation, "(?<=\\()-?\\d+(\\.\\d+)"),
         lon = str_extract(Geolocation, "-?\\d+(\\.\\d+)(?=\\))"),
         facility.type.code = `Facility Type`,
         facility.type = factor(`Facility Type NAME`)) %>% 
  mutate_each(funs(as.numeric), lat, lon) %>% 
  filter(facility.type.code %in% c(1, 3, 4, 6, 8, 9),
         !is.na(lat), !is.na(lon)) %>% 
  as.data.frame 

# ke.health.facilities$facility.type %>% table

sp.ke.health.facilities <- ke.health.facilities
coordinates(sp.ke.health.facilities) <- ~ lon + lat
proj4string(sp.ke.health.facilities) <- proj4string(subcounties.adm.data)

sp.ke.health.facilities@data %<>% 
  select(`Facility Name`, District, County, facility.type) %>% 
  rename(facility.name = `Facility Name`,
         county = County) %>% 
  as.data.frame

hc.mapping.data <- haven::read_dta("~/Data/KE Health Facilities Mapping/KenyaMappingData.dta") %>% 
  # select(County, District, ClinicName, FPrimary, Province, GPSlatitude1, GPSlongitude1, ends_with("1_corrected"), FacilityType) %>% 
  select(County, District, ClinicName, GPSlatitude1, GPSlongitude1, FacilityType) %>% 
  mutate(FacilityType = as_factor(FacilityType)) %>% 
  rename(lon = GPSlongitude1,
         lat = GPSlatitude1, 
         facility.name = ClinicName, 
         facility.type = FacilityType,
         county = County) %>% 
  filter(!is.na(lon), !is.na(lat)) %>% 
  as.data.frame %>% 
  `coordinates<-`(c("lon", "lat")) %>% 
  `proj4string<-`(wgs.84)

# hc.mapping.data %>% 
#   spTransform(kenya.proj4) %>% 
#   magrittr::extract(drop(gWithin(., rct.areas, byid = TRUE)), ) %>%
#   # magrittr::extract(.$County == "KAKAMEGA", ) %>% 
#   spTransform(wgs.84) %>% 
#   coordinates %>%
#   as.data.frame %>% {
#     # ggmap(get_map(location = make_bbox(lon, lat, data = .), source = "google")) +
#     ggmap(get_map(location = make_bbox(x, y, data = t(spTransform(rct.areas, wgs.84)@bbox)), source = "google")) +
#       geom_point(aes(lon, lat), alpha = 0.5, size = 1, data = .) +
#       geom_point(aes(lon, lat), alpha = 0.5, size = 1, color = "red", data = sp.ke.health.facilities %>% magrittr::extract(!is.na(.$County) & .$County == "SIAYA", ) %>% coordinates %>% as.data.frame) 
#   }

rct.health.facilities <- sp.ke.health.facilities %>% 
  magrittr::extract(!is.na(.$County) & .$County == "SIAYA", ) %>% 
  rbind(hc.mapping.data) %>% 
  magrittr::extract(drop(gWithin(spTransform(., kenya.proj4), rct.areas, byid = TRUE)), ) 

# rct.health.facilities$cluster.id <- paste0("hf-", seq_len(nrow(rct.health.facilities)))
rct.health.facilities$cluster.id <- as.character(seq_len(nrow(rct.health.facilities)))
          
#   
#   sp.ke.health.facilities[rowSums(gContains(subcounties.adm.data, sp.ke.health.facilities, byid = TRUE)) == 1, ]
# 
# rct.health.facilities$closest.dist <- gDistance(spTransform(rct.health.facilities, CRS(kenya.proj4)), byid = TRUE) %>% 
#   alply(1, function(dist.row) min(dist.row[dist.row > 0])) %>% 
#   unlist %>% 
#   divide_by(1000)

# ggplot(rct.health.facilities@data) +
#   geom_density(aes(closest.dist))

# rct.health.facil.buffers <- gBuffer(spTransform(rct.health.facilities, CRS(kenya.proj4)), byid = TRUE, width = 5000) %>% 
#   spTransform(CRS(proj4string(subcounties.adm.data)))

# Schools -----------------------------------------------------------------

ke.schools.data <- read_csv("~/Data/TakeUp/Kenya_Primary_Schools.csv") %>% 
  mutate(lat = str_extract(Geolocation, "(?<=\\()-?\\d+(\\.\\d+)"),
         lon = str_extract(Geolocation, "-?\\d+(\\.\\d+)(?=\\))")) %>% 
  rename(school.name = `Name of School`) %>% 
  mutate_each(funs(as.numeric), lat, lon) %>%
  filter(!is.na(County), County != "KAKAMEGA") %>% 
  bind_rows(kakamega.sch.geo.data %>% 
              mutate(school.name = ifelse(OQ3_name_correct == "Yes", SI4_schl, OQ4_name_enterd),
                     school.name = gsub("_", " ", school.name),
                     County = "KAKAMEGA",
                     Geolocation = sprintf("(%f, %f)", lat, lon)) %>% 
              select(-selected) %>% 
              rename(District = SI2_district, 
                     Ward = SI3_ward)) %>% 
  as.data.frame

sp.ke.schools.data <- ke.schools.data
coordinates(sp.ke.schools.data) <- ~ lon + lat
proj4string(sp.ke.schools.data) <- proj4string(subcounties.adm.data)

rct.schools.data <- sp.ke.schools.data[rowSums(gContains(subcounties.adm.data, sp.ke.schools.data, byid = TRUE)) == 1, ] %>% 
  magrittr::extract(.$County %in% c("KAKAMEGA", "SIAYA", "BUSIA"), )
  # magrittr::extract(.$`Status of School` == "PUBLIC", )
rct.schools.data@data$cluster.id <- as.character(seq_len(nrow(rct.schools.data)))

write_rds(rct.schools.data, "takeup_rct_schools.rds")

# Locate schools based on ADM boundary data -------------------------------

school.ward.dict <- rct.schools.data %>% 
  spTransform(kenya.proj4) %>% 
  gWithin(spTransform(rct.ward.lvl.adm.data, kenya.proj4), byid = TRUE) %>% 
  adply(1, . %>% magrittr::extract(rct.schools.data$cluster.id, .) %>% data_frame(cluster.id = .), .id = "ward.id")

rct.schools.data@data %<>% 
  left_join(school.ward.dict, by = "cluster.id") %>% 
  left_join(rct.ward.lvl.adm.data@data, by = "ward.id")

# Getting Ward info for Kakamega County
# rct.schools.data@data <- spTransform(kakamega.wards, kenya.proj4) %>% 
#   gWithin(spTransform(rct.schools.data, kenya.proj4), ., byid = TRUE) %>% 
#   alply(1, . %>% magrittr::extract(rct.schools.data$cluster.id, .) %>% data_frame(cluster.id = .), .dims = TRUE) %>% 
#   map_df(~ ., .id = "ward.id") %>% 
#   merge(kakamega.wards, by.x = "ward.id", by.y = "row.names") %>% 
#   mutate(county = str_to_upper(county)) %>% 
#   left_join(rct.schools.data@data, ., by = c("cluster.id", "County" = "county")) %>% 
#   select(-ward.id) %>% 
#   mutate(District = ifelse(County == "KAKAMEGA", subcounty, District))

# Don't uncomment unless you want to overwrite saved schools, especially their cluster IDs!
# write_rds(rct.schools.data, "rct_schools_data2.rds")

# rct.school.buffers <- gBuffer(spTransform(rct.schools.data, CRS(kenya.proj4)), byid = TRUE, width = 5000) %>% 
#   spTransform(CRS(proj4string(subcounties.adm.data)))

# Plot all RCT schools ----------------------------------------------------

subcounties.adm.data %>% 
  spTransform(wgs.84) %>% 
  tidy %>% 
  make_bbox(long, lat, data = .) %>% 
  get_map(location = ., source = "google") %>% 
  ggmap() +
  geom_point(aes(lon, lat), alpha = 0.5, size = 1, data = ke.schools.data) +
  geom_point(aes(lon, lat), alpha = 0.5, size = 1, color = "red", data = tidy(rct.schools.data))

# School centered clusters ------------------------------------------------

# Move to separate file, for profiling using package lineprof
source("takeup_rct_assign_clusters.R")

# Let's take out pilot PoT and villages
pilot.locations <- village.data %>% 
  select(lat, lon) %>%
  mutate(type = "village") %>% 
  bind_rows(pot.geoloc.data %>% distinct(lat, lon) %>% mutate(type = "pot")) %>% 
  as.data.frame

sp.exclude <- pilot.locations 

coordinates(sp.exclude) <- sp.exclude[, c("lon", "lat")] 
proj4string(sp.exclude) <- CRS(wgs.84)

sp.exclude %<>% buffer.clusters(3000) %>% 
  gUnaryUnion

rct.areas <- subcounties.adm.data %>% 
  spTransform(kenya.proj4) %>% 
  gUnaryUnion #%>% 
  #gDifference(sp.exclude)  

# rct.schools.data %<>%
#   spTransform(kenya.proj4) %>% 
#   magrittr::extract(.$County == "KAKAMEGA" | drop(gWithin(., rct.areas, byid = TRUE)), )

targetable.schools <- rct.schools.data %>% 
  spTransform(kenya.proj4) %>% 
  magrittr::extract(drop(!gWithin(., sp.exclude, byid = TRUE)), )

# prof <- lineprof::lineprof(xx <- get.rct.clusters(rct.schools.data, rct.areas, rct.schools.data,
#                    num.clusters = c(control.ink = 75, bracelet.airtime = 75),
#                    ca.outer.radius = c(control.ink = 2500, bracelet.airtime = 4000),
#                    # ca.outer.radius = c(control.ink = 2500, bracelet = 4000, airtime = 4000),
#                    ca.inner.radius = 2500, # c(control.ink = 2500, bracelet = 2500, airtime = 2500),
#                    cluster.group.order = "random",
#                    cluster.size.tester = generate.cluster.pop.area.tester(),
#                    plot.rct.clusters = FALSE))

zz7 <- foreach(simul.id = 1:10, .errorhandling = "remove") %dopar% {
  school.buffer.radius <- 1000
  min.area.frac <- 0.6
  school.area <- pi * (school.buffer.radius^2)
  
  rct.clusters <- get.rct.clusters(rct.schools.data, NULL, rct.schools.data,
  # rct.clusters <- get.rct.clusters(rct.health.facilities, NULL, rct.schools.data, 
                   num.clusters = c(control.ink = 79, bracelet.airtime = 79),
                   ca.outer.radius = c(control.ink = 3000, bracelet.airtime = 4000), 
                   ca.inner.radius = 2500, 
                   cluster.group.order = "random",
                   cluster.size.tester = generate.cluster.num.schools.tester2(targetable.schools, 
                                                                              school.buffer.radius = school.buffer.radius, 
                                                                              min.num.schools = 1,
                                                                              min.area.frac = min.area.frac),
                   plot.rct.clusters = FALSE) %>% 
    magrittr::extract(.$selected, )
  
  return(rct.clusters)
  
}

# names(zz4) <- 1:10
# cdc.data <- ldply(.data = zz4, .id = "simul.id", .fun = . %>% get.cluster.villages.data %>% group_by(pot.cluster.id, cluster.group) %>% summarize(cdc = first(cluster.dist.cat)) %>% count(cluster.group, cdc))
# 
# cdc.data %>% 
#   group_by(simul.id, cluster.group) %>% 
#   mutate(bal.limit = ceiling(sum(n)/2)) %>% 
#   ungroup %>% 
#   filter(cdc != "mixed") %>% 
#   group_by(simul.id) %>% 
#   summarize(balanced = all(n <= bal.limit))

# zz.plots <- zz4[1:4] %>% purrr::map(~ ggplot.clusters(., pilot.locations = pilot.loc.data)) 

# RCT Assignment ----------------------------------------------------------

rct.cluster.selection <- read_rds("rct_cluster_selection_2.0.rds")

rct.targetable.schools <- get.cluster.villages.data(rct.cluster.selection) %>% 
  left_join(distinct(., pot.cluster.id) %>% 
              select(-c(targeted.cluster.id, village.dist.cat, dist)) %>% 
              left_join(count(., cluster.group, cluster.dist.cat) %>% 
                          spread(cluster.dist.cat, n)) %>% 
              group_by(cluster.group) %>% 
              mutate(total.grp.clusters = n()) %>% 
              group_by(cluster.dist.cat, add = TRUE) %>% 
              mutate(alloc.far = (total.grp.clusters %/% 2) - first(far),
                     alloc.close = (total.grp.clusters %/% 2) + (total.grp.clusters %% 2) - first(close)) %>% 
              mutate(assigned.dist.cat = ifelse(cluster.dist.cat != "mixed",  
                                                cluster.dist.cat, 
                                                sample(c(rep("far", first(alloc.far)), rep("close", first(alloc.close)))))) %>% 
              select(pot.cluster.id, assigned.dist.cat)) %>% 
  left_join(filter(., village.dist.cat == assigned.dist.cat) %>% 
              group_by(pot.cluster.id) %>% 
              sample_n(1) %>% 
              select(targeted.cluster.id) %>% 
              mutate(selected.targeted = TRUE)) %>% 
  mutate(selected.targeted = ifelse(!is.na(selected.targeted), selected.targeted, FALSE))

#write_rds(rct.targetable.schools, "rct_targetable_schools_2.0.rds")
rct.targetable.schools <- read_rds("rct_targetable_schools_2.0.rds")

# rct.cluster.selection@data %<>% left_join(targetable.schools %>% 
#                                             spTransform(wgs.84) %>% 
rct.cluster.selection@data %<>% left_join(rct.schools.data %>% 
                                            spTransform(wgs.84) %>% 
                                            as.data.frame %>% 
                                            select(cluster.id, lon, lat)) %>% 
  rename(pot.lon = lon, pot.lat = lat)

rct.schools.data@data %<>% left_join(rct.targetable.schools, by = c("cluster.id" = "targeted.cluster.id"))

# Select backup clusters --------------------------------------------------

backup.per.strata <- 2

rct.backup.clusters <- rct.cluster.selection@data %>%  
  left_join(rct.schools.data@data %>% 
              filter(selected.targeted) %>% 
              select(pot.cluster.id, village.dist.cat),
            by = c("cluster.id" = "pot.cluster.id")) %>% 
  group_by(cluster.group, village.dist.cat) %>% 
  sample_n(backup.per.strata) %>% 
  mutate(intra.strata.order = sample.int(backup.per.strata)) %>% 
  select(cluster.id, cluster.group, village.dist.cat, intra.strata.order)

# write_rds(rct.backup.clusters, "takeup_rct_backup_clusters_2.0.rds")
rct.backup.clusters <- read_rds("takeup_rct_backup_clusters_2.0.rds")

# Generate PoT and target village survey documents ------------------------

# doParallel::registerDoParallel(cores = 6)

rct.cluster.selection %>% {
    foreach(cluster.index = .$cluster.id) %do% {
      current.cluster <- magrittr::extract(., .$cluster.id == cluster.index, )
      rct.cluster.pdf.filename <- sprintf("%s/rct_cluster_survey_2.1/rct_clust_%s.pdf", getwd(), cluster.index)
      rct.cluster.intermed.dir <- sprintf("%s/rct_cluster_survey_2.1/rct_clust_%s_intermed/", getwd(), cluster.index)
      
      cluster.envir <- new.env() 
      cluster.envir$subcounty.polygon <- subcounties.adm.data %>% 
        magrittr::extract(.$NAME_1 == str_to_upper(current.cluster$County, ))
      cluster.envir$sp.cluster <- current.cluster 
      cluster.envir$target.school <- rct.schools.data %>% 
        magrittr::extract(!is.na(.$pot.cluster.id) & .$pot.cluster.id == cluster.index & .$selected.targeted, )
      cluster.envir$intermed.dir <- rct.cluster.intermed.dir 
      
      rmarkdown::render("takeup_cluster_survey.Rmd", output_file = rct.cluster.pdf.filename, intermediates_dir = rct.cluster.intermed.dir, envir = cluster.envir)
    }
  }

# Generate CSV for selected clusters --------------------------------------

rct.cluster.selection@data %>% 
  left_join(rct.schools.data@data, by = c("cluster.id" = "pot.cluster.id")) %>% 
  filter(selected.targeted) %>% 
  set_names(sub("\\.x", ".pot", names(.))) %>% 
  set_names(sub("\\.y", ".anchor", names(.))) %>% 
  select(starts_with("school.name"), starts_with("cluster.id"), matches("^(county|constituency|ward|Geolocation)", ignore.case = FALSE)) %>% 
  mutate(pot.anchor.same = cluster.id == cluster.id.anchor) %>% 
  select(-cluster.id.anchor) %>% 
  mutate_each(funs(gsub("\\(|\\)", "", .)), starts_with("Geolocation")) %>% 
  separate(Geolocation.pot, c("lat.pot", "lon.pot"), ", ", convert = TRUE) %>% 
  separate(Geolocation.anchor, c("lat.anchor", "lon.anchor"), ", ", convert = TRUE) %>% 
  right_join(rct.schools.data %>% 
               magrittr::extract(!is.na(.$selected.targeted) & .$selected.targeted, ) %>% 
               buffer.clusters(1000) %>% 
               spTransform(wgs.84) %>% 
               adply(1, function(school.buff) t(school.buff@bbox) %>% tidy %>% rename(lon = x, lat = y) %>% mutate(cluster.id = school.buff$pot.cluster.id), .id = NULL) %>% 
               gather(variable, value, -c(cluster.id, .rownames)) %>% 
               unite(temp, .rownames, variable) %>% 
               spread(temp, value) %>% 
               select(cluster.id, starts_with("max_"), starts_with("min_"))) %>% 
  write_csv("takeup_rct_clusters_2.1.csv")

# Satellite view of school regions ----------------------------------------

test.region <- school.regions[sample.int(length(school.regions), 1), ] %>% 
  spTransform(CRS(wgs.84))

get_map(location=coordinates(gCentroid(test.region)), zoom=15, maptype="satellite") %>% 
  ggmap + 
  geom_polygon(aes(x = long, y = lat), color = "black", fill = "green", alpha = 0, data = broom::tidy(test.region))

# Pilot village satellite view ----------------------------------------

cluster.village.plots <- foreach(cluster.index = 1:10) %do% {
  pot.geoloc.data %>% 
    filter(Cluster == cluster.index) %>% {
      pot.data <- .
      # vill.data <- village.data %>% filter(Cluster == cluster.index) %>% as.data.frame
      vill.data <- pilot.loc.data %>% filter(cluster == cluster.index) %>% as.data.frame
     
      ret.plot <- get_map(location = if (!empty(pot.data)) pot.data[1, c("lon", "lat")] else vill.data[1, c("lon", "lat")], 
                    zoom=15, 
                    maptype="satellite") %>% 
        ggmap + 
        geom_point(aes(x = lon, y = lat), shape = 1, color = "blue", data = vill.data) 
      
      if (!empty(pot.data)) 
        ret.plot <- ret.plot + geom_point(aes(x = lon, y = lat), shape = 1, color = "red", data = pot.data) 
      
      return(ret.plot)
    }
}  
  
do.call(gridExtra::grid.arrange, cluster.village.plots)
  
# Testing cluster selection -----------------------------------------------

# candidate.cluster.alloc <- foreach(1:10, .errorhandling = "remove") %dopar% {
  yy<-get.rct.clusters(rct.schools.data, subcounties.adm.data, num.clusters = 100, min.cluster.size = 2, cluster.size.fun = generate.cluster.num.schools.fun(rct.schools.data))
# }
  
test.primary.schools <- rct.schools.data %>%
  magrittr::extract(!duplicated(coordinates(.)), ) %>% 
  spTransform(CRS(kenya.proj4)) %>% 
  gBuffer(byid = TRUE, width = 1000, joinStyle = "ROUND")

zz <- foreach (ii=seq_len(nrow(yy))) %do% {
  gDifference(yy[ii, ], yy[-ii, ]) %>% gContains(test.primary.schools, byid = T) %>% sum
  # plot(gDifference(yy[ii, ], yy[-ii, ])) 
  # plot(test.primary.schools, add = T, col=alpha("blue", 0.1))
}

# Health facilities -------------------------------------------------------

# Subcounty plots ---------------------------------------------------------

# plot(counties.adm.data, border = "red")
# plot(rct.schools.data, pch = 20, col = "grey", add = TRUE)
# rct.health.facil.buffers %>% gUnaryUnion %>% gIntersection(subcounties.adm.data) %>% plot(col = alpha("grey", 0.4), border = "transparent")

rct.health.facil.buffers %>% gIntersection(subcounties.adm.data, byid = TRUE) %>% plot(col = alpha("red ", 0.1), border = "transparent")
rct.health.facil.buffers %>% gUnaryUnion %>% gDifference(subcounties.adm.data, ., byid = TRUE) %>% plot(col = alpha("black", 0.2), border = "transparent", add = TRUE)
plot(subcounties.adm.data, lwd = 2, add = TRUE)
plot(rct.health.facilities, pch = 20, col = "red", add = TRUE)
plot(rct.health.facil.buffers, lty = 3, add = TRUE)

rct.school.buffers %>% gIntersection(subcounties.adm.data, byid = TRUE) %>% plot(col = alpha("blue", 0.05), border = "transparent")
plot(subcounties.adm.data, lwd = 2, add = TRUE)
plot(rct.schools.data, pch = 20, col = "blue", add = TRUE)
plot(rct.school.buffers, lty = 3, add = TRUE)

# Cluster survey pilot maps -----------------------------------------------

survey.pilot.data <- read_csv("~/Downloads/Background PoT_June 2016.xlsx - Copy of Sheet1.csv") %>% 
  select(-Degrees) %>% 
  mutate(location.type = ifelse(grepl("PoT", Village), "PoT", "Village"),
         location.type = ifelse(grepl("alt", Village, ignore.case = TRUE), "Alternative PoT", location.type),
         location.type = ifelse(grepl("anchor", Village, ignore.case = TRUE), "Anchor PoT", location.type),
         Cluster = as.character(Cluster)) %>% 
  rename(lon = `GPS  East`,
         lat = `GPS North`,
         loc.name = Village) %>% 
  mutate(lat = if_else(grepl("'", lat), as.numeric(char2dms(sprintf("0d%s\" N", lat))), as.numeric(lat)),
         lon = if_else(grepl("'", lon), as.numeric(char2dms(sprintf("34d%s\" E", lon))), as.numeric(lon))) %>% 
  left_join(rct.cluster.selection@data[, c("cluster.id", "pot.lon", "pot.lat")], by = c("Cluster" = "cluster.id")) %>% 
  left_join(rct.schools.data %>% 
              tidy %>% 
              filter(selected.targeted) %>% 
              select(pot.cluster.id, lon, lat), by = c("Cluster" = "pot.cluster.id"), suffix = c(".survey", ".original.village"))

plot.cluster.survey.locations <- function(cluster.locations) {
  onekm.buffers <- cluster.locations %>% 
    filter(location.type %in% c("PoT", "Anchor PoT")) %>% 
    `coordinates<-`(~ lon.survey + lat.survey) %>% 
    `proj4string<-`(wgs.84) %>% 
    buffer.clusters(.width = 1000) %>% 
    spTransform(wgs.84) %>% 
    tidy
  
  ggmap(get_googlemap(center = c(cluster.locations$pot.lon[1], cluster.locations$pot.lat[1]), maptype = "hybrid", zoom = 14, scale = 2, key = config$google_api_key)) +
    geom_polygon(aes(x = long, y = lat, group = group), linetype = "dotted", alpha = 0, color = "green", size = 1, data = onekm.buffers) +
    geom_point(aes(x = pot.lon, y = pot.lat), shape = 3, color = "red", size = 4, stroke = 1, data = cluster.locations[1, ]) +
    geom_point(aes(x = lon.original.village, y = lat.original.village), shape = 3, color = "green", size = 4, stroke = 1, data = cluster.locations[1, ]) +
    geom_point(aes(x = lon.survey, y = lat.survey), color = "red", alpha = 0.4, size = 1, stroke = 2, data = cluster.locations %>% filter(location.type == "PoT")) +
    geom_point(aes(x = lon.survey, y = lat.survey), color = "green", alpha = 0.4, size = 1, stroke = 2, data = cluster.locations %>% filter(location.type == "Anchor PoT")) +
    geom_point(aes(x = lon.survey, y = lat.survey), color = "blue", alpha = 0.4, size = 1, stroke = 2, data = cluster.locations %>% filter(location.type == "Alternative PoT")) +
    # geom_label_repel(aes(x = lon.survey, y = lat.survey, label = loc.name), force = 3, color = "pink", segment.color = "pink", data = cluster.locations %>% filter(location.type %in% c("PoT", "Anchor PoT"))) +
    geom_point(aes(x = lon.survey, y = lat.survey), color = "white", alpha = 0.4, size = 1, stroke = 2, data = cluster.locations %>% filter(location.type == "Village")) +
    # geom_label_repel(aes(x = lon.survey, y = lat.survey, label = loc.name), force = 3, segment.color = "white", data = cluster.locations %>% filter(location.type == "Village")) +
    labs(x = "", y = "") +
    ggtitle(sprintf("Cluster %s", cluster.locations$Cluster[1]))
}

cluster.survey.plots <- survey.pilot.data %>% 
  group_by(Cluster) %>% 
  do(cluster.location.plot = plot.cluster.survey.locations(.)) 


# Pre-census --------------------------------------------------------------

pre.census.data <- read_csv("~/Code/dewormtheworld/Pre-Census.csv") %>% 
  left_join(read_csv("~/Code/dewormtheworld/Pre-Census-household.csv"), c("KEY" = "PARENT_KEY"), suffix = c(".cluster", ".hh")) %>% 
  rename(lat = `GPS-Latitude`,
         lon = `GPS-Longitude`)

rct.villages <- read_rds("rct_target_villages_2.0.rds") 
second.rct.village <- read_rds("rct_target_villages_2.0-4.rds")

rct.villages %<>% 
  mutate(new.villages = FALSE) %>% 
  bind_rows(mutate(second.rct.village, new.villages = TRUE))

# Determine actual cluster ID for surveyed households
pre.census.data %<>%
  left_join(filter(., !is.na(lon), !is.na(lat)) %>% {
    transmute(., act.cluster.id = gWithinDistance(convert.to.sp(., ~ lon + lat, wgs.84) %>% 
                                                 spTransform(kenya.proj4), 
                                               convert.to.sp(rct.villages, ~ target.lon + target.lat, wgs.84) %>% 
                                                 spTransform(kenya.proj4), 
                                               dist = 1500, byid = TRUE) %>% 
             alply(2, . %>% 
                     magrittr::extract(rct.villages$cluster.id, .) %>% {
                       ifelse(length(.) == 0, NA, .)
                     }) %>% 
             unlist,
             KEY.hh)  
  },
  by = "KEY.hh")

pre.census.data %<>%
  group_by(act.cluster.id) %>%
  mutate(num.cluster.hhs = n()) %>%
  ungroup 

pre.census.data %>% 
  write_rds("pre.census.processed.rds")

gen.label.data <- . %>% 
  mutate(alt.pot.lon = gps2longitude,   
         alt.pot.lat = gps2latitude) %>% 
  select(target.lon, target.lat, matches("^(alt\\.)?pot\\.(lon|lat)$"), cluster.id, target.village.id) %>% 
  gather(key, value, -c(cluster.id, target.village.id)) %>% 
  tidyr::extract(key, c("loc.type", "coord.type"), "^((?:(?:alt\\.)?pot)|target)\\.(lat|lon)") %>% 
  spread(coord.type, value) %>% 
  mutate(target.village.id = ifelse(loc.type %in% c("pot", "alt.pot"), NA, target.village.id))
  
label.data <- gen.label.data(rct.villages) %>% 
  mutate(rct.village = TRUE) %>% 
  bind_rows(gen.label.data(unused.villages) %>% 
              mutate(rct.village = FALSE)) %>% 
  filter(loc.type == "target" | rct.village) %>% 
  distinct(cluster.id, target.village.id, loc.type, .keep_all = TRUE)

plot.pre.census <- function(.data, .zoom = 13, no.labels = FALSE, drop.hhs = FALSE, drop.other.hhs = FALSE, drop.cluster.range = FALSE, .cluster.id = .data$act.cluster.id[1]) {
  label.data %<>% 
    filter(., loc.type %in% c("pot", "alt.pot") | rct.village | cluster.id == .cluster.id)
  
  get.center.coord <- . %>% 
    make_bbox(lon, lat, data = .) %>% 
    matrix(nrow = 2) %>% 
    rowMeans()
  
  # plot.obj <- tryCatch(ggmap(get_map(location = make_bbox(lon, lat, data = .data) %>% matrix(nrow = 2) %>% rowMeans(), zoom = .zoom, maptype = "hybrid", scale = 2, source = "google", api_key = config$google_api_key)) +
  plot.obj <- ggmap(get_googlemap(center = if(!is.null(.data)) get.center.coord(.data) else filter(label.data, cluster.id == .cluster.id, loc.type == "alt.pot") %>% 
                                    select(lon, lat) %>% 
                                    unlist, 
                                  maptype = "hybrid", 
                                  zoom = .zoom, 
                                  scale = 2, 
                                  key = config$google_api_key)) 
  
  if (!is.null(.data)) {
    plot.obj <- plot.obj + 
      geom_point(aes(lon, lat), shape = 15, alpha = 0.25, size = 2, stroke = 1, color = "blue", data = filter(pre.census.data, is.na(act.cluster.id))) 
    
    if (!drop.hhs) {
      plot.obj <- plot.obj + 
        geom_point(aes(lon, lat), shape = 15, alpha = 0.25, size = 2, stroke = 1, color = "green", data = .data) 
    }
    
    if (!drop.other.hhs) {
      plot.obj <- plot.obj +
        geom_point(aes(lon, lat), shape = 15, alpha = 0.2, size = 2, stroke = 1, color = "yellow", data = filter(pre.census.data, act.cluster.id != .data$act.cluster.id[1])) 
    }
   
    plot.obj <- plot.obj +
      ggtitle(sprintf("Cluster %s (%d households)", .cluster.id, nrow(.data)))
  } else {
    plot.obj <- plot.obj +
      ggtitle(sprintf("Cluster %s", .cluster.id))
  }
  
  if (!no.labels) {
    plot.obj <- plot.obj +
      geom_label_repel(aes(lon, lat, color = loc.type, label = cluster.id), segment.color = "white", data = label.data) +
      scale_color_manual("", values = c("red", "darkgrey", "darkgreen"))
  }
  
  plot.obj <- plot.obj + geom_point(aes(lon, lat, fill = loc.type, shape = factor(rct.village, levels = c(TRUE, FALSE))), color = "white", size = 2, stroke = 1, data = label.data)
  
  if (!drop.cluster.range) {
    pot.range.data <- label.data %>% 
      filter(loc.type %in% c("alt.pot"), cluster.id == .cluster.id, !is.na(lat), !is.na(lon))  
    
    if (!empty(pot.range.data)) {
      plot.obj <- plot.obj + geom_polygon(aes(long, lat, group = group), color = "white", linetype = "dotted", alpha = 0,
                                          data = pot.range.data %>% 
                                            convert.to.sp(~ lon + lat, wgs.84) %>% 
                                            buffer.clusters(3000) %>% 
                                            spTransform(wgs.84) %>% 
                                            tidy)
    }
  }
 
  plot.obj <- plot.obj + 
    scale_shape_manual("", values = c(21, 24), drop = FALSE) +
    scale_fill_manual("", values = c("red", "darkgrey", "darkgreen")) +
    labs(x = "", y = "") +
    theme(legend.position = "none") 
  
  
  
  return(plot.obj)
}

pre.census.plots <- pre.census.data %>%
  filter(!is.na(act.cluster.id)) %>%
  # filter(num.cluster.hhs < 100) %>%
  dlply(.(act.cluster.id), function(.sub.data) { tryCatch(plot.pre.census(.sub.data), error = function(err) browser()) }) 

# alply(c("129", "459", "657"), 1, . %>% plot.pre.census(.data = NULL, .cluster.id = .)) %>% 
#   marrangeGrob(ncol = 1, nrow = 1) %>% 
#   ggsave("bad.alt.pot.pdf", ., scale = 2)

pre.census.plots %>%
  compact %>% 
  marrangeGrob(ncol = 1, nrow = 1) %>% 
  ggsave("pre.census.hhs.pdf", ., scale = 2)

# Select secondary villages -----------------------------------------------

second.rct.village <- pre.census.data %>% 
  filter(!is.na(act.cluster.id)) %>%
  count(act.cluster.id) %>% 
  filter(n < 100) %>% 
  anti_join(second.rct.village, by = c("act.cluster.id" = "cluster.id")) %>% 
  inner_join(filter(label.data, loc.type == "target", !rct.village), ., by = c("cluster.id" = "act.cluster.id")) %>% 
  group_by(cluster.id) %>% 
  sample_n(1) %>% 
  ungroup %>% 
  select(cluster.id, target.village.id) %>% 
  left_join(cluster.survey.data) %>% 
  bind_rows(second.rct.village)

# write_rds(second.rct.village, "rct_target_villages_2.0-4.rds")
second.rct.village <- read_rds("rct_target_villages_2.0-4.rds")

second.rct.village %>% write.rct.village.info("takeup_rct_villages-4.csv")

# Baseline backcheck cluster selection ------------------------------------



