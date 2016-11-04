library(tidyr)
library(stringr)
library(readr)
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggmap)
library(jsonlite)
library(yaml)
library(sp)
library(rgeos)
library(rgdal)

config <- yaml.load_file("local_config.yaml")

pilot.selected.sch.data <- read_csv("~/Data/TakeUp/pilot_selected_schools.csv") %>% 
  mutate(selected=TRUE)

kakamega.health.facil.data <- read_csv("~/Data/TakeUp/Distribution_of_Health_Facilities_in_Kakamega_County.csv") %>% 
  filter(County == "KAKAMEGA") %>% 
  mutate(lat=str_extract(Geolocation, "(?<=\\()\\d+(\\.\\d+)"),
         lon=str_extract(Geolocation, "\\d+(\\.\\d+)(?=\\))")) %>% 
  mutate_each(funs(as.numeric), lat, lon) %>% 
  rename(facil.type=`Facility Type NAME`)

sp.kakamega.health.facil.data <- as.data.frame(kakamega.health.facil.data)
coordinates(sp.kakamega.health.facil.data) <- ~ lon + lat

kakamega.sch.geo.data <- read_csv("~/Data/TakeUp/Kakamega schools geolocation.csv") %>% 
  rename(lat=OQ2_sch_gps_latitude,
         lon=OQ2_sch_gps_longitude,
         school_id=SI5_schlid) %>% 
  mutate(school_id=gsub("_", "-", school_id)) %>% 
  left_join(pilot.selected.sch.data %>% select(school_id, selected)) %>% 
  mutate(selected=ifelse(!is.na(selected), selected, FALSE)) 

sp.kakamega.sch.geo.data <- as.data.frame(kakamega.sch.geo.data)
coordinates(sp.kakamega.sch.geo.data) <- ~ lon + lat

kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

proj4string(sp.kakamega.sch.geo.data) <- CRS(wgs.84)
sp.kakamega.sch.geo.data %<>% spTransform(CRS(kenya.proj4))

proj4string(sp.kakamega.health.facil.data) <- CRS(wgs.84)
sp.kakamega.health.facil.data %<>% spTransform(CRS(kenya.proj4))

min.max.sch.dist.km <- gDistance(sp.kakamega.sch.geo.data, byid=TRUE) %>% 
  aaply(.margins=1, function(sch.dist) { c(sort(sch.dist, partial=2)[2], max(sch.dist)) }) %>%
  divide_by(1000)

min.max.dist.km <- gDistance(sp.kakamega.sch.geo.data, sp.kakamega.health.facil.data, byid=TRUE) %>% 
  aaply(.margins=1, function(sch.dist) c(min(sch.dist), max(sch.dist))) %>%
  divide_by(1000)

min.max.dist.km.interquart <- quantile(min.max.dist.km[, 1], c(0.25, 0.5, 0.75))
min.max.sch.dist.km.interquart <- quantile(min.max.sch.dist.km[, 1], c(0.25, 0.5, 0.75))

qplot(min.max.dist.km[, 1], geom="density") + 
  geom_vline(xintercept=min.max.dist.km.interquart, linetype="dashed") + 
  labs(x="Distance to nearest health facility (km)", y="Density")

qplot(min.max.sch.dist.km[, 1], geom="density") + 
  geom_vline(xintercept=min.max.dist.km.interquart, linetype="dashed") + 
  labs(x="Distance between schools (km)", y="Density") + 
  coord_cartesian(xlim=c(0, 2))

# get_map("kakamega county kenya", zoom=10, maptype="roadmap") %>% 
get_map(location=c(34.75, 0.5), zoom=10, maptype="roadmap", scale=2) %>%
  ggmap + 
  geom_point(aes(lon, lat), color="red", size=2.5, data=kakamega.sch.geo.data %>% filter(selected)) +
  geom_point(aes(lon, lat), color="black", size=2.5, data=kakamega.sch.geo.data %>% filter(!selected)) +
  geom_point(aes(lon, lat), color="darkgrey", size=1, data=kakamega.sch.geo.data) +
  geom_point(aes(lon, lat), shape=17, color="blue", size=2.75, data=kakamega.health.facil.data) +
  geom_point(aes(lon, lat), shape=17, color="darkgrey", size=1, data=kakamega.health.facil.data) +
  labs(x="", y="")

# get_map(location=c(34.49744, 0.2187383), zoom=19, maptype="satellite", scale=2) %>% 
#   ggmap + 
#   geom_point(x=34.49744, y=0.2187383, shape=3, color="red", size=2.5) 
  
  

place.search <- . %>% 
  gsub(" ", "+", ., fixed=true) %>% 
  sprintf("https://maps.googleapis.com/maps/api/place/textsearch/json?query=%s&key=%s", ., config$google_api_key) %>% 
  fromjson(flatten=true)

clean.sch.name <- . %>% 
  gsub(sprintf("\\s+(%s)\\.?(?=\\s+|$)", paste(c("PRY", "PRI(M(ARY)?)?", "SCH(OOL)?"), collapse="|")), "", ., perl=TRUE)

ke.moe07.data <- read_csv("~/Data/TakeUp/Kenya_Primary_Schools.csv") %>% 
  set_names(gsub("\\s+", ".", tolower(names(.)))) %>% 
  mutate(name.of.school=clean.sch.name(name.of.school))

kakamega.moe07.data <- ke.moe07.data %>% 
  filter(county == "KAKAMEGA") 
 
ke.y3.treat.data <- read_csv("~/Data/TakeUp/Enrolled_Treated_Y3.csv") %>% 
  set_names(gsub("\\s+", ".", tolower(names(.)))) %>% 
  mutate(a_school_name=clean.sch.name(a_school_name)) %>% 
  mutate(geocode.name=sprintf("%s school %s Kenya", a_school_name, county)) %>% 
  mutate_each(funs(toupper), county, district_name, division_name) 
    
kakamega.y3.treat.data <- ke.y3.treat.data %>% 
  filter(county == "KAKAMEGA")

kakamega.health.facil.data <- read_csv("~/Data/TakeUp/Distribution_of_Health_Facilities_in_Kakamega_County.csv")

# kakamega.y3.treat.data %>% anti_join(kakamega.moe07.data, by=c("a_school_name"="name.of.school")) #, "division_name"="division"))

# kakamega.y3.treat.data %<>% 
#   bind_cols(geocode(.$geocode.name, output="more"))
# 
# xx <- kakamega.y3.treat.data %>% 
#   filter(is.na(query)) %>% 
#   rowwise %>% 
#   do(place.info=place.search(.$geocode.name))

pilot.schools <- kakamega.y3.treat.data %>% 
  filter(!is.na(total.registered)) %>% 
  filter(total.registered >= quantile(total.registered, probs=0.25),
         total.registered <= quantile(total.registered, probs=0.75)) %>% # Consider those in interquartile range
  inner_join(group_by(., district_name) %>% 
               summarize(wt=var(total.registered)) %>% 
               sample_n(2, weight=wt)) %>% # Pick two random subcounties weighted by variance of school size 
  # mutate(adult.child.ratio=max(adults.dewormed.sth, adults.dewormed.schisto, na.rm=TRUE)/total.registered) %>% 
  # inner_join(group_by(., division_name) %>% 
  #              summarize(wt=var(total.registered)) %>% 
  #              sample_n(5, weight=wt)) %>% 
  mutate(adult.child.ratio=adults.dewormed.sth/total.registered) %>% 
  # arrange(adult.child.ratio) %>% 
  # slice(c(seq_len(40), seq(from=n() - 40 + 1, to=n()))) %>% 
  # mutate(takeup.type=ifelse(row_number() <= 40, "low", "high")) %>% 
  mutate(takeup.type=ifelse(adult.child.ratio <= median(adult.child.ratio), "low", "high")) %>% 
  group_by(takeup.type) %>% 
  # group_by(district_name, takeup.type) %>% 
  sample_n(20)
  
pilot.loc.data <- read_csv("~/Data/TakeUp/pilot_locations_gps.csv") %>% 
  set_names(gsub("\\s+", ".", names(.)) %>% tolower) %>% 
  mutate(admin.level=tolower(admin.level) %>% sub("\\d+$", "", .) %>% str_trim) %>% 
  select(-c(phone, degrees.east, elevation, accuracy, two.treatment, comments)) %>% 
  filter(!is.na(gps.north)) %>% 
  separate(gps.east, c("gps.east.min", "gps.east.sec"), "'") %>% 
  separate(gps.north, c("gps.north.min", "gps.north.sec"), "'") %>% 
  mutate_each(funs(sub("['\"]+$", "", .)), ends_with("sec")) %>%
  filter(!is.na(gps.east.min) & !is.na(gps.north.min)) %>% 
  mutate(lat=sprintf("0d%s'%s\"N", gps.north.min, gps.north.sec)) %>% 
  mutate(lon=sprintf("34d%s'%s\"E", gps.east.min, gps.east.sec)) %>% 
  mutate_each(funs(as.numeric(char2dms(.))), lat, lon) %>% 
  mutate_each(funs(factor), admin.level, cluster) %>% 
  filter(!admin.level %in% c("sub chief", "contact person", "chief", "kulumbeni")) %>% 
  group_by(cluster) %>% 
  mutate(loc.id=paste(cluster, 1:n(), sep="_")) %>% 
  ungroup

sp.pilot.loc.data <- as.data.frame(pilot.loc.data)
coordinates(sp.pilot.loc.data) <- sp.pilot.loc.data[, c("lon", "lat")] #~ lon + lat

proj4string(sp.pilot.loc.data) <- CRS(wgs.84)

map.cluster <- function(.sp.data, cluster.id, .zoom=15) {
 get_map(location=bbox(.sp.data[.sp.data$cluster %in% cluster.id, ]) %>% as.vector, zoom=.zoom, maptype="hybrid", scale=2) %>%
    ggmap + 
    geom_point(aes(lon, lat, color=admin.level, shape=cluster), alpha=0.5, size=3, stroke=2, data=.sp.data@data %>% filter(cluster %in% cluster.id)) 
  }

find.closest.treat.loc <- function(cluster.data) {
  cluster.data %<>% as.data.frame
  coordinates(cluster.data) <- cluster.data[, c("lon", "lat")]
  proj4string(cluster.data) <- CRS(wgs.84)
  cluster.data %<>% spTransform(CRS(kenya.proj4))  
  
  villages <- cluster.data[cluster.data$admin.level == "village", ]
  treat.loc <- cluster.data[cluster.data$admin.level != "village", ]
    
  dist <- gDistance(villages, treat.loc, byid=TRUE) %>%
    divide_by(1000) %>% 
    adply(.margins=2, function(vill.dist) data.frame(treat.loc.id=which.min(vill.dist), min.dist.km=min(vill.dist))) %>% 
    rename(village.id=X1) %>% 
    mutate(village.id=villages$loc.id[village.id],
           treat.loc.id=treat.loc$loc.id[treat.loc.id])
  
  return(cluster.data@data %>% left_join(dist, by=c("loc.id"="village.id")))
}

vill.treat.min.dist.data <- pilot.loc.data %>% 
  group_by(cluster) %>% 
  do(find.closest.treat.loc(.)) %>% 
  ungroup %>% 
  select(-starts_with("gps")) %>% 
  arrange(cluster, admin.level, min.dist.km)

write_csv(vill.treat.min.dist.data, "pilot_villages.csv", na="")

