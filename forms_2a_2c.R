library(readr)
library(plyr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(lubridate)
library(magrittr)
library(sp)
library(rgeos)
library(rgdal)

source("takeup_pilot_util.R")

pot.geoloc.data <- read_csv("~/Data/TakeUp/pilot/pot_geoloc.csv") %>% 
  mutate_each(funs(sub("\\.", "d", 
                       sub("^([NE])(.+\")", "\\2 \\1", 
                           sub("\'{2}", "\"", .), .))), lat, lon) %>% 
  mutate_each(funs(as.numeric(char2dms(.))), lat, lon)

form.status.data <- read_csv("~/Data/TakeUp/pilot/form_status.csv") %>% 
  set_names(c("filename", "PoT", "status", "digitized")) %>% 
  filter(!grepl("D2D$", PoT)) #%>% 
  # left_join(pot.geoloc.data, by="PoT")

form.2a.fn.dict <- read_csv("~/Data/TakeUp/pilot/clust_10_fn_dict.csv")

usable.clusters <- 1:10 

village.codes <- read_csv("~/Data/TakeUp/pilot/village_codes.csv") %>%
  mutate_each(funs(str_to_upper), Village) %>% 
  rename(village.code = Code)

village.data <- read_csv("~/Downloads/Census Data.2016.3.4 - Copy of Sheet1.csv") %>%
  mutate_each(funs(str_to_upper), Village) %>% 
  partial.match.village.code("Village", village.codes) %>% 
  left_join(village.codes[, c("village.code", "pop")], by = "village.code") %>% 
  rename(village.pop = pop) %>% 
  left_join(pilot.loc.data %>% 
              left_join(village.codes[, c("Village", "village.code")], by = c("location" = "Village")) %>% 
              filter(admin.level == "village") %>% 
              select(village.code, lat, lon), 
            by = "village.code") %>% 
  distinct(village.code) %>% {
    sp.village <- as.data.frame(select(., village.code, lat, lon))
    coordinates(sp.village) <- sp.village[, c("lon", "lat")] 
    
    sp.pot <- pot.geoloc.data %>% distinct(lat, lon) %>% as.data.frame
    coordinates(sp.pot) <- sp.pot[, c("lon", "lat")] 
    
    proj4string(sp.village) <- CRS(wgs.84)
    proj4string(sp.pot) <- CRS(wgs.84)
    
    sp.village %<>% spTransform(CRS(kenya.proj4))  
    sp.pot %<>% spTransform(CRS(kenya.proj4))  
    
    mutate(., min.dist.pot = unlist(alply(gDistance(sp.village, sp.pot, byid = TRUE), 2, min)) / 1000)
  }

clust.pop.data <- village.data %>% 
  left_join(village.codes[, c("village.code", "pop")]) %>% 
  group_by(Cluster) %>% 
  # summarize(pop = sum(Population)) %>% 
  summarize(cluster.pop = sum(pop)) %>% 
  arrange(Cluster) %>% 
  mutate(max.days = c(8, 10, 1, 8, 8, 10, 3, 1, 3, 8),
         incentive = c("none", "ink", "airtime", "bracelet", "ink", "bracelet", "bracelet", "none", "ink", "airtime"),
         mda.type = c("central.location", "d2d", "school", "market", "market+school", "d2d", "market+school", "clinic", "clinic", "central.location"),
         unique.pot = c(NA, NA, NA, "Lwandeti MKT", NA, NA, NA, NA, NA, NA),
         d2d = mda.type == "d2d")

process.village <- function(.data) {
  vill.matches <- str_match(.data$Village, "(.+?)(\\s+(\\w2?))?$")
  
  .data %>% 
    mutate(village.name = vill.matches[, 2],
           village.code = vill.matches[, 4])
}

fill.missing.from.before <- function(.col, na = "") {
  if (length(.col) > 1) {
    for (col.index in 2:length(.col)) {
      if (.col[col.index] %in% na) {
        .col[col.index] <- .col[col.index - 1]
      }
    }
  } 
  
  return(.col)
}

process.filename <- function(.data, drop.unmatched = TRUE) {
    fn.matches <- str_match(.data$filename, "^Cluster\\s*(\\d{1,2})_([^\\(]+)(\\(.+\\))?_((\\d{1,2})-(\\d{1,2})(-(\\d{2,4}))?)")
    
    .data %<>% filter(!is.na(fn.matches[, 1]))
    fn.matches <- fn.matches[!is.na(fn.matches[, 1]), ]
   
    .data %>%  
      mutate(cluster = fn.matches[, 2],
             filename.pot = fn.matches[, 3],
             origin.deworming.date = fn.matches[, 4], 
             deworming.date = sprintf("2016-3-%d", pmax(as.numeric(fn.matches[, 6]), as.numeric(fn.matches[, 7]))) %>% ymd, 
             deworming.day = difftime(deworming.date, ymd("2016-3-7"), units = "days") + 1)
}

read.forms <- function(file.pattern, col.types = NULL, replace.names = NULL) 
  dir("~/Data/TakeUp/pilot/", file.pattern, full.names = TRUE) %>% 
    map(~ read_csv(., 
                   na = c("", "NA", "--blank--", "-", "--impossible--"), 
                   col_types = col.types)) %>% 
    bind_rows() %>% {
      if (!is.null(replace.names)) {
        plyr::rename(., replace.names)
      } else return(.)
    } %>% 
    separate(name, into = c("filename", "row.num"), sep = ": row") %>% 
    separate(filename, into = c("filename", "page.num"), sep = "\\s+p") %>% 
    mutate_each(funs(tolower(.)), starts_with("Name")) %>% 
    mutate_each(funs(as.numeric(.)), row.num, page.num) %>% 
    mutate(Gender = factor(toupper(str_sub(Gender, 1, 1)), levels = c("M", "F")),
           phone.num = str_replace(phone.num, "\\D", ""))

form.2a.data <- read.forms("Form2A_03-19.+.csv", 
                           col.types = cols(`Phone Numbers` = col_character(),
                                            Age = col_character()),
                           replace.names = c(Names = "Name",
                                             `Phone Numbers` = "phone.num",
                                             Village = "village.code")) %>% 
  process.filename %>% 
  mutate_each(funs(as.numeric(.)), cluster, deworming.day) %>% 
  filter(!is.na(Name)) %>%
  mutate(Age = str_extract(Age, "\\d{2}") %>% as.integer(),
         village.code = ifelse(village.code == "1.00E+06", "IE6", village.code),
         valid.code = aaply(village.code, 1, function(code) code %in% village.codes$village.code),
         code.length = aaply(village.code, 1, nchar),
         village.code = ifelse(valid.code, village.code, ifelse(code.length < 2, NA, village.code))) %>% 
  partial.match.village.code("village.code", village.codes, max.dist = 2) %>% 
  left_join(village.data[, c("village.code", "village.pop", "Cluster", "lat", "lon")], by = "village.code") %>% 
  rename(resp.cluster = Cluster,
         village.lat = lat,
         village.lon = lon) %>% 
  mutate_each(funs(ifelse(nchar(.) == 0, NA, .)), village.code) %>% 
  mutate(valid.code = aaply(village.code, 1, function(code) code %in% village.codes$village.code),
         targeted.village = !is.na(resp.cluster) & (resp.cluster == cluster)) %>% 
  select(-c(code.length, old.village.code)) %>% 
  distinct(cluster, Name, Age, Gender, phone.num) %>% 
  mutate(id = 1:n()) %>% 
  left_join(clust.pop.data, by=c("cluster" = "Cluster")) %>% 
  left_join(form.status.data %>% select(filename, PoT), by = "filename") %>% {
    mutate(., 
           PoT = ifelse(is.na(PoT), unique.pot, PoT),
           PoT = ifelse(is.na(PoT), filename.pot, PoT))
  } %>% 
  left_join(form.2a.fn.dict, by = c("filename" = "old.filename")) %>% 
  mutate(filename = ifelse(!is.na(new.filename), new.filename, filename),
         PoT = ifelse(!is.na(new.filename), PoT.y, PoT.x)) %>% 
  select(-new.filename, -PoT.x, -PoT.y) %>% 
  left_join(pot.geoloc.data %>% select(PoT, lat, lon), by = "PoT") %>% 
  rename(pot.lat = lat,
         pot.lon = lon) %>% 
  left_join(village.data[, c("Cluster", "village.code", "min.dist.pot")], by = c("cluster" = "Cluster", "village.code")) %>% 
  mutate(cluster.targeted.village = !is.na(min.dist.pot)) %>% 
  select(-min.dist.pot)
  

form.2c.data <- read.forms("Form2C_03-16+.csv",
                           replace.names = c(`Mobile Number` = "phone.num")) %>% 
  filter(!is.na(Name)) %>% 
  mutate(Village = str_replace_all(Village, "'|\"|(BLANK)|(\\)\\))", "") %>% 
           str_replace("^\\d+$", "") %>%
           str_replace_all("\\.|-", " ") %>% 
           str_replace_all("\\s+", " ") %>% 
           str_to_upper %>% 
           str_trim) %>% 
  mutate_each(funs(ifelse(nchar(.) == 0, NA, .)), Village) %>% 
  partial.match.village.code("Village", village.codes, max.dist = 2) %>% 
  mutate_each(funs(ifelse(nchar(.) == 0, NA, .)), village.code) %>% 
  group_by(filename, page.num) %>% 
  mutate(village.code = ifelse(length(unique(village.code)) == 2, 
                               table(village.code) %>% magrittr::extract(which.max(.)) %>% names, 
                               village.code)) %>% 
  ungroup %>% 
  left_join(village.data %>% select(village.code, Cluster), by = "village.code") %>% 
  group_by(Cluster) %>% 
  mutate(min.name.dist = adist(Name, form.2a.data$Name[form.2a.data$cluster == first(Cluster) & !is.na(form.2a.data$Name)]) %>% 
           alply(1, min) %>% 
           unlist, 
         min.phone.dist = ifelse(is.na(phone.num), 
                                 NA,
                                 adist(phone.num, 
                                       form.2a.data$phone.num[form.2a.data$cluster == first(Cluster) & !is.na(form.2a.data$phone.num)]) %>% 
                                   alply(1, min) %>% 
                                   unlist)) %>% 
  ungroup %>% 
  distinct(Cluster, Name, Gender, phone.num) %>% 
  mutate(id = 1:n())

# find.closest.2a.match <- function(2c.data, 2a.data, max.distance, join.by, 2c.col, 2a.col = 2c.col) {
#   left_join(2a.data, 2c.data, by = join.by)
# }

save(form.2a.data, form.2c.data, village.codes, village.data, clust.pop.data, file = "pilot_forms.RData")
