
# Clean data --------------------------------------------------------------

rct.cluster.selection <- read_rds("rct_cluster_selection_2.0.rds")

rct.cluster.selection@data %<>% left_join(rct.schools.data %>% 
                                            spTransform(wgs.84) %>% 
                                            as.data.frame %>% 
                                            select(cluster.id, lon, lat)) %>% 
  rename(pot.lon = lon, pot.lat = lat)

read.cluster.survey.data <- function(data.filename) {
  switch(tools::file_ext(data.filename),
         csv = read_csv(data.filename),
         dta = read_dta(data.filename)) %>% {
    if (!any(grepl("^anchor_gps", names(.)))) mutate(., anchor_gpslatitude = NA, anchor_gpslongitude = NA) else .
  } %>% {
    if (!"correct_ward" %in% names(.)) mutate(., correct_ward = NA, correction = NA) else .
  } %>% #{
  #   if (is.character(.$surveydate)) {
  #     mutate(., surveydate = sub("^(\\d+)/(\\d+)/(\\d+).*", "\\3-\\2-\\1", surveydate) %>% 
  #              sub("^(\\d+)jul2016", "2016-07-\\1", .)) 
  #   } else return(.) 
  # } %>% 
    rename(cluster.id = clusterid) %>% 
    select(cluster.id, found, surveydate, new_name, manual_lat, manual_long,
           correct_ward, correction, # Ward updates
           boda, boda_contact, alt_boda, alt_boda_contact,
           matches("^gps[2-3]?\\.?l(at|on)"), starts_with("gps_work"),
           alt_name, location_type, other_type,
           # starts_with("market_days"),
           starts_with("anchor_gps"),
           matches("^(village_name|(second_)?(elder|chv)\\d?(_(contact|phone2|2nd_phone|2nd_contact|name))?|manual_.+|(?#anchor_(big|local)_market_[1-7]|)pop_households|pop_individuals)_?\\d"),
           matches("^pop_(households|individuals)_?\\d")) %>% 
    rename(anchor.gps.work = gps_work2,
           anchor.manual.lat = manual_lat2,
           anchor.manual.lon = manual_long2) %>%
    set_names(sub("^anchor_(.+\\d)$", "target.\\1", names(.))) %>% 
    set_names(sub("(.+)(\\d)(\\d)$", "\\1\\2_\\3", names(.))) %>% 
    set_names(sub("(village_name|gps3.+[^_]|check|pop_.*[^_])_?(\\d)$", "\\1_\\2", names(.))) %>% 
    set_names(ifelse(grepl("^(?!market|target).+_?\\d{1,2}$", names(.), perl = TRUE), paste0("target.", names(.)), names(.))) %>% 
    gather(key, val, starts_with("target.")) %>% 
    tidyr::extract(key, c("var", "target.village.id"), "(.+?)_?(\\d+)$") %>% 
    spread(var, val) %>% 
    filter(!is.na(target.village_name), nchar(target.village_name) > 0) %>% 
    mutate_at(vars(matches("target.+(long|lat|pop_(households|individuals))")), funs(as.numeric)) %>% 
    mutate(data.filename = data.filename,
           surveydate = lubridate::as_date(surveydate),
           corrected.ward = ifelse(correct_ward == 1, correction, NA),
           corrected.pot.school.name = ifelse(nchar(new_name) > 0, new_name, NA)) %>% 
    select(-c(correct_ward, correction, new_name)) %>% 
    rename(alt.pot.name = alt_name,
           alt.pot.location.type = location_type,
           alt.pot.location.type.other = other_type)
}

cluster.survey.data <- ldply(c("Cluster Survey V1.dta", "Cluster Survey 27.07.2016.dta"), read.cluster.survey.data) %>% 
  group_by(cluster.id) %>% 
  mutate(target.village.id = seq_len(n())) %>% 
  ungroup %>% 
  mutate(cluster.id = as.character(cluster.id),
         target.lat = ifelse(target.gps_work3 == 1 | target.gps_work3 == "Yes", target.gps3latitude, target.manual_lat3),
         target.lon = ifelse(target.gps_work3 == 1 | target.gps_work3 == "Yes", target.gps3longitude, target.manual_long3),
         found = factor(found, 1:4, labels = c("Found", "Name Changed", "Closed", "Not Found"))) %>% 
  left_join(rct.cluster.selection@data[, c("cluster.id", "pot.lon", "pot.lat")], by = "cluster.id") 

cluster.survey.data %<>% filter(!is.na(cluster.id))

cluster.survey.data %>% 
  check_that(gps_work > 0 | (!is.na(manual_lat) & !is.na(manual_long)),
             gps_work > 0 | (manual_lat == manual_lat_check & manual_long == manual_long_check),
             anchor.gps.work > 0 | (!is.na(anchor.manual.lat) & !is.na(anchor.manual.lon)),
             anchor.gps.work > 0 | (anchor.manual.lat == manual_lat2_check & anchor.manual.lon == manual_long2_check),
             target.gps_work3 > 0 | (!is.na(target.manual_lat3) & !is.na(target.manual_long3)),
             target.gps_work3 > 0 | (target.manual_lat3 == target.manual_lat3_check & target.manual_long3 == target.manual_long3_check)) %>% 
  summary

cluster.survey.data %>% 
  count(cluster.id, target.village.id) %>%
  ungroup %>% 
  check_that(n == 1) %>% 
  summary

anchor.buffers <- rct.schools.data %>% 
  magrittr::extract(!is.na(.$selected.targeted) & .$selected.targeted, ) %>% 
  buffer.clusters(.width = 1000)

cluster.survey.data %<>% left_join(filter(., !is.na(target.lon), !is.na(target.lat)) %>% 
                                     convert.to.sp(~ target.lon + target.lat, wgs.84) %>% 
                                     spTransform(kenya.proj4) %>% {
                                       transmute(.@data, 
                                                 cluster.id, 
                                                 target.village.id, 
                                                 target.anchor.zone.cluster.id = gWithin(., anchor.buffers, returnDense = FALSE, byid = TRUE) %>% 
                                                   laply(. %>% { ifelse(is.null(.), NA, anchor.buffers$pot.cluster.id[.]) }),
                                                 target.in.anchor.zone = !is.na(target.anchor.zone.cluster.id) & cluster.id == target.anchor.zone.cluster.id)
                                     }, by = c("cluster.id", "target.village.id")) %>% 
  group_by(cluster.id) %>% 
  mutate(num.valid.villages = sum(target.in.anchor.zone, na.rm = TRUE)) %>% 
  ungroup

#write_rds(cluster.survey.data, "takeup_cluster_survey.rds")

# Survey status plot -----------------------------------------------------------

anchor.buffers %>% 
  spTransform(wgs.84) %>% 
  tidy(region = "pot.cluster.id") %>% 
  left_join(distinct(cluster.survey.data, cluster.id, num.valid.villages), by = c("id" = "cluster.id")) %>% 
  mutate(cluster.survey.status = case_when(is.na(.$num.valid.villages) ~ 1, #"Not Surveyed",
                                           .$num.valid.villages > 0 ~ 2, #"Surveyed - Complete",
                                           .$num.valid.villages == 0 ~ 3) %>%  #"Surveyed - Incomplete") 
           factor(levels = 1:3, labels = c("Not Surveyed", "Surveyed - Complete", "Surveyed - Incomplete"))) %>% {
    ggmap(get_map(make_bbox(long, lat, data = .), maptype = "toner")) +
      geom_polygon(aes(long, lat, group = group, fill = cluster.survey.status), color = "black", alpha = 0.5, data = .) +
      labs(x = "", y = "") +
      scale_fill_manual("", values = c("red", "green", "yellow"), drop = FALSE) +
      theme(legend.position = "bottom")
 }

# Plot individual clusters ------------------------------------------------
  
plot.cluster.survey.locations <- function(cluster.locations, .maptype = "terrain", .zoom = 14) {
  onekm.buffers <- anchor.buffers %>% 
    magrittr::extract(.$pot.cluster.id %in% cluster.locations$cluster.id, ) %>% 
    spTransform(wgs.84) %>% 
    tidy
  
  ggmap(get_googlemap(center = c(mean(cluster.locations$pot.lon), mean(cluster.locations$pot.lat)), maptype = .maptype, zoom = .zoom, scale = 2, key = config$google_api_key)) +
    geom_polygon(aes(x = long, y = lat, group = group), linetype = "dotted", alpha = 0, color = "darkgreen", size = 1, data = onekm.buffers) +
    geom_point(aes(x = pot.lon, y = pot.lat), shape = 1, color = "red", size = 3, stroke = 1, data = cluster.locations %>% distinct(cluster.id, .keep_all = TRUE)) +
    geom_point(aes(x = gpslongitude, y = gpslatitude), shape = 4, color = "red", size = 4, stroke = 1, data = cluster.locations %>% distinct(cluster.id, .keep_all = TRUE)) +
    geom_point(aes(x = gps2longitude, y = gps2latitude), shape = 2, color = "red", size = 3, stroke = 1, data = cluster.locations %>% distinct(cluster.id, .keep_all = TRUE)) +
    geom_point(aes(x = anchor_gpslongitude, y = anchor_gpslatitude), shape = 4, color = "darkgreen", size = 3, stroke = 1, data = cluster.locations %>% distinct(cluster.id, .keep_all = TRUE)) +
    geom_point(aes(x = target.lon, y = target.lat), shape = 3, color = "darkgreen", size = 3, stroke = 1, 
               data = cluster.locations %>% filter(target.in.anchor.zone)) +
    geom_point(aes(x = target.lon, y = target.lat), shape = 3, color = "darkgrey", size = 3, stroke = 1, 
               data = cluster.locations %>% filter(is.na(target.anchor.zone.cluster.id))) +
    geom_point(aes(x = target.lon, y = target.lat), shape = 3, color = "blue", size = 3, stroke = 1, 
               data = cluster.locations %>% filter(!is.na(target.anchor.zone.cluster.id) & cluster.id != target.anchor.zone.cluster.id)) +
    labs(x = "", y = "") +
    ggtitle(sprintf("Cluster%s %s (%s)", 
                    ifelse(length(unique(cluster.locations$cluster.id)) > 1, "s", ""), 
                    paste(unique(cluster.locations$cluster.id), collapse = ", "),
                    paste(as.character(unique(cluster.locations$surveydate)), collapse = ", "))) 
}

# missing.cluster.survey.plots <- cluster.survey.data %>% 
cluster.survey.plots <- cluster.survey.data %>%
  arrange(surveydate) %>% 
  # filter(num.valid.villages < 1) %>%
  # filter(cluster.id %in% missing.clusters) %>% 
  dlply(.(cluster.id), .parallel = FALSE,  plot.cluster.survey.locations) 

cluster.survey.plots %>%
# missing.cluster.survey.plots %>%  
  marrangeGrob(ncol = 2, nrow = 2) %>% 
  ggsave("cluster_survey.pdf", ., scale = 2)
  # ggsave("cluster_survey_incomplete.pdf", ., scale = 2)

# Cluster stats -----------------------------------------------------------

cluster.survey.data %>% 
  count(cluster.id, num.valid.villages) %>% 
  arrange(as.numeric(cluster.id)) %>% 
  rename(total.villages = n) %>% 
  ungroup %>% 
  as.data.frame

# Village Population ------------------------------------------------------

cluster.survey.data %>% 
  filter(target.in.anchor.zone) %>% { 
    ggplot(.) +
      geom_histogram(aes(as.numeric(target.pop_households)),  alpha = 0.75, binwidth = 10) + 
      geom_vline(xintercept = median(as.numeric(.$target.pop_households), na.rm = TRUE), color = "red", linetype = "dashed") + 
      scale_x_continuous(breaks = seq(0, 1000, 100)) +
      coord_cartesian(xlim = c(0, 1000)) + 
      labs(x = "Village Households", y = "")
  # geom_histogram(aes(as.numeric(target.pop_individuals)), binwidth = 10)
  }

cluster.survey.data %>% 
  filter(target.in.anchor.zone) %>% { 
    ggplot(.) +
      geom_histogram(aes(as.numeric(target.pop_individuals)),  alpha = 0.75, binwidth = 50) + 
      geom_vline(xintercept = median(as.numeric(.$target.pop_individuals), na.rm = TRUE), color = "red", linetype = "dashed") + 
      scale_x_continuous(breaks = seq(0, 1000, 500)) +
      coord_cartesian(xlim = c(0, 2000)) + 
      labs(x = "Village Individuals", y = "")
  # geom_histogram(aes(as.numeric(target.pop_individuals)), binwidth = 10)
  }

cluster.survey.data %>% 
  filter(target.in.anchor.zone) %>% 
  mutate(target.individ.hh.ratio = as.numeric(target.pop_individuals)/as.numeric(target.pop_households)) %>% {
    ggplot(.) +
      geom_histogram(aes(target.individ.hh.ratio),  alpha = 0.75, binwidth = 50) + 
      geom_vline(xintercept = median(.$target.individ.hh.ratio, na.rm = TRUE), color = "red", linetype = "dashed") + 
      # scale_x_continuous(breaks = seq(0, 1000, 500)) +
      # coord_cartesian(xlim = c(0, 2000)) + 
      labs(x = "Village Individuals/HH", y = "")
  # geom_histogram(aes(as.numeric(target.pop_individuals)), binwidth = 10)
  }

# Alt Pot data ------------------------------------------------------------

rct.villages %>% 
  arrange(as.numeric(cluster.id)) %>% 
  select(cluster.id, alt.pot.name, alt.pot.location.type, alt.pot.location.type.other) %>% 
  mutate(alt.pot.location.type = factor(alt.pot.location.type, levels = 1:5, labels = c("Clinic", "Church", "Market", "Home", "Other"))) %T>% 
  write_csv("alt.pot.csv") %T>% {
    print(ggplot(.) + geom_bar(aes(alt.pot.location.type)))
  } %T>% {
    filter(., nchar(alt.pot.location.type.other) > 0) %>% {
      ggplot(.) + 
        geom_bar(aes(alt.pot.location.type.other)) +
        theme(axis.text.x = element_text(angle=90, vjust = 0.5))
    } %>% 
      print
  }