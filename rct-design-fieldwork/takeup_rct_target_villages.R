
# Identify valid target villages and clusters -----------------------------

rct.cluster.selection <- read_rds("rct_cluster_selection_2.0.rds")
rct.targetable.schools <- read_rds("rct_targetable_schools_2.0.rds")

clusters.to.drop <- c("99", "691", "892", "1293", "239") %>% # Problems with clusters at cluster survey stage; dropped from study
  c("201", "853", "1402") # These are are study cluster that are no longer valid based on distance to other PoT

rct.cluster.selection %<>% magrittr::extract(!.$cluster.id %in% clusters.to.drop, )

rct.cluster.selection@data %<>% left_join(rct.schools.data %>% 
                                            spTransform(wgs.84) %>% 
                                            as.data.frame %>% 
                                            select(cluster.id, lon, lat)) %>% 
  rename(pot.lon = lon, pot.lat = lat)

rct.cluster.selection %<>% 
  `@`(data) %>% 
  left_join(cluster.survey.data %>% 
              distinct(cluster.id, found, num.valid.villages, gpslatitude, gpslongitude, gps2latitude, gps2longitude) %>% 
              filter(!duplicated(cluster.id)),
            by = "cluster.id") %>% 
  mutate(alt.pot.lat = gps2latitude,
         alt.pot.lon = gps2longitude,
         pot.lat = ifelse(!is.na(gpslatitude), #& found %in% c("Found", "Name Changed"), 
                          gpslatitude, 
                          ifelse(!is.na(gps2latitude), #& found %in% c("Closed", "Not Found"), 
                                 alt.pot.lat, 
                                 pot.lat)),
         pot.lon = ifelse(!is.na(gpslongitude), #& found %in% c("Found", "Name Changed"), 
                          gpslongitude, 
                          ifelse(!is.na(gps2longitude), #& found %in% c("Closed", "Not Found"), 
                                 alt.pot.lon, 
                                 pot.lon))) %>%
  as.data.frame %>%
  convert.to.sp(~ pot.lon + pot.lat, wgs.84) %>% 
  spTransform(kenya.proj4)

cluster.survey.data %<>%
  filter(!is.na(target.lat), !is.na(target.lon)) %>% 
  arrange(cluster.id, target.village.id) %>% 
  left_join(filter(., !is.na(target.lon), !is.na(target.lat)) %>% 
              ddply(.(cluster.id), function(cluster.villages) {
                convert.to.sp(cluster.villages, ~ target.lon + target.lat, wgs.84) %>%
                  spTransform(kenya.proj4) %>%  
                  gDistance(rct.cluster.selection, byid = TRUE) %>% 
                  set_colnames(cluster.villages$target.village.id) %>% 
                  magrittr::extract(, !duplicated(colnames(.)), drop = FALSE) %>% 
                  tibble::as_tibble() %>% 
                  mutate(dist.from.cluster.id = rct.cluster.selection$cluster.id,
                         cluster.id = cluster.villages$cluster.id[1],
                         cluster.group = rct.cluster.selection$cluster.group) %>% 
                  tidyr::gather(target.village.id, dist, -c(dist.from.cluster.id, cluster.id, cluster.group)) %>% 
                  left_join(., 
                            select(., dist.from.cluster.id, target.village.id, dist), 
                            by = c("cluster.id" = "dist.from.cluster.id", "target.village.id"),
                            suffix = c(".to.other.pot", ".to.pot")) %>%  
                  filter(dist.from.cluster.id != cluster.id) %>% 
                  select(-dist.from.cluster.id) %>% 
                  group_by(target.village.id, cluster.group) %>% 
                  filter(row_number(dist.to.other.pot) == 1) %>% 
                  ungroup %>%
                  spread(cluster.group, dist.to.other.pot) %>% 
                  mutate(valid.target.village = bracelet.airtime >= 3860 & control.ink >= 3000,
                         target.village.id = as.integer(target.village.id)) 
              }),
            by = c("cluster.id", "target.village.id")) %>% 
  left_join(filter(rct.targetable.schools, selected.targeted) %>% 
              select(pot.cluster.id, village.dist.cat, cluster.dist.cat),
            by = c("cluster.id" = "pot.cluster.id")) %>% 
  group_by(cluster.id) %>% {
    if (first(.$village.dist.cat) == "far") arrange(., desc(dist.to.pot)) else arrange(., dist.to.pot)
  } %>% 
  filter(!duplicated(target.village_name)) %>%
  ungroup %>%
  mutate(target.actual.village.dist.cat = ifelse(dist.to.pot <= (2500/2), "close", "far"),
         target.valid.dist.to.pot = target.actual.village.dist.cat == village.dist.cat,
         target.in.range = dist.to.pot <= 2500,
         dist.switchable = !target.valid.dist.to.pot & valid.target.village,
         valid.target.village = valid.target.village & target.in.range, 
         target.village.pop.size = ifelse(target.pop_households <= median(target.pop_households, na.rm = TRUE), "small", "big")) %>% 
  filter(!is.na(target.village.pop.size)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster.pop.size.strata = unique(target.village.pop.size) %>% { ifelse(length(.) == 2, "mixed", .) },
         actual.cluster.dist.cat = unique(target.actual.village.dist.cat) %>% { ifelse(length(.) == 2, "mixed", .) }) %>% 
  ungroup

rct.cluster.selection@data %<>% left_join(cluster.survey.data %>% 
                                            mutate(backup.cluster = FALSE) %>% # cluster.id %in% rct.backup.clusters$cluster.id) %>% 
                                            group_by(cluster.id, backup.cluster, cluster.pop.size.strata, actual.cluster.dist.cat) %>% 
                                            summarize(any.buffer.valid.village = any(valid.target.village),
                                                      dist.switchable = any(dist.switchable)) %>% 
                                            ungroup,
                                          by = "cluster.id") %>% 
  left_join(rct.schools.data@data %>% 
              filter(selected.targeted) %>% 
              select(pot.cluster.id, village.dist.cat, cluster.dist.cat),
            by = c("cluster.id" = "pot.cluster.id")) 

rct.cluster.selection@data %<>%
  mutate(village.dist.cat = ifelse(dist.switchable & !backup.cluster, ifelse(village.dist.cat == "close", "far", "close"), village.dist.cat),
         any.buffer.valid.village = (dist.switchable & !backup.cluster) | any.buffer.valid.village)


# Display RCT clusters ----------------------------------------------------
 
rct.cluster.selection@data %>%
  filter(!any.buffer.valid.village) %>%
  select(cluster.id, cluster.group) %>%
  left_join(cluster.survey.data) %>%
  select(cluster.id, matches("valid"), bracelet.airtime, control.ink, dist.to.pot, village.dist.cat, cluster.group, dist.switchable)

# Village selection -------------------------------------------------------

rct.villages <- cluster.survey.data %>% 
  filter(valid.target.village, cluster.id %in% rct.cluster.selection$cluster.id) %>% 
  group_by(cluster.id) %>% 
  sample_n(1) %>% 
  ungroup
 
# write_rds(rct.villages, "rct_target_villages_2.0.rds") 
rct.villages <- read_rds("rct_target_villages_2.0.rds") 

unused.villages <- cluster.survey.data %>% anti_join(rct.villages, by = c("cluster.id", "target.village.id"))

# reassign.pop.size.strata <- . %>% {
#   if (first(.$cluster.pop.size.strata) == "mixed" & first(.$actual.cluster.dist.cat != "mixed")) {
#     stopifnot(first(.$missing.big.pop.size) <= nrow(.))
#     
#     new.strata.assignment <- sample(c(rep("big", first(.$missing.big.pop.size)), rep("small", nrow(.) - first(.$missing.big.pop.size))))
#   
#     mutate(., assigned.pop.size.strata = new.strata.assignment)
#   } else return(.)
# }
# 
# reassign.dist.cat <- . %>% {
#   if (first(.$cluster.pop.size.strata) != "mixed" & first(.$actual.cluster.dist.cat == "mixed")) {
#     stopifnot(first(.$missing.far.size) <= nrow(.))
#     
#     new.dist.assignment <- sample(c(rep("far", first(.$missing.far.size)), rep("close", nrow(.) - first(.$missing.far.size))))
#   
#     mutate(., assigned.cluster.dist.cat = new.dist.assignment)
#   } else return(.)
# }
# 
# reassign.mixed.mixed <- . %>% {
#   
# }
#   
# rct.cluster.selection@data %>% 
#   filter(any.buffer.valid.village) %>% 
#   # Assign strata to "mixed" clusters
#   mutate(assigned.pop.size.strata = cluster.pop.size.strata,
#          assigned.cluster.dist.cat = actual.cluster.dist.cat) %>% 
#   group_by(cluster.group) %>% 
#   mutate(cluster.group.size = n()) %>% 
#   group_by(actual.cluster.dist.cat, add = TRUE) %>% 
#   mutate(assigned.big.pop.size = n() %/% 2 + sample(seq(0, n() %% 2), 1),
#          missing.big.pop.size = max(0, first(assigned.big.pop.size) - sum(cluster.pop.size.strata == "big"))) %>% 
#   group_by(cluster.pop.size.strata, add = TRUE) %>% 
#   do(reassign.pop.size.strata(.)) %>% 
#   group_by(cluster.group, cluster.pop.size.strata) %>%
#   mutate(assigned.far.size = n() %/% 2 + sample(seq(0, n() %% 2), 1),
#          missing.far.size = max(0, first(assigned.far.size) - sum(actual.cluster.dist.cat == "far"))) %>%
#   group_by(actual.cluster.dist.cat, add = TRUE) %>%
#   do(reassign.dist.cat(.)) %>%
#   group_by(cluster.group) %>% 
#   ungroup %>% 
#   count(cluster.group, assigned.cluster.dist.cat, assigned.pop.size.strata)

# RCT target villages info table ------------------------------------------

write.rct.village.info <- function(.data, .file.name) {
  .data %>%
    select(cluster.id, target.village_name, target.lon, target.lat, 
           matches("boda"), matches("target.*elder"), matches("target.*chv")) %>% 
    left_join(rct.cluster.selection %>% 
                tidy %>% 
                select(cluster.id, school.name, county, constituency, ward, pot.lon, pot.lat),
              by = "cluster.id") %>%
    arrange(as.integer(cluster.id)) %>% 
    write_csv(.file.name, na = "")
  
}

rct.villages %>% write.rct.village.info("takeup_rct_villages.csv")

# Plot RCT clusters and villages ------------------------------------------

rct.cluster.selection %>%
  magrittr::extract(.$any.buffer.valid.village, ) %>%
  spTransform(wgs.84) %>% 
  tidy %>%
  transmute(lon = pot.lon, lat = pot.lat, type = "PoT", cluster.id) %>% 
  bind_rows(rct.villages %>% 
              transmute(lon = target.lon, lat = target.lat, type = "Target Village", cluster.id)) %>% {
    spatial.joint.data <- filter(., type == "PoT") %>%
      convert.to.sp(~ lon + lat, wgs.84) %>% {
      rbind(mutate(tidy(spTransform(buffer.clusters(., .width = 3000), wgs.84)), cluster.group = "control.ink"), 
            mutate(tidy(spTransform(buffer.clusters(., .width = 4000), wgs.84)), cluster.group = "bracelet.calendar"))
    } %>% 
      mutate(id = paste(cluster.group, id, sep = "-"))
                            
    ggmap(get_map(make_bbox(lon, lat, data = .), maptype = "toner")) +
      geom_point(aes(lon, lat, color = type), shape = 1, stroke = 1.5, alpha = 0.5, data = .) +
      geom_line(aes(lon, lat, group = cluster.id), data = .) +
      geom_polygon(aes(long, lat, group = id, fill = cluster.group), alpha = 0.25, data = spatial.joint.data) +
      labs(x = "", y = "") +
      scale_color_manual("", values = c("red", "darkgreen")) +
      theme(legend.position = "bottom")
              }

# Intercluster distance ---------------------------------------------------

rct.villages %>%
  left_join(rct.cluster.selection@data %>% select(cluster.id, cluster.group)) %>% 
  group_by(cluster.group) %>% 
  do(dist.mat = convert.to.sp(., ~ target.lon + target.lat, wgs.84) %>% 
       spTransform(kenya.proj4) %>% 
       gDistance(rct.cluster.selection, byid = TRUE))
  