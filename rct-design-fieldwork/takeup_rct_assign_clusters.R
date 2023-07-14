buffer.clusters <- function(cluster.points.data, .width, proj4 = kenya.proj4) {
  cluster.points.data %>%
    magrittr::extract(!duplicated(coordinates(.)), ) %>% 
    spTransform(CRS(proj4)) %>% # Transform in order to calculate distances and areas
    gBuffer(byid = TRUE, width = .width, joinStyle = "ROUND")
}

generate.cluster.pop.area.tester <- function(min.pop.prop = 2, school.buffer.radius = 500) { 
  school.pop.area <- pi * (school.buffer.radius^2)
  
  function(clusters, ...) {
    if (is.null(clusters) || empty(clusters)) return(FALSE)
    
    gArea(clusters, byid = TRUE) >= min.pop.prop * school.pop.area
  }
}

generate.cluster.num.schools.tester <- function(primary.schools, min.num.schools = 2) {
  inner.primary.schools <- buffer.clusters(primary.schools, 1000) 
  
  function(clusters, ...) {
    contain.matrix <- gContains(clusters, inner.primary.schools, byid = TRUE) 
    
    colSums(contain.matrix) %>% 
      magrittr::is_weakly_greater_than(min.num.schools)
  }
}

generate.cluster.num.schools.tester2 <- function(primary.schools, min.num.schools = 1, min.area.frac = 0.5, school.buffer.radius = 500) {
  inner.primary.schools <- buffer.clusters(primary.schools, school.buffer.radius) 
  school.pop.area <- pi * (school.buffer.radius^2)
  
  function(cluster, cluster.schools.ids, ...) {
    length(cluster.schools.ids) > 0 && (inner.primary.schools %>% 
                                                   magrittr::extract(.$cluster.id %in% cluster.schools.ids & str_to_upper(.$county) == str_to_upper(cluster$county), ) %>% {
                                                     if (!empty(.)) gIntersection(., cluster, byid = TRUE)
                                                   } %>% {
                                                     !is.null(.) && (gArea(., byid = TRUE) %>% 
                                                                       is_weakly_greater_than(min.area.frac * school.pop.area) %>% 
                                                                       sum %>% 
                                                                       is_weakly_greater_than(min.num.schools))
        })
    }
}

generate.cluster.schools.pop.tester <- function(primary.schools, min.pop.prop = 2, school.buffer.radius = 500) {
  inner.primary.schools <- buffer.clusters(primary.schools, school.buffer.radius) %>% 
    gUnaryUnion
  
  school.pop.area <- pi * (school.buffer.radius^2)
  
  function(clusters, ...) {
    gIntersection(clusters, inner.primary.schools, byid = TRUE) %>% {
      if (is.null(.)) {
        return(FALSE)
      }
      
      gArea(., byid = TRUE) %>% 
        magrittr::is_weakly_greater_than(min.pop.prop * school.pop.area)
    }
  }
}

ggplot.clusters <- function(selected.clusters, 
                            .rct.schools.data = rct.schools.data, 
                            pilot.locations = NULL, 
                            proj4 = kenya.proj4, 
                            maptype = "roadmap", 
                            source = "google",
                            include.cluster.ids = TRUE, 
                            suppress.selected.clusters = FALSE,
                            caption = NULL,
                            ...) {
  all.clusters <- gUnaryUnion(selected.clusters) 
 
  if (!is.null(.rct.schools.data)) {
    targetable.schools <- .rct.schools.data %>% 
      spTransform(proj4) %>% 
      gWithin(all.clusters, byid = TRUE) %>% 
      drop
  } else {
    targetable.schools <- NULL
  }
  
  all.clusters %<>% spTransform(wgs.84)
  
  map.obj <- all.clusters %>% 
    tidy %>% 
    make_bbox(long, lat, data = .) %>% 
      # get_map(maptype = maptype, source = source, ...) %>% 
      # get_openstreetmap(...) %>% 
      get_stamenmap(maptype = maptype, ...) %>% 
      ggmap()

  if (!suppress.selected.clusters) {
    cluster_polygons <- selected.clusters |>
      spTransform(wgs.84) |>  
      sf::st_as_sf()
    
    map.obj <- map.obj + geom_sf(aes(color = county), alpha = 0.5, linetype = "dashed", inherit.aes = FALSE, data = cluster_polygons)
  }
  
  cluster_center <- spTransform(rct.schools.data, wgs.84)
  
  if (!suppress.selected.clusters) {
    cluster_center %<>% magrittr::extract(., .$cluster.id %in% selected.clusters$cluster.id,)
    
    if (include.cluster.ids) {
      map.obj <- map.obj + geom_text(aes(lon, lat, label = cluster.id), data = as.data.frame(cluster_center))
    } else {
      map.obj <- map.obj + geom_sf(inherit.aes = FALSE, shape = 3, data = sf::st_as_sf(cluster_center)) 
    }
  } else {
    map.obj <- map.obj + geom_sf(inherit.aes = FALSE, shape = 3, color = alpha("darkred", 0.5), data = sf::st_as_sf(cluster_center))
  }
  
  map.obj <- map.obj + geom_point(aes(lon, lat), shape = 3, size = 1, data = pilot.locations) +
    labs(x = "", y = "", caption = caption) +
    scale_color_discrete("County") +
    theme(legend.position = "bottom")
  
  return(map.obj)
}

plot.clusters <- function(selected.clusters, selected.buffers, cluster.ids = FALSE, school.radius = 500) {
  all.clusters <- gUnaryUnion(selected.clusters)
  
  rct.schools.data %>% 
    spTransform(CRS(kenya.proj4)) %>% {
      # buffer.clusters(., school.radius) %>% {
        # magrittr::extract(., !.$cluster.id %in% selected.clusters$cluster.id, ) %>% {
          magrittr::extract(., !drop(gWithin(., all.clusters, byid = TRUE)), ) %>% 
            buffer.clusters(school.radius) %>% 
            plot(border="grey")
          
          magrittr::extract(., drop(gWithin(., all.clusters, byid = TRUE)), ) %>%
            buffer.clusters(school.radius) %>% 
            plot(add = TRUE, border="blue")
      #   }
      # }
      
      magrittr::extract(., .$cluster.id %in% selected.clusters$cluster.id, ) %>% { 
        if (!cluster.ids) plot(., add = TRUE, col = "red")
        if (cluster.ids) text(coordinates(.), labels = .$cluster.id, col = alpha("blue", 0.8), lwd = 2)
      } 
    }

  spTransform(subcounties.adm.data, CRS(kenya.proj4)) %>% 
    gUnaryUnion %>% 
    gIntersection(selected.clusters, byid = TRUE) %>% 
    plot(add = TRUE, col=alpha("red", 0.2))
  
  if (!missing(selected.buffers)) {
    selected.buffers %>% 
      plot(add = TRUE, lty = "dotted", border = alpha("black", 0.5))
  }
}

tu.plot.clusters <- function(yy, ...) {
  yy %>%
    magrittr::extract(.$selected, ) %>%
    plot.clusters(...)
  
  spTransform(rct.schools.data, CRS(kenya.proj4)) %>% {
    schools.data <- .
    purrr::pmap(list(c("bracelet.airtime", "control.ink"), 
                     # c(4000, 2500),
                     attr(yy, "ca.outer.radius")[c("bracelet.airtime", "control.ink")],
                     c("darkgreen", "red")),
                function(grp, radius, color) { 
                  try(schools.data[schools.data$cluster.id %in% yy$cluster.id[yy$cluster.group == grp], ] %>%
                        buffer.clusters(radius) %>%
                        plot(add = TRUE, lty = "dotted", border = alpha(color, 1)), silent = TRUE) 
                })
  }
  
  invisible(NULL)
}

get.rct.clusters <- function(sp.pot.data, rct.area, schools.data, 
                             num.clusters = NULL, # c(control = 30, social.incentive = 60, airtime = 30)
                             ca.outer.radius = 5000, # c(control = 2000, social.incentive = 4000, airtime = 5000) 
                             ca.inner.radius = ca.outer.radius, # c(control = 2000, social.incentive = 2000, airtime = 4000) 
                             cluster.group.order = c("random", "desc.outer.radius"),
                             #school.radius = 500,
                             cluster.size.tester = function(cluster, ...) cluster %>% gArea(byid = TRUE) %>% is_weakly_greater_than(pi * (ca.inner.radius^2) * 0.5),
                             proj4 = kenya.proj4,
                             verify.shapes = FALSE,
                             plot.rct.clusters = FALSE) {
  # Unique IDs
  if (is.null(sp.pot.data$cluster.id)) {
    sp.pot.data$cluster.id <- factor(seq_len(nrow(sp.pot.data)))
  }
  
  rand.cluster.group.order <- switch(match.arg(cluster.group.order), random = TRUE, desc.outer.radius = FALSE)
  
  stopifnot(length(dplyr::setdiff(names(ca.outer.radius), names(num.clusters))) == 0)
  stopifnot(all(ca.outer.radius >= ca.inner.radius))
  
  ca.outer.radius %<>% sort(decreasing = TRUE)
  num.clusters  <- num.clusters[names(ca.outer.radius)]
  
  create.buffers <- . %>% 
    purrr::map(~ buffer.clusters(sp.pot.data, ., proj4)) %>% 
    purrr::map(function(buffers) { rownames(buffers@data) <- buffers$cluster.id; buffers })
  
  outer.rct.buffers <- create.buffers(ca.outer.radius)

  inner.rct.buffers <- create.buffers(ca.inner.radius)[[1]]  
  
  polygon.ids <- . %>% `@`(polygons) %>% purrr::map(~ .@ID) %>% unlist
  
  # Inner function for possible splitting by (sub)counties 
  inner.get.rct.clusters <- function(inner.buffers, outer.buffers) {
    
    if (plot.rct.clusters && !is.null(rct.area)) {
      rct.area %>% spTransform(CRS(kenya.proj4)) %>% plot
    }
    
    prep.inner.buffers <- function(.inner.buffers) {
      # avail.area <- buffer.clusters(schools.data, school.radius) %>% 
      #   gUnaryUnion %>% 
      #   gIntersection(rct.area) 
      
      if (!is.null(rct.area)) {
        avail.area <- rct.area
        
        .inner.buffers <- .inner.buffers[c(gIntersects(avail.area, .inner.buffers, byid = TRUE)), ]
        
        .inner.buffers@polygons <- gIntersection(.inner.buffers, avail.area, byid = TRUE)@polygons 
        
        stopifnot(gIsValid(.inner.buffers))
      }
      
      .inner.buffers$selected <- FALSE
      .inner.buffers$cluster.group <- NA
      
      return(.inner.buffers)
    }
   
    inner.buffers %<>% prep.inner.buffers 
    
    rct.clusters.history <- NULL
    dropped.clusters.history <- NULL
    cluster.index <- num.clusters
    
    inner.buffers.intersect.mat <- t(gIntersects(outer.buffers[[1]], inner.buffers, byid = TRUE))
    inner.buffers.intersect.mat <- inner.buffers.intersect.mat[outer.buffers[[1]]$cluster.id %in% inner.buffers$cluster.id, ]
    colnames(inner.buffers.intersect.mat) <- inner.buffers$cluster.id
    rownames(inner.buffers.intersect.mat) <- inner.buffers$cluster.id
    diag(inner.buffers.intersect.mat) <- FALSE
   
    if (!is.null(schools.data)) {
      contained.schools.mat <- schools.data %>% 
        spTransform(CRS(kenya.proj4)) %>% # Transform in order to calculate distances and areas
        gIntersects(inner.buffers, byid = TRUE) %>% 
        magrittr::extract(, schools.data$cluster.id %in% inner.buffers$cluster.id)
      
      colnames(contained.schools.mat) <- inner.buffers$cluster.id
      # colnames(contained.schools.mat) <- schools.data$cluster.id
      rownames(contained.schools.mat) <- inner.buffers$cluster.id
    } else {
      contained.schools.mat <- NULL
    }
    
    while (is.null(num.clusters) || (sum(cluster.index) > 0)) {
     
      if (is.null(cluster.index)) {
        current.cluster.group <- NULL
      } else if (rand.cluster.group.order) {
        current.cluster.group <- sample(length(cluster.index), 1, prob = cluster.index * ca.outer.radius)
      } else {
        current.cluster.group <- min(which(cluster.index > 0)) 
      }
      
      new.rct.cluster <- inner.buffers %>% 
        magrittr::extract(!.$selected, ) %>% 
        use_series(cluster.id) %>% {
          if (length(.) > 1) {
            sample(., 1)
          } else if (length(.) == 1) {
            return(.)
          } else {
            return(NULL)
          }
        }
      
      if (is.null(new.rct.cluster)) {
        break
      }
    
      new.cluster <- inner.buffers %>% 
        magrittr::extract(.$cluster.id == new.rct.cluster, )
      
      if (!cluster.size.tester(new.cluster, colnames(contained.schools.mat)[contained.schools.mat[new.rct.cluster,]])) {
        inner.buffers %<>% magrittr::extract(.$cluster.id != new.rct.cluster, )
        
        dropped.clusters.history %<>% c(., new.rct.cluster)  
        
        next  
      }
      
      original.new.cluster <- outer.buffers[[if (!is.null(current.cluster.group)) current.cluster.group else 1]] %>% 
        magrittr::extract(.$cluster.id == new.rct.cluster, )
      
      stopifnot(!is.null(original.new.cluster))
    
      intersected.inner.buffers <- colnames(inner.buffers.intersect.mat)[inner.buffers.intersect.mat[new.rct.cluster, ]] %>% 
        dplyr::intersect(inner.buffers$cluster.id) #%>% 
        # dplyr::setdiff(new.rct.cluster)
      
      if (any(inner.buffers$cluster.id %in% intersected.inner.buffers)) {
        old.inner.polygons <- inner.buffers@polygons[inner.buffers$cluster.id %in% intersected.inner.buffers]
        
        saved.intersect.ids <- Kmisc::duplicate(intersected.inner.buffers)
        
        differenced.polygons <- inner.buffers[inner.buffers$cluster.id %in% intersected.inner.buffers, ] %>% 
          gDifference(original.new.cluster, byid = TRUE, id = saved.intersect.ids, checkValidity = verify.shapes)  
        
        if (!is.null(differenced.polygons)) {
          differenced.polygon.ids <- polygon.ids(differenced.polygons)
          inner.buffers@polygons[inner.buffers$cluster.id %in% differenced.polygon.ids] <- differenced.polygons@polygons
          
          buffers.to.remove <- dplyr::setdiff(intersected.inner.buffers, differenced.polygon.ids)
        } else {
          buffers.to.remove <- intersected.inner.buffers
        }
      }
      
      any.rct.population.lost <- inner.buffers %>% 
        magrittr::extract(.$selected & .$cluster.id %in% intersected.inner.buffers, ) %>% {
          for (intersected.index in seq_len(nrow(.))) {
            if (!cluster.size.tester(.[intersected.index, ], colnames(contained.schools.mat)[inner.buffers.intersect.mat[.$cluster.id[intersected.index],]])) {
              return(TRUE)
            }
          } 
          
          return(FALSE)
        } 
    
      # Check if the new cluster would disqualify any of the already selected clusters 
      if (any(rct.clusters.history %in% buffers.to.remove) || any.rct.population.lost) { 
        inner.buffers@polygons[inner.buffers$cluster.id %in% intersected.inner.buffers] <- old.inner.polygons
        
        # Only exclude if there isn't a potentially smaller cluster in other smaller out radius cluster groups
        if (is.null(cluster.index) || max(which(cluster.index > 0)) <= current.cluster.group) { 
          inner.buffers %<>% magrittr::extract(.$cluster.id != new.rct.cluster, )
          
          dropped.clusters.history %<>% c(., new.rct.cluster)  
        }
        
        next  
      } 
      
      inner.buffers$selected[inner.buffers$cluster.id == new.rct.cluster] <- TRUE
      
      inner.buffers <- inner.buffers[!inner.buffers$cluster.id %in% buffers.to.remove, ]
      dropped.clusters.history %<>% c(., buffers.to.remove)
      
      if (!is.null(current.cluster.group)) {
        inner.buffers$cluster.group[inner.buffers$cluster.id == new.rct.cluster] <- names(cluster.index)[current.cluster.group]
        cluster.index[current.cluster.group] %<>% subtract(1)
      }
      
      if (plot.rct.clusters) {
        inner.buffers %>% magrittr::extract(.$selected, ) %>% plot(add = TRUE, col = scales::alpha("grey", 0.2))
        plot(new.cluster, add = TRUE, border = "red", col = scales::alpha("red", 0.2))
      }
     
      rct.clusters.history %<>% c(., new.rct.cluster)

    }
    
    stopifnot(!anyDuplicated(rct.clusters.history))
    
    attr(inner.buffers, "ca.inner.radius") <- ca.inner.radius
    attr(inner.buffers, "ca.outer.radius") <- ca.outer.radius
    attr(inner.buffers, "cluster.group.order") <- cluster.group.order
    # attr(inner.buffers, "school.radius") <- school.radius
 
    return(inner.buffers) 
  }
  
  inner.get.rct.clusters(inner.rct.buffers, outer.rct.buffers)
}

get.cluster.villages.data <- function(rct.clusters, school.buffer.radius = 1000, min.area.frac = 0.6, .rct.schools.data = rct.schools.data) {
  school.area <- pi * (school.buffer.radius^2)
  
  all.rct.clusters <- gUnaryUnion(rct.clusters)
  
  schools.within <- .rct.schools.data %>% 
    spTransform(CRS(kenya.proj4)) %>% 
    gWithin(rct.clusters, byid = TRUE) %>% 
    t %>% 
    set_rownames(.rct.schools.data$cluster.id) %>% 
    set_colnames(rct.clusters$cluster.id) %>% 
    magrittr::extract(rowSums(.) > 0, )
  
  pot.pts <- .rct.schools.data %>% 
    magrittr::extract(.$cluster.id %in% rct.clusters$cluster.id, ) %>% 
    spTransform(CRS(kenya.proj4)) 
 
  school.regions <- .rct.schools.data %>% 
    magrittr::extract(.$cluster.id %in% rownames(schools.within), ) %>% 
    buffer.clusters(school.buffer.radius) %>% 
    gIntersection(all.rct.clusters, byid = TRUE, id = .$cluster.id) 
  
  ret.data <- school.regions %>% 
    gArea(byid = TRUE) %>% 
    magrittr::is_weakly_greater_than(min.area.frac * school.area) %>% 
    magrittr::extract(., .) %>%
    names %>% 
    magrittr::extract(schools.within, ., ) %>% 
    adply(1, . %>% { data.frame(pot.cluster.id = names(.)[.]) }, .id = "targeted.cluster.id") %>% 
    left_join(rct.clusters@data[, c("cluster.id", "cluster.group", "county")], by = c("pot.cluster.id" = "cluster.id")) %>% 
    left_join(.rct.schools.data@data[, c("cluster.id", "county")], by = c("targeted.cluster.id" = "cluster.id")) %>% 
    filter(county.x == county.y) %>% 
    select(-starts_with("county"))
  
  dist.mat <- spTransform(.rct.schools.data, kenya.proj4) %>% 
    magrittr::extract(.$cluster.id %in% ret.data$targeted.cluster.id,) %>% { 
      set_rownames(set_colnames(gDistance(., pot.pts, byid = T), .$cluster.id), pot.pts$cluster.id)
    }
  
  ret.data %>% 
    group_by(pot.cluster.id) %>% 
    mutate(dist = dist.mat[first(pot.cluster.id), as.character(targeted.cluster.id)],
           village.dist.cat = (ifelse(dist <= (2500/2), "close", "far")),
           cluster.dist.cat = ifelse(all(village.dist.cat == "close"), "close", ifelse(all(village.dist.cat == "far"), "far", "mixed"))) %>% 
    ungroup
}