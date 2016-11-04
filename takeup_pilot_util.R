kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

partial.match.village.code <- function(.data, col.name, village.codes.data, max.dist = Inf, costs = NULL) {
  new.village.code <- adist(.data[[col.name]], 
                              village.codes.data[[col.name]], 
                              ignore.case = TRUE, 
                              costs = costs) %>% 
    alply(1, function(dist.row) { 
      if(all(is.na(dist.row)) || min(dist.row) > max.dist) {
        return("")
      } else {
        village.codes.data$village.code[which.min(dist.row)] 
      }
    }) %>% 
    unlist
  
  .data %<>% 
    mutate(targeted.village = nchar(new.village.code) > 0)

  if ("village.code" %in% names(.data)) {
    .data %<>% 
      mutate(old.village.code = village.code,
             village.code = ifelse(targeted.village, 
                                   new.village.code, 
                                   sub(sprintf(".*(%s).*", paste(village.codes.data$village.code, collapse = "|")), "\\1", old.village.code)))
  } else {
    .data %<>%
      mutate(village.code = new.village.code)
  }
  
  return(.data)
}