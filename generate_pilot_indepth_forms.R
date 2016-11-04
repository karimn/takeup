library(readr)
library(plyr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(magrittr)
library(doParallel)

tryCatch({
  config <- yaml::yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores)
}, error=function(err) {
  registerDoSEQ()
})

load("pilot_forms.RData")

deworm.lists <- form.2c.data %>% 
  filter(!is.na(Cluster)) %>% 
  group_by(village.code) %>% 
  filter(min.name.dist <= 1 | min.phone.dist <= 2)

no.deworm.lists <- form.2c.data %>% 
  group_by(village.code) %>% 
  filter(min.name.dist > 4, min.phone.dist > 4 | is.na(min.phone.dist))

deworm.lists %>% 
  count(Cluster, village.code) %>% 
  rename(n.deworm = n) %>% 
  left_join(no.deworm.lists %>% 
              count(Cluster, village.code) %>% 
              rename(n.no.deworm = n)) %>% 
  as.data.frame

indepth.villages <- c("MY9", "NA3", "CH5", "RS6", "MG5", "EB7", "KW1", "EM1", "MK1", "MB2")

all.deworm.lists <- deworm.lists %>% 
  mutate(dewormed = TRUE) %>% 
  bind_rows(no.deworm.lists %>% 
              mutate(dewormed = FALSE))

generate.social.know.table <- function(num.vill.interviews = 8, num.sample.per.grp = 5, knit.list.files = TRUE, individ.interview.callback = NULL) {
  foreach(vill = indepth.villages) %do% {
    vill.deworm.list <-  all.deworm.lists %>% filter(dewormed, village.code == vill) %>% sample_n(nrow(.))
    vill.no.deworm.list <- all.deworm.lists %>% filter(!dewormed, village.code == vill) %>% sample_n(nrow(.))
    
    list.env <- new.env()
    list.env$cluster <-  vill.deworm.list %$% Cluster %>% first
    list.env$village.code <- vill
    list.env$village.name <- village.codes %>% filter(village.code == vill) %$% Village %>% first
    list.env$deworm.names <- vill.deworm.list %$% Name %>% str_to_upper 
    list.env$deworm.ids <- vill.deworm.list %$% id 
    list.env$no.deworm.names <- vill.no.deworm.list %$% Name %>% str_to_upper
    list.env$no.deworm.ids <- vill.no.deworm.list %$% id 
    
    md.filename <- sprintf("individ_list_clust%d_%s.md", list.env$cluster, vill)
   
    if (knit.list.files) {
      knitr::knit("individ_cluster_list.Rmd", md.filename, envir = list.env)
      
      knitr::pandoc(md.filename, format = "latex")
    }
    
    if(empty(vill.deworm.list)) return(NULL)

    stopifnot(!empty(vill.no.deworm.list))
    
    foreach(interview.index = seq_len(num.vill.interviews)) %do% {
      table.env <- new.env()
      table.env$cluster <- list.env$cluster
      table.env$village.code <- vill
      table.env$village.name <- list.env$village.name
      table.env$table.num <- interview.index
      
      sampled <- vill.deworm.list %>% sample_n(num.sample.per.grp) %>% 
        bind_rows(vill.no.deworm.list %>% sample_n(num.sample.per.grp)) %>% 
        sample_n(nrow(.))
      
      table.env$obs.names <- sampled$Name %>% str_to_upper
      table.env$ids <- sampled$id
      table.env$dewormed <- sampled$dewormed
      
      table.md.filename <- sprintf("observability_clust%d_%s_%d.md", table.env$cluster, vill, interview.index)
     
      if (knit.list.files) {
        knitr::knit("observability_table.Rmd", 
                    table.md.filename,
                    envir = table.env)
        
        knitr::pandoc(table.md.filename, format = "latex")
      }
      
      if (!is.null(individ.interview.callback)) {
        individ.interview.callback(table.env)
      }
    }
  }
}

soc.know.csv.factory <- function() {
  csv.data <- NULL
  
  list(callback = function(table.env) {
    invisible(csv.data <<- bind_rows(csv.data, with(table.env, data.frame(cluster, village.code, village.name, ids, obs.names, dewormed))))
  },
  print = function(filename = paste0(tempfile(), ".csv")) {
    write_csv(csv.data, filename)
    
    invisible(filename)
  })
}

rct.instr.pilot.data.manager <- soc.know.csv.factory()
generate.social.know.table(num.vill.interviews = 1, knit.list.files = FALSE, individ.interview.callback = rct.instr.pilot.data.manager$callback)

#save(all.deworm.lists, file = "all.deworm.RData")