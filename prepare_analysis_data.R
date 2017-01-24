library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(lubridate)

source("analysis_util.R")

# This data was prepared in takeup_field_notebook.Rmd
census.data <- read_rds(file.path("data", "takeup_census.rds"))
baseline.data <- read_rds(file.path("data", "takeup_baseline_data.rds"))
takeup.data <- read_rds(file.path("data", "takeup.rds"))
all.endline.data <- read_rds(file.path("data", "all_endline.rds"))
reconsent.data <- read_rds(file.path("data", "reconsent.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))
sms.content.data <- read_rds(file.path("data", "takeup_sms_treatment.rds"))

endline.data <- prepare.endline.data(all.endline.data)
consent.dewormed.reports <- prepare.consent.dewormed.data(all.endline.data, reconsent.data)
analysis.data <- prepare.analysis.data(census.data, takeup.data, endline.data, consent.dewormed.reports, cluster.strat.data)

old.baseline.data <- baseline.data
baseline.data %<>% prepare.baseline.data

cluster.takeup.data <- prepare.cluster.takeup.data(analysis.data)

save(all.endline.data, endline.data, consent.dewormed.reports, analysis.data, baseline.data, cluster.takeup.data, 
     census.data, reconsent.data, takeup.data, sms.content.data,
     file = file.path("data", "analysis.RData"))