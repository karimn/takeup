
library(magrittr)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggmap)
library(ggstance)
library(gridExtra)
library(cowplot)
library(rgeos)
library(sp)
library(knitr)
library(modelr)
library(car)
library(rstan)
library(latex2exp)
library(ggthemes)

library(sf)
library(econometr)
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))

fit_version <- 66
default_top_levels = c("Bracelet", "Combined")

# 66 ed fit
# 60 Karim fit
# 62 also Karim fit


model_fit_by = if_else(fit_version %in% c(60, 62), "Karim", "Ed")


models_we_want = c(
  "STRUCTURAL_LINEAR_U_SHOCKS"
)


quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

output_basepath = str_glue("temp-data/output_dist_fit{fit_version}")

canva_palette_vibrant <- "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))

load(file.path("data", "takeup_village_pot_dist.RData"))

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data

cluster.strat.data %>%
  colnames()




## Fit Loading

load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))

dist_fit_data %<>%
  left_join(tribble(
    ~ fit_type,        ~ model_color,
      "fit",           "black", 
      "prior-predict", "darkgrey",
  ), by = "fit_type")

rf_analysis_data <- dist_fit_data %>% 
  filter(
    fct_match(model_type, "reduced form"),
    fct_match(fit_type, "fit"),
  ) %$% 
  stan_data[[1]]$analysis_data 
delta <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))





takeup_data = dist_fit_data %>%
  filter( 
    fct_match(fit_type, "fit") &
    fct_match(model, 'STRUCTURAL_LINEAR_U_SHOCKS')
  ) %>%
  select(cluster_takeup_prop) %>%
  unnest(cluster_takeup_prop) 

takeup_data %>%
  filter(assigned_treatment == "control") %>%
  select(roc_distance, mean_est) %>%
  unique()

control_data = takeup_data %>%
  filter(assigned_treatment == "control")
  

struct_control_spline = 
  splinefun(
    x = control_data$roc_distance,
    y = control_data$mean_est, 
    method = "natural"
  )

bracelet_data = takeup_data %>%
  filter(assigned_treatment == "bracelet") 
struct_bracelet_spline = 
  splinefun(
    x = bracelet_data$roc_distance,
    y = bracelet_data$mean_est, 
    method = "natural"
  )


lm_fit = analysis_data %>%
  lm(
    data = ., 
    dewormed ~ cluster.dist.to.pot*factor(assigned.treatment) + factor(assigned.treatment)*cluster.dist.to.pot^2
  )



approx_df = tibble(
  roc_distance = seq(from = 0, to = 5000, by = 100)
) %>% 
  mutate(takeup_approx = struct_bracelet_spline(roc_distance))

p = approx_df %>%
    ggplot(aes(
        x = roc_distance, 
        y = takeup_approx
    )) +
    geom_point()


## Grabbing Data


cluster_ids_actually_used = analysis_data %>%
  select(cluster.id) %>%
  pull() %>%
  unique()

rct_cluster_df = st_as_sf(rct.cluster.selection) 

set.seed(100)

rct_pot_cluster_ids = rct_cluster_df %>%
  filter(selected == TRUE) %>%
  filter(cluster.id %in% cluster_ids_actually_used) %>%
  pull(cluster.id)

n_orig_pots = length(rct_pot_cluster_ids)

rct_school_df = st_as_sf(rct.schools.data)

unselected_pots = setdiff(
  rct_school_df$cluster.id,
  rct_pot_cluster_ids)

# additional_pots = sample(unselected_pots, round(n_orig_pots*1.1), replace = FALSE)
additional_pots = NULL



pot_df = rct_school_df %>%
  filter(cluster.id %in% c(additional_pots, rct_pot_cluster_ids)) %>%
  mutate(id = 1:n())


village_df = village.centers %>%
  st_as_sf(coords = c("lon", "lat"), crs = wgs.84) %>%
  mutate(id = 1:n())



distance_matrix = st_distance(
  village_df,
  pot_df
)


long_distance_mat = distance_matrix %>%
    as_tibble() %>%
    mutate(index_i = 1:n()) %>%
    gather(variable, dist, -index_i) %>%
    mutate(index_j = str_extract(variable, "\\d+") %>% as.numeric())  %>%
    select(index_i, index_j, dist) %>%
    mutate(
        dist = as.numeric(dist),
        dist_km = as.numeric(dist/1000))


create_approx_demand = function(dist_m, demand_func) {
  demand_df = dist_m %>%
    mutate(
        pred_takeup = demand_func(dist), 
        pred_takeup = pmax(pred_takeup, 0),
        pred_takeup = replace_na(pred_takeup, 0)
    ) %>%
    rename(village_i = index_i, pot_j = index_j, demand = pred_takeup)
  return(demand_df)
}


bracelet_demand_df = create_approx_demand(long_distance_mat, struct_bracelet_spline)
control_demand_df = create_approx_demand(long_distance_mat, struct_control_spline)

lm_preds = long_distance_mat %>%
  mutate(
    assigned.treatment = "bracelet", 
    cluster.dist.to.pot = dist
  ) %>%
  predict(lm_fit, newdata = .)

reduced_bracelet_demand_df = long_distance_mat %>%
  mutate(pred_takeup = lm_preds, 
        pred_takeup = pmax(pred_takeup, 0),
        pred_takeup = replace_na(pred_takeup, 0)) %>%
  rename(
    village_i = index_i, 
    pot_j = index_j, 
    demand = pred_takeup) %>%
  group_by(village_i) %>%
  mutate(closest_pot = dist == min(dist)) %>%
  ungroup()


reduced_bracelet_demand_df %>%
    write_csv("optim/data/approx-reduced-bracelet-demand.csv")
control_demand_df %>%
    write_csv("optim/data/approx-control-demand.csv")
bracelet_demand_df %>%
    write_csv("optim/data/approx-bracelet-demand.csv")


pot_df %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]) %>%
    as_tibble() %>%
    write_csv("optim/data/pot-df.csv")

village_df %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]) %>%
    as_tibble() %>%
    write_csv("optim/data/village-df.csv")
