#!/usr/bin/Rscript
script_options = docopt::docopt(
  stringr::str_glue("Usage:
  balance.R [options]
  
Options:
  --num-boot-draws=<num-boot-draws>  Number of bootstrap draws to use in wild sub-cluster bootstrap [default: 1000]
  --output-path=<path>  Where to save output files [default: {file.path('temp-data')}]
  --cts-interval=<cts-interval>  Interval for continuous distance binning [default: 200]
"),
  args = if (interactive()) "
    --num-boot-draws=999 \
    --output-path=temp-data \
    --cts-interval=200
    " else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "takeup cv --models=REDUCED_FORM_NO_RESTRICT --cmdstanr --include-paths=stan_models --update --output-path=data/stan_analysis_data --outputname=test --folds=2 --sequential" else commandArgs(trailingOnly = TRUE)

) 


set.seed(12932)

library(tidyverse)
library(marginaleffects)
library(fwildclusterboot) # I hope you have Julia installed Ed
library(broom)
library(knitr)
library(kableExtra)
library(ggthemes)
library(fixest)
library(magrittr)


script_options$num_boot_draws = as.numeric(script_options$num_boot_draws)
script_options$cts_interval = as.numeric(script_options$cts_interval)


source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))

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

#| balance-setup
baseline.data = baseline.data %>%
  mutate(
    baseline_neighbours_worm_knowledge = case_when(
      neighbours_worms_affect == "yes" ~ TRUE, 
      neighbours_worms_affect == "no" ~ FALSE
    )
  )


# Create floor quality variable
baseline.data = baseline.data %>%
  mutate(
    floor_tile_cement = floor == "Cement" | floor == "Tiles"
  ) %>%
  group_by(
    cluster.id
  ) %>%
  mutate(
    frac_floor_tile_cement = mean(floor_tile_cement, na.rm = TRUE)
  ) %>%
  ungroup()

# Create years of schooling and completley primary variables
baseline.data = baseline.data %>%
  mutate(
    completed_primary = (school == "Primary 8" | str_detect(school, "Secondary|College|University"))
  ) %>%
  mutate(
    schooling_years_plus = case_when(
      str_detect(school, "Primary") ~ 0, 
      str_detect(school, "Secondary") ~ 8, 
      str_detect(school, "College") ~ 16, 
      str_detect(school, "University") ~ 16
    ), 
    digits_schooling = str_extract(school, "\\d+") %>% as.numeric(), 
    years_schooling = digits_schooling + schooling_years_plus, 
    years_schooling = if_else(school == "Never gone to school", 0, years_schooling), 
    years_schooling = if_else(str_detect(school, "College|University"), 16, years_schooling)
  ) %>%
  select(-digits_schooling, schooling_years_plus)

# Months since an individual took deworming treatment at baseline (i.e. independent of the campaign)
baseline.data = baseline.data %>%
  mutate(
    treated_digit = str_extract(treated_when, "\\d+") %>% as.numeric, 
    treated_months = case_when(
      str_detect(treated_when, "year") ~ 12, 
      str_detect(treated_when, "mon") ~ 1, 
      TRUE ~ NA_real_
    )
    ) %>%
    mutate(
      months_since_treatment = treated_digit*treated_months
    ) %>%
    select(-treated_digit, -treated_months) %>%
    # has someone been dewormed in the last 12 months
    mutate(
        dewormed_last_12 = case_when(
            str_detect(treated_when, "mon|(1 year)") ~ TRUE,
            is.na(treated_when) ~ NA,
            TRUE ~ FALSE
        )
    )

baseline.data = baseline.data %>%
  mutate(
    have_phone_lgl = case_when(
      have_phone == "Yes" ~ TRUE, 
      have_phone == "No" ~ FALSE, 
      TRUE ~ NA
    ), 
    treated_lgl = case_when(
      treated == "yes" ~ TRUE, 
      treated == "no" ~ FALSE, 
      TRUE ~ NA
    ), 
  )
baseline.data = baseline.data %>%
  # these are nested lists of responses so we map_lgl and use any()
  mutate(
    all_can_get_worms = map_lgl(who_worms, ~any(str_detect(.x, "everyone") | (str_detect(.x, "adult") & str_detect(.x, "child")))), 
    correct_when_treat = map_lgl(when_treat, ~any(.x == "every 6 months")), 
    know_deworming_stops_worms = map_lgl(stop_worms, ~any(.x == "medicine"))
  ) 


# creating a single treat x distance variable for balance testing
cluster_treat_df = analysis_data %>%
  mutate(
      treat_dist = paste0(
      "treat: ", 
      assigned.treatment,
      ", dist: ", dist.pot.group
      ) %>% factor()
  ) %>%
  select(cluster.id, treat_dist) %>%
  unique()



know_balance_data = analysis_data %>%
  nest_join(
    endline.know.table.data %>% 
      filter(fct_match(know.table.type, "table.A")),
    by = "KEY.individ", 
    name = "knowledge_data"
  ) %>% 
  mutate(
    map_dfr(knowledge_data, ~ {
      tibble(
        obs_know_person = sum(.x$num.recognized),
        obs_know_person_prop = mean(.x$num.recognized),
        knows_other_dewormed = sum(fct_match(.x$dewormed, c("yes", "no")), na.rm = TRUE),
        knows_other_dewormed_yes = sum(fct_match(.x$dewormed, "yes"), na.rm = TRUE),
        thinks_other_knows = sum(fct_match(.x$second.order, c("yes", "no")), na.rm = TRUE),
        thinks_other_knows_yes = sum(fct_match(.x$second.order, "yes"), na.rm = TRUE),
      )
    }
  )) %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  )


know_vars = c("obs_know_person")
know_balance_fit = feols(
    data = know_balance_data, 
    .[know_vars] ~ 0 + treat_dist, 
    ~county
    ) 


## Baseline Balance
baseline_balance_data = baseline.data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  )

baseline_balance_data = baseline_balance_data %>%
  left_join(
    analysis_data %>%
      select(cluster.id, cluster.dist.to.pot) %>%
      unique(), 
    by = "cluster.id"
  ) 

baseline_vars = c(
  "years_schooling", 
  "know_deworming_stops_worms",
  "treated_lgl", 
  "dewormed_last_12", 
  "floor_tile_cement",
  "all_can_get_worms",
  "correct_when_treat",
  "baseline_neighbours_worm_knowledge"
)


baseline_balance_fit = feols(
    data = baseline_balance_data, 
    .[baseline_vars] ~ 0 + treat_dist, 
    ~county
    ) 

# PoT level balance variables
balance_variables = c(
  "cluster.dist.to.pot"
)
# Indiv level balance variables
indiv_balance_vars = c(
  "female", 
  "phone_owner", 
  "age"
)

rct_school_df = rct.schools.data %>% 
    as_tibble()
# Adding school (PoT) data to analysis df
analysis_school_data = left_join(
    analysis_data,
    rct_school_df %>% mutate(cluster.id = as.numeric(cluster.id)) ,
    by = "cluster.id"
)


analysis_school_data = analysis_school_data %>%
mutate(
    treat_dist = paste0(
    "treat: ", 
    assigned.treatment,
    ", dist: ", dist.pot.group
    ) %>% factor()
)  %>%
mutate(
  female = fct_match(gender, "female")
)

  
indiv_balance_fit = feols(
    data = analysis_school_data, 
    .[indiv_balance_vars] ~ 0 + treat_dist,
    cluster = ~county
    ) 

# Probably a better way to get the schools in the sample
school_treat_df = analysis_school_data %>%
  filter(!is.na(assigned.treatment)) %>%
  select(any_of(colnames(rct_school_df)), treat_dist, cluster.dist.to.pot, constituency, county) %>%
  unique()


school_balance_fit = feols(
    data = school_treat_df, 
    .[balance_variables] ~ 0 + treat_dist,
    ~county

  )

#### Endline
endline_balance_data = endline.data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  mutate(
    endline_neighbours_worm_knowledge = case_when(
      neighbours_worms_affect == "yes" ~ TRUE, 
      neighbours_worms_affect == "no" ~ FALSE, 
      TRUE ~ NA
    )
  ) %>%
  mutate(
    all_can_get_worms = map_lgl(who_worms, ~any(str_detect(.x, "everyone") | (str_detect(.x, "adult") & str_detect(.x, "child")))), 
    correct_when_treat = map_lgl(when_treat, ~any(.x == "every 6 months")), 
    know_deworming_stops_worms = map_lgl(stop_worms, ~any(.x == "medicine"))
  ) 
  


endline_and_baseline_data = bind_rows(
  endline_balance_data %>%
    select(
      treat_dist, 
      neighbours_worm_knowledge = endline_neighbours_worm_knowledge, 
      all_can_get_worms,
      correct_when_treat, 
      know_deworming_stops_worms,
      constituency, county) %>%
    mutate(
      type = "endline"
    ),
  baseline_balance_data %>%
    select(
      treat_dist, 
      neighbours_worm_knowledge = baseline_neighbours_worm_knowledge, 
      all_can_get_worms,
      correct_when_treat, 
      know_deworming_stops_worms,
      constituency, county) %>%
    mutate(
      type = "baseline"
    )
) %>%
  na.omit()

endline_vars = c(
  "endline_neighbours_worm_knowledge", 
  "all_can_get_worms", 
  "correct_when_treat", 
  "know_deworming_stops_worms"
  )
endline_balance_fit = feols(
    data = endline_balance_data, 
    .[endline_vars] ~ 0 + treat_dist, 
    cluster = ~county
    ) 

# put all the baseline balance fits into a list we can map over
balance_fits = c(
  indiv_balance_fit,
  list("lhs: cluster.dist.to.pot" = school_balance_fit),
  baseline_balance_fit, 
  list("lhs: num_recognised" = know_balance_fit)
)


create_balance_comparisons = function(fit) {
  comp_df = avg_comparisons(
    fit,
    variables = list("treat_dist" = "all")
    ) %>%
    as_tibble()

  subset_comp_df = comp_df %>%
    mutate(
      lhs_treatment = str_extract(contrast, "(?<=^treat: )\\w+"), 
      rhs_treatment = str_extract(contrast, "(?<=- treat: )\\w+"), 
      lhs_dist = str_extract(contrast, "(?<=dist: )\\w+"),
      rhs_dist = str_extract(contrast, "(?<=, dist: )\\w+$")
    ) %>%
    filter(
      lhs_dist == rhs_dist
    ) %>%
    filter(rhs_treatment == "control" | lhs_treatment == "control") %>%
    filter(lhs_treatment != rhs_treatment)  


    rhs_control_comp_df = subset_comp_df %>%
      filter(rhs_treatment == "control") 

    lhs_control_comp_df = subset_comp_df %>%
      filter(rhs_treatment != "control")

    lhs_control_comp_df = lhs_control_comp_df %>%
      mutate(
        new_estimate = estimate*-1, 
        new_statistic = statistic*-1, 
        new_conf.low = conf.high*-1,
        new_conf.high = conf.low*-1,
        new_lhs_treatment = rhs_treatment,
        new_rhs_treatment = lhs_treatment
      )  %>%
      mutate(
        estimate = new_estimate, 
        statistic = new_statistic,
        conf.low = new_conf.low,
        conf.high = new_conf.high, 
        lhs_treatment = new_lhs_treatment,
        rhs_treatment = new_rhs_treatment
      ) %>%
      select(-contains('new_'))

    rearranged_comp_df = bind_rows(
      lhs_control_comp_df, 
      rhs_control_comp_df
   ) %>%
   select(-contrast)


    control_mean_df = fit %>%
      tidy(conf.int = TRUE) %>%
      filter(str_detect(term, "control")) %>%
      mutate(
        lhs_treatment = "control", rhs_treatment = NA, 
        lhs_dist = if_else(str_detect(term, "close"), "close", "far"), 
        rhs_dist = lhs_dist
      ) %>%
      select(
        -term
      )

  rearranged_comp_df = rearranged_comp_df %>%
    bind_rows(
      control_mean_df
    )

    return(rearranged_comp_df)

}



comp_balance_tidy_df = balance_fits %>%
  map_dfr(
    create_balance_comparisons, 
    .id = "lhs"
  )  %>%
  mutate(
      lhs = str_remove(lhs, "lhs: ")
  ) %>%
  mutate(
    lhs = str_replace_all(lhs, "\\.", " ") %>% str_to_title()
  )



balance_tidy_df = balance_fits %>%
    map_dfr(tidy, .id = "lhs") %>%
    mutate(
        lhs = str_remove(lhs, "lhs: ")
    ) %>%
    select(
        lhs, term, estimate, std.error
    )   %>%
    mutate(
      lhs = str_replace_all(lhs, "\\.", " ") %>% str_to_title()
    )

comp_balance_tidy_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "comp_balance_tidy_df.csv"
        )
    )

balance_tidy_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "balance_tidy_df.csv"
        )
    )
#### Continuous Balance
#| cts-dist-balance
dist_balance_cts_fun = function(rhs_var) {
  ## Indiv
  indiv_balance_cts_fit = feols(
      data = analysis_school_data %>%
        mutate(dist_measure = {{ rhs_var }}/1000), 
      .[indiv_balance_vars] ~  dist_measure + dist_measure^2,
      cluster = ~county
      ) 

  ## Know
  know_balance_cts_fit = feols(
      data = know_balance_data %>%
        mutate(dist_measure = {{ rhs_var }}/1000), 
      .[know_vars] ~ dist_measure + dist_measure^2, 
      ~county
      ) 

  ## baseline
  baseline_balance_cts_fit = feols(
      data = baseline_balance_data %>%
        mutate(dist_measure = {{ rhs_var }}/1000), 
      .[baseline_vars] ~  dist_measure + dist_measure^2, 
      ~county
      ) 

    return(
      c(
        indiv_balance_cts_fit,
        list("lhs: num_recognised" = know_balance_cts_fit),
        baseline_balance_cts_fit
      )
    )
}

dist_balance_disc_fun = function(rhs_var, interval_length) {
  analysis_school_data = analysis_school_data %>%
    mutate(
      dist_measure = cut_interval({{ rhs_var }}, length = interval_length )
    )

  know_balance_data = know_balance_data %>%
    mutate(
      dist_measure = cut_interval({{ rhs_var }}, length = interval_length )
    )
  baseline_balance_data = baseline_balance_data %>%
    mutate(
      dist_measure = cut_interval({{ rhs_var }}, length = interval_length )
    )

  ## Indiv
  indiv_balance_disc_fit = feols(
      data = analysis_school_data,
      .[indiv_balance_vars] ~  0 + dist_measure,
      cluster = ~county
      ) 

  ## Know
  know_balance_disc_fit = feols(
      data = know_balance_data,
      .[know_vars] ~ 0 + dist_measure, 
      ~county
      ) 

  ## baseline
  baseline_balance_disc_fit = feols(
      data = baseline_balance_data,
      .[baseline_vars] ~ 0 + dist_measure, 
      ~county
      ) 

    return(
      c(
        indiv_balance_disc_fit,
        list("lhs: num_recognised" = know_balance_disc_fit),
        baseline_balance_disc_fit
      )
    )
}


fully_cts_dist_balance = dist_balance_cts_fun(cluster.dist.to.pot)
disc_dist_balance = dist_balance_disc_fun(cluster.dist.to.pot, interval_length = 500)


construct_joint_test_m = function(object) {
  n_coef = length(coef(object))
  diag_m = diag(n_coef - 1)
  neg_1_m = matrix(-1, nrow = n_coef - 1, ncol = 1)

  hyp_m = cbind(neg_1_m, diag_m)
  return(hyp_m)
}



#### Joint Tests ####
#| joint-tests

## We want to test for balance across all conditions and balance within distance condition
## I don't know how to do such a joint test in R easily so we setup the test matrix 
## manually for the wald test

# Number of dist groups x treatment
n_variables = 8
# matrix R for test 
hyp_matrix = cbind(
  matrix(-1, nrow = n_variables - 1, ncol = 1 ), 
  diag(x = 1, nrow = n_variables - 1, ncol = n_variables)[, 1:(n_variables - 1)]
)

zero_matrix = matrix(0, nrow = 3, ncol = n_variables - 1) 
part_hyp_matrix = zero_matrix
for (i in 1:3) {
  part_hyp_matrix[i, 2*i] = 1
}

hyp_matrix_close = cbind(
  matrix(-1, nrow = 3, ncol = 1), 
  part_hyp_matrix
)

hyp_matrix_far = cbind(
  matrix(0, nrow = 3, ncol = 1),
  matrix(-1, nrow = 3, ncol = 1), 
  part_hyp_matrix[, 1:(ncol(part_hyp_matrix) - 1)]
)


perform_balance_cluster_boot = function(data, var, joint_R, close_R, far_R) {
  subset_data = data %>%
    select(all_of(var), treat_dist, county, constituency) %>%
    na.omit()


  fml = as.formula(paste0(var, " ~ 0 + treat_dist"))
  fit = lm(
    formula = fml, 
    data = subset_data
  )
  joint_boot_output = mboottest(
    object = fit, 
    clustid = c("county", "constituency"), 
    bootcluster = "constituency", 
    B = script_options$num_boot_draws, 
    R = joint_R
  )
  close_boot_output = mboottest(
    object = fit, 
    clustid = c("county", "constituency"), 
    bootcluster = "constituency", 
    B = script_options$num_boot_draws, 
    R = close_R
  )
  far_boot_output = mboottest(
    object = fit, 
    clustid = c("county", "constituency"), 
    bootcluster = "constituency", 
    B = script_options$num_boot_draws, 
    R = far_R
  )
  return(lst(
    joint_pval = joint_boot_output$p_val, 
    close_pval = close_boot_output$p_val, 
    far_pval = far_boot_output$p_val))
}

indiv_boot = map(
  indiv_balance_vars, 
  ~perform_balance_cluster_boot(
    data = analysis_school_data, 
    var = .x, 
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
     )
)


know_boot = map(
  know_vars, 
  ~perform_balance_cluster_boot(
    data = know_balance_data,
    var = .x, 
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

baseline_boot = map(
  baseline_vars, 
  ~perform_balance_cluster_boot(
    data = baseline_balance_data,
    var = .x, 
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

school_boot = map(
  balance_variables,
  ~perform_balance_cluster_boot(
    data = school_treat_df,
    var = .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)

endline_boot = map(
  endline_vars,
  ~perform_balance_cluster_boot(
    data = endline_balance_data,
    var = .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)


boot_fits = c(
  indiv_boot, 
  school_boot, 
  baseline_boot, 
  know_boot
)


boot_fits %>%
    saveRDS(
        file.path(
            script_options$output_path,
            "wild_boot_balance_fits.rds"
        )
    )




#| baseline-learning

# Here we test if there's a difference in learning between baseline and endline 
# across a range of alternative mechanisms.
comp_endline_vars = endline_vars %>%
  str_remove(., "endline_")
comp_endline_vars = comp_endline_vars[comp_endline_vars != "know_deworming_stops_worms"]
baseline_endline_externality_fit = feols(
      data = endline_and_baseline_data, 
      .[comp_endline_vars] ~ 0 + treat_dist:type, 
      ~county
      ) 

perform_externality_cluster_boot = function(data, var, R) {
  subset_data = data %>%
    select(all_of(var), treat_dist, type, county, constituency) %>%
      na.omit()


  fml = as.formula(paste0(var, " ~ 0 + treat_dist:type"))
  fit = lm(
    formula = fml, 
    data = subset_data
  )

  if (is.matrix(R)) {
    joint_boot_output = mboottest(
      object = fit, 
      clustid = c("county", "constituency"), 
      bootcluster = "constituency", 
      B = script_options$num_boot_draws, 
      R = R
    )
  } else {
    joint_boot_output = boottest(
      object = fit, 
      clustid = c("county", "constituency"),
      param = names(coef(fit)), 
      bootcluster = "constituency", 
      B = script_options$num_boot_draws, 
      R = R
    )
  }
  return(joint_boot_output$p_val)
}



indiv_externality_comp_test = baseline_endline_externality_fit %>%
  map(
    ~comparisons(
      .x,
      variable = list(type = "reference"), 
      newdata = datagrid(
        treat_dist = unique(baseline_balance_data$treat_dist)
      )
    )
  )

indiv_externality_comp_df = tibble(
  lhs = names(indiv_externality_comp_test), 
  p.value = map(indiv_externality_comp_test, "p.value"), 
  treat_dist = list(unique(baseline_balance_data$treat_dist)) 
) %>%
  mutate(lhs = str_remove(lhs, "lhs: ")) %>%
    unnest(c(p.value, treat_dist))

#' Another way to generate the hypothesis matrix - slightly more general
generate_joint_externality_hyp_m = function(fit, treat_term, dist_term) {
  hyp_df = fit %>%
    tidy() %>%
    select(term) %>% 
    mutate(
      treat = str_extract(term, "(?<=treat: ).*(?=,)"),
      dist = str_extract(term, "(?<=dist: ).*(?=:)"), 
      type = str_extract(term, "(?<=type).*$")
    ) %>%
    mutate(
      val = 0,
      val = if_else(
        treat == treat_term &
        dist == dist_term &
        type == "baseline", 
        -1, 
        val
        ),
      val = if_else(
        treat == treat_term &
        dist == dist_term &
        type == "endline", 
        1, 
        val
        )
    )
  return(hyp_df$val)
}

dist_treat_grid = expand_grid(
  treat = c("bracelet", "calendar", "ink", "control"), 
  dist = c("close", "far")
) %>%
  arrange(dist)

externality_joint_hyp_matrix = map2(
  dist_treat_grid$treat, 
  dist_treat_grid$dist, 
  ~generate_joint_externality_hyp_m(fit = baseline_endline_externality_fit[[1]], .x, .y )
) %>%
  do.call(rbind, .)


joint_externality_p_value =  map_dbl(
  comp_endline_vars,
  ~perform_externality_cluster_boot(
    data = endline_and_baseline_data, 
    var = .x, 
    R = externality_joint_hyp_matrix
  )
)

gen_close_p_val = function(x){
  map(
      c(split(externality_joint_hyp_matrix[1:4, ], 1:4), list(externality_joint_hyp_matrix[1:4, ])),
      ~perform_externality_cluster_boot(endline_and_baseline_data, x, .x)
  )
}

indiv_close_externality_p_value =  map(
  comp_endline_vars, 
  gen_close_p_val
)

gen_far_p_val = function(x) {
  map(
    c(split(externality_joint_hyp_matrix[5:8, ], 1:4), list(externality_joint_hyp_matrix[5:8, ])),
    ~perform_externality_cluster_boot(endline_and_baseline_data,x,  .x)
)
}

indiv_far_externality_p_value = map(comp_endline_vars, gen_far_p_val)



endline_p_val_df = map(1:length(comp_endline_vars), ~c(comp_endline_vars[.x], indiv_close_externality_p_value[[.x]] %>% unlist(), indiv_far_externality_p_value[[.x]] %>% unlist(), joint_externality_p_value[[.x]])) %>%
  bind_rows()

treat_levels_c = c("control", "ink", "calendar", "bracelet")
treat_levels = c("ink", "calendar", "bracelet")
col_order = c(
  "lhs", 
  paste0(treat_levels_c, "_close"),
  "close_joint_p",
  paste0(treat_levels_c, "_far"),
  "far_joint_p",
  "joint_p"
)

colnames(endline_p_val_df) = col_order

endline_p_val_df = endline_p_val_df %>%
  mutate(across(c(everything(), -lhs), as.numeric)) %>%
  mutate(fit_type = "pval") 


endline_tidy_df = endline_balance_fit %>%
    map_dfr(tidy, .id = "lhs") %>%
    mutate(
        lhs = str_remove(lhs, "lhs: ")
    ) %>%
    select(
        lhs, term, estimate, std.error
    )   %>%
    mutate(
      lhs = str_replace_all(lhs, "\\.", " ") %>% str_to_title()
    ) 


endline_tidy_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "endline_balance_tidy_df.csv"
        )
    )

endline_p_val_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            "endline_balance_p_val_df.csv"
        )
    )

#### Continuous Distance Tests ####


fully_cts_dist_balance = dist_balance_cts_fun(cluster.dist.to.pot)


fully_cts_dist_balance %>%
  saveRDS(
    "temp-data/fully_cts_dist_balance.rds"
  )

disc_dist_balance = dist_balance_disc_fun(
    cluster.dist.to.pot, 
    interval_length = script_options$cts_interval)


construct_joint_test_m = function(object) {
  n_coef = length(coef(object))
  diag_m = diag(n_coef - 1)
  neg_1_m = matrix(-1, nrow = n_coef - 1, ncol = 1)

  hyp_m = cbind(neg_1_m, diag_m)
  return(hyp_m)
}


perform_cts_distance_boot = function(data, y_var, dist_var, interval_length, suppress_intercept = FALSE) {
    sym_dist_var = sym(dist_var)
    subset_data = data %>%
        select(all_of(y_var), all_of(dist_var), county, constituency) %>%
        na.omit() 
    if (!is.null(interval_length)) {
        subset_data = subset_data %>%
            mutate(
            dist_measure = cut_interval({{ sym_dist_var }}, length = interval_length )
            )
    } else {
        subset_data = subset_data %>%
            mutate(dist_measure = {{ sym_dist_var }})
    }
    if (suppress_intercept) {
        fml = as.formula(paste0(y_var, " ~ 0 + dist_measure"))
    } else {
        fml = as.formula(paste0(y_var, " ~  dist_measure"))
    }
    fit = lm(fml, data = subset_data)
    R = construct_joint_test_m(fit)
    
  coef(fit) %>% length()

    boot_p = mboottest(
        object = fit, 
        clustid = c("county", "constituency"),
        bootcluster = "constituency", 
        B = script_options$num_boot_draws, 
        R = R
    )$p_val
    return(lst(boot_p))
}


cts_boots = function(interval) {
  cts_indiv_boot = map(
    indiv_balance_vars, 
    ~perform_cts_distance_boot(
      data = analysis_school_data, 
      y_var = .x, 
      dist_var = "cluster.dist.to.pot", 
      interval_length = interval, 
      suppress_intercept = TRUE
      )
  )
  cts_know_boot = map(
    know_vars, 
    ~perform_cts_distance_boot(
      data = know_balance_data,
      y_var = .x, 
      dist_var = "cluster.dist.to.pot",
      interval_length = interval, 
      suppress_intercept = TRUE
    )
  )
  cts_baseline_boot = map(
    baseline_vars, 
    ~perform_cts_distance_boot(
      data = baseline_balance_data,
      y_var = .x, 
      dist_var = "cluster.dist.to.pot",
      interval_length = interval, 
      suppress_intercept = TRUE
    )
  )

  cts_boot_pvals = c(
      # indiv
      unlist(cts_indiv_boot),
      # know 
      unlist(cts_know_boot),
      # balance
      unlist(cts_baseline_boot)
  )

  names(cts_boot_pvals) = names(fully_cts_dist_balance)


  tidy_cts_boot_pvals_df = enframe(
      cts_boot_pvals
  ) %>%
  mutate(
    interval = interval
  )
  return(tidy_cts_boot_pvals_df)
}

tidy_cts_boot_pvals_df = cts_boots(script_options$cts_interval)

tidy_cts_boot_pvals_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            str_glue("cts_distance_binned_pvals_INTERVAL_{script_options$cts_interval}.csv")
        )
    )


tidy_cts_boot_many_pvals_df = map_dfr(
  c(100, 200, 300, 400, 500, 1000), 
  cts_boots
) 

tidy_cts_boot_many_pvals_df %>%
  write_csv(
        file.path(
            script_options$output_path,
            str_glue("cts_distance_binned_pvals_many_intervals.csv")
        )
  )



balance_data = lst(
  analysis_school_data,
  analysis_data, 
  endline_balance_data, 
  baseline_balance_data,
  endline_and_baseline_data,
  endline_vars, 
  baseline_vars,
  know_vars,
  indiv_balance_vars
)


saveRDS(
  balance_data, 
  file.path(
    script_options$output_path,
    "saved_balance_data.rds"
  )
)
