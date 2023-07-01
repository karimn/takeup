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
    --num-boot-draws=200 \
    --output-path=temp-data \
    --cts-interval=200
    " else commandArgs(trailingOnly = TRUE)
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
    fully_aware_externalities = case_when(
      neighbours_worms_affect == "yes" & worms_affect == "yes" ~ TRUE, 
      is.na(neighbours_worms_affect) | is.na(worms_affect) ~ NA,
      TRUE ~ FALSE
    ),
    partially_aware_externalities = case_when(
      neighbours_worms_affect == "yes" | worms_affect == "yes" ~ TRUE,
      is.na(neighbours_worms_affect) | is.na(worms_affect) ~ NA,
      TRUE ~ FALSE
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
  select(-digits_schooling, -schooling_years_plus)
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
    all_can_get_worms = map2_lgl(
      who_worms, 
      who_worms_other,
      ~any(
        "everyone" %in% .x | 
        str_detect(str_to_lower(.y), "any") | 
        (("adult" %in% .x) & ("child" %in% .x)) |
        ("adult" %in% .x | str_detect(str_to_lower(.y), "adult|man|woman|men|women|person")) & ("child" %in% .x | str_detect(str_to_lower(.y), "child|under|young|teenager|below"))
    )
    ), 
    correct_when_treat = map_lgl(when_treat, ~any(.x == "every 6 months")), 
    know_deworming_stops_worms = map2_lgl(
      stop_worms, 
      stop_worms_other,
      ~any(.x %in% c(
        "medicine", 
        "wearing shoes", 
        "using toilets", 
        "wash hands") | str_detect(.y, "cooked|prepar|cook")))
  ) 

baseline.data %>%
  select(who_worms) %>%
  mutate(
    ad_and_child = map_lgl(who_worms, ~any(("adult" %in% .x) & ("child" %in% .x)))
  ) %>%
  filter(ad_and_child == TRUE)

baseline.data %>%
  select(who_worms_other) %>%
  count(who_worms_other) %>%
  arrange(-n) %>%
  print(n = 100)

# creating a single treat x distance variable for balance testing
cluster_treat_df = analysis_data %>%
  mutate(
      treat_dist = paste0(
      "treat: ", 
      assigned.treatment,
      ", dist: ", dist.pot.group
      ) %>% factor()
  ) %>%
  select(cluster.id, treat_dist, dist_to_pot = cluster.dist.to.pot) %>%
  unique()





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



# Remove surveys that took less than 10 minutes or more than 3 hours
# (long right tail of survey times past 500 mins, not sure where we should cut here
# if at all)
baseline_balance_data = baseline_balance_data %>%
  mutate(time = lubridate::mdy_hms(endtime) - lubridate::mdy_hms(starttime)) %>%
  filter(time > 10,  time < 3*60) 



baseline_vars = c(
  "completed_primary", 
  "know_deworming_stops_worms",
  "treated_lgl", 
  "floor_tile_cement",
  "all_can_get_worms",
  "correct_when_treat",
  "fully_aware_externalities",
  "partially_aware_externalities"
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
) %>%
  rename(dist_to_pot = dist.to.pot)

  

# Probably a better way to get the schools in the sample
school_treat_df = analysis_school_data %>%
  filter(!is.na(assigned.treatment)) %>%
  select(treat_dist, cluster.dist.to.pot = standard_cluster.dist.to.pot,  cluster.id, county) %>%
  unique()



#### Endline
endline_balance_data = endline.data %>%
  inner_join(
    cluster_treat_df, 
    by = "cluster.id"
  ) %>%
  mutate(
    fully_aware_externalities = case_when(
      neighbours_worms_affect == "yes" & worms_affect == "yes" ~ TRUE, 
      is.na(neighbours_worms_affect) | is.na(worms_affect) ~ NA,
      TRUE ~ FALSE
    ),
    partially_aware_externalities = case_when(
      neighbours_worms_affect == "yes" | worms_affect == "yes" ~ TRUE,
      is.na(neighbours_worms_affect) | is.na(worms_affect) ~ NA,
      TRUE ~ FALSE
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
      all_can_get_worms,
      correct_when_treat, 
      know_deworming_stops_worms,
      constituency, cluster.id, 
      county,
      fully_aware_externalities,
      partially_aware_externalities
      ) %>%
    mutate(
      type = "endline"
    ),
  baseline_balance_data %>%
    select(
      treat_dist, 
      all_can_get_worms,
      correct_when_treat, 
      know_deworming_stops_worms,
      constituency, cluster.id,
      county,
      fully_aware_externalities,
      partially_aware_externalities
      ) %>%
    mutate(
      type = "baseline"
    )
) %>%
  na.omit()

endline_vars = c(
  "fully_aware_externalities",
  "partially_aware_externalities",
  "all_can_get_worms", 
  "correct_when_treat", 
  "know_deworming_stops_worms"
  )
## Fits
endline_balance_fit = feols(
    data = endline_balance_data, 
    .[endline_vars] ~ 0 + treat_dist + i(county, ref = "Busia"), 
    cluster = ~cluster.id
    ) 

baseline_balance_fit = feols(
    data = baseline_balance_data, 
    .[baseline_vars] ~ 0 + treat_dist + i(county, ref = "Busia"), 
    ~cluster.id
    ) 

indiv_balance_fit = feols(
    data = analysis_school_data, 
    .[indiv_balance_vars] ~ 0 + treat_dist + i(county, ref = "Busia"),
    cluster = ~cluster.id
    ) 

school_balance_fit = feols(
    data = analysis_school_data %>%
      select(any_of(balance_variables), treat_dist, county), 
    .[balance_variables] ~ 0 + treat_dist + i(county, ref = "Busia"),
    vcov = "HC"
  )

# put all the baseline balance fits into a list we can map over
balance_fits = c(
  indiv_balance_fit,
  list("lhs: cluster.dist.to.pot" = school_balance_fit),
  baseline_balance_fit
)


create_balance_comparisons = function(fit) {
  comp_df = avg_comparisons(
    fit,
    variables = list("treat_dist" = "all")
    ) %>%
    as_tibble()


  comp_df = comp_df %>%
    mutate(
      lhs_treatment = str_extract(contrast, "(?<=^treat: )\\w+"), 
      rhs_treatment = str_extract(contrast, "(?<=- treat: )\\w+"), 
      lhs_dist = str_extract(contrast, "(?<=dist: )\\w+"),
      rhs_dist = str_extract(contrast, "(?<=, dist: )\\w+$")
    ) 


  same_dist_subset_comp_df = comp_df  %>%
    filter(
      lhs_dist == rhs_dist,
      rhs_treatment == "control" | lhs_treatment == "control",
      lhs_treatment != rhs_treatment
    ) 

  same_dist_bra_cal_comp_df = comp_df %>%
    filter(lhs_dist == rhs_dist) %>%
    filter(str_detect(contrast, "bracelet") & str_detect(contrast, "calendar"))

    rhs_control_comp_df = same_dist_subset_comp_df %>%
      filter(rhs_treatment == "control") %>%
      bind_rows(
        same_dist_bra_cal_comp_df %>%
          filter(rhs_treatment == "calendar")
      )

    lhs_control_comp_df = same_dist_subset_comp_df %>%
      filter(rhs_treatment != "control") %>%
      bind_rows(
        same_dist_bra_cal_comp_df %>%
          filter(rhs_treatment == "bracelet")
      )

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
    ) %>%
    mutate(comp_type = "treatment")

    ## Now within treatment across distances

    dist_control_mean_df = fit %>%
      tidy(conf.int = TRUE) %>%
      filter(str_detect(term, "close")) %>%
      mutate(
        lhs_treatment = str_extract(
          term, 
          "(?<=treat: )\\w+"),
        rhs_treatment = NA,
        lhs_dist = "close",
        rhs_dist = NA
      )

    dist_comp_df = comp_df %>%
      filter(
        rhs_dist != lhs_dist,
        lhs_treatment == rhs_treatment
      )  %>%
      select(-contrast) %>%
      bind_rows(
        dist_control_mean_df
      ) %>%
      mutate(comp_type = "distance")

    final_clean_comp_df = bind_rows(
      rearranged_comp_df,
      dist_comp_df
    )
    return(final_clean_comp_df)

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
        lhs, term, estimate, std.error, p.value
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


perform_balance_joint_test = function(fit, var, joint_R, close_R, far_R) {
  county_0_mat = matrix(
    0,
    nrow = max(nrow(joint_R), nrow(close_R), nrow(far_R)),
    ncol = coef(fit) %>% length() - 8
    )

  resid_df = fixest::degrees_freedom(fit, type = "resid")
  close_test = car::lht(
    fit,
    cbind(close_R, county_0_mat[1:nrow(close_R), ]),
    error.df = resid_df,
    test = "F"
  )

  far_test = car::lht(
    fit,
    cbind(far_R, county_0_mat[1:nrow(far_R), ]),
    error.df = resid_df,
    test = "F"
  )

  joint_test = car::lht(
    fit,
    cbind(joint_R, county_0_mat[1:nrow(joint_R),]),
    error.df = resid_df,
    test = "F"
  )


  pvals = lst(
    joint_pval = joint_test$`Pr(>F)`[2],
    far_pval = far_test$`Pr(>F)`[2],
    close_pval = close_test$`Pr(>F)`[2]
  ) 

  return(pvals)
}

balance_joint_tests = map(
  balance_fits,
  ~perform_balance_joint_test(
    .x,
    joint_R = hyp_matrix,
    close_R = hyp_matrix_close,
    far_R = hyp_matrix_far
  )
)



balance_joint_tests %>%
    saveRDS(
        file.path(
            script_options$output_path,
            "cluster_balance_fits.rds"
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
      ~cluster.id
      ) 

perform_externality_cluster_boot = function(data, var, R) {
  subset_data = data %>%
    select(all_of(var), treat_dist, type, cluster.id, constituency) %>%
      na.omit()


  fml = as.formula(paste0(var, " ~ 0 + treat_dist:type"))
  fit = lm(
    formula = fml, 
    data = subset_data
  )

  if (is.matrix(R)) {
    joint_boot_output = mboottest(
      object = fit, 
      clustid = c("cluster.id"), 
      bootcluster = "cluster.id", 
      B = script_options$num_boot_draws, 
      R = R
    )
  } else {
    joint_boot_output = boottest(
      object = fit, 
      clustid = c("cluster.id"),
      param = names(coef(fit)), 
      bootcluster = "cluster.id", 
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
        lhs, term, estimate, std.error, p.value
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
balance_data = lst(
  analysis_school_data,
  analysis_data, 
  endline_balance_data, 
  baseline_balance_data,
  endline_and_baseline_data,
  endline_vars, 
  baseline_vars,
  indiv_balance_vars
)


saveRDS(
  balance_data, 
  file.path(
    script_options$output_path,
    "saved_balance_data.rds"
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
      cluster = ~cluster.id
      ) 

  ## Know

  ## baseline
  baseline_balance_cts_fit = feols(
      data = baseline_balance_data %>%
        mutate(dist_measure = {{ rhs_var }}/1000), 
      .[baseline_vars] ~  dist_measure + dist_measure^2, 
      ~cluster.id
      ) 

    return(
      c(
        indiv_balance_cts_fit,
        baseline_balance_cts_fit
      )
    )
}

dist_balance_disc_fun = function(rhs_var, interval_length) {
  analysis_school_data = analysis_school_data %>%
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
      cluster = ~cluster.id
      ) 

  ## baseline
  baseline_balance_disc_fit = feols(
      data = baseline_balance_data,
      .[baseline_vars] ~ 0 + dist_measure, 
      ~cluster.id
      ) 

    return(
      c(
        indiv_balance_disc_fit,
        baseline_balance_disc_fit
      )
    )
}


# fully_cts_dist_balance = dist_balance_cts_fun(cluster.dist.to.pot)
# disc_dist_balance = dist_balance_disc_fun(cluster.dist.to.pot, interval_length = 500)

# fully_cts_dist_balance = dist_balance_cts_fun(cluster.dist.to.pot)


# fully_cts_dist_balance %>%
#   saveRDS(
#     "temp-data/fully_cts_dist_balance.rds"
#   )

# disc_dist_balance = dist_balance_disc_fun(
#     cluster.dist.to.pot, 
#     interval_length = script_options$cts_interval)


# construct_joint_test_m = function(object) {
#   n_coef = length(coef(object))
#   diag_m = diag(n_coef - 1)
#   neg_1_m = matrix(-1, nrow = n_coef - 1, ncol = 1)

#   hyp_m = cbind(neg_1_m, diag_m)
#   return(hyp_m)
# }


# perform_cts_distance_boot = function(data, y_var, dist_var, interval_length, suppress_intercept = FALSE) {
#     sym_dist_var = sym(dist_var)
#     subset_data = data %>%
#         select(all_of(y_var), all_of(dist_var), cluster.id, constituency) %>%
#         na.omit() 
#     if (!is.null(interval_length)) {
#         subset_data = subset_data %>%
#             mutate(
#             dist_measure = cut_interval({{ sym_dist_var }}, length = interval_length )
#             )
#     } else {
#         subset_data = subset_data %>%
#             mutate(dist_measure = {{ sym_dist_var }})
#     }
#     if (suppress_intercept) {
#         fml = as.formula(paste0(y_var, " ~ 0 + dist_measure"))
#     } else {
#         fml = as.formula(paste0(y_var, " ~  dist_measure"))
#     }
#     fit = lm(fml, data = subset_data)
#     R = construct_joint_test_m(fit)
    
#   coef(fit) %>% length()

#     boot_p = mboottest(
#         object = fit, 
#         clustid = c("cluster.id"),
#         bootcluster = "cluster.id", 
#         B = script_options$num_boot_draws, 
#         R = R
#     )$p_val
#     return(lst(boot_p))
# }


# cts_boots = function(interval) {
#   cts_indiv_boot = map(
#     indiv_balance_vars, 
#     ~perform_cts_distance_boot(
#       data = analysis_school_data, 
#       y_var = .x, 
#       dist_var = "cluster.dist.to.pot", 
#       interval_length = interval, 
#       suppress_intercept = TRUE
#       )
#   )
#   cts_know_boot = map(
#     know_vars, 
#     ~perform_cts_distance_boot(
#       data = know_balance_data,
#       y_var = .x, 
#       dist_var = "cluster.dist.to.pot",
#       interval_length = interval, 
#       suppress_intercept = TRUE
#     )
#   )
#   cts_baseline_boot = map(
#     baseline_vars, 
#     ~perform_cts_distance_boot(
#       data = baseline_balance_data,
#       y_var = .x, 
#       dist_var = "cluster.dist.to.pot",
#       interval_length = interval, 
#       suppress_intercept = TRUE
#     )
#   )

#   cts_boot_pvals = c(
#       # indiv
#       unlist(cts_indiv_boot),
#       # know 
#       unlist(cts_know_boot),
#       # balance
#       unlist(cts_baseline_boot)
#   )

#   names(cts_boot_pvals) = names(fully_cts_dist_balance)


#   tidy_cts_boot_pvals_df = enframe(
#       cts_boot_pvals
#   ) %>%
#   mutate(
#     interval = interval
#   )
#   return(tidy_cts_boot_pvals_df)
# }

# tidy_cts_boot_pvals_df = cts_boots(script_options$cts_interval)

# tidy_cts_boot_pvals_df %>%
#     write_csv(
#         file.path(
#             script_options$output_path,
#             str_glue("cts_distance_binned_pvals_INTERVAL_{script_options$cts_interval}.csv")
#         )
#     )


# tidy_cts_boot_many_pvals_df = map_dfr(
#   c(100, 200, 300, 400, 500, 1000), 
#   cts_boots
# ) 

# tidy_cts_boot_many_pvals_df %>%
#   write_csv(
#         file.path(
#             script_options$output_path,
#             str_glue("cts_distance_binned_pvals_many_intervals.csv")
#         )
#   )
