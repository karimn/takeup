# ---- name.match.monitored
# Name matching function. Given the census data and take-up from a particular attempt to 
# find individuals whose names were recorded at the PoT.
name.match.monitored <- function(census.data, 
                                 takeup.data, 
                                 max.cost = 1,
                                 suffix = NULL) { # This is the maximum number of "edits" or difference allowed between names
  dist.mat <- adist(census.data$name1st, takeup.data$name1st, ignore.case = TRUE) +
    adist(census.data$last_name, takeup.data$last_name, ignore.case = TRUE)
  
  census.data %>% 
    mutate(min.name.match.dist = aaply(dist.mat, 1, . %>% min(na.rm = TRUE)) %>% na_if(Inf),
           which.min.name.match.dist = ifelse(!is.na(min.name.match.dist), 
                                              aaply(dist.mat, 1, . %>% 
                                                      which.min %>% 
                                                      magrittr::extract(takeup.data$KEY.survey.individ, .)),
                                              NA),
           dewormed.matched = !is.na(which.min.name.match.dist) & min.name.match.dist <= max.cost,
           which.min.name.match.dist = ifelse(dewormed.matched, which.min.name.match.dist, NA)) %>% 
    # rename_(.dots = c(str_subset(names(.), "min\\.name\\.match\\.dist$"), "dewormed.matched") %>% 
    #           setNames(paste0(., suffix))) %>% 
    select(KEY.individ, starts_with("dewormed.matched"), contains("min.name.match.dist"))
}
# ---- end
# Misc Functions and Constants ----

reg.covar <- c("school", "floor", "ethnicity", "sms.ctrl.subpop", "age", "age_squared", "gender")
census.reg.covar <- c("gender", "age.census", "age.census_squared")

prepare.consent.dewormed.data <- function(.all.endline.data, .reconsent.data) {
  list(endline.survey = .all.endline.data, 
       reconsent = .reconsent.data %>% mutate_at(vars(monitor.consent), as.logical)) %>% 
    map_df(. %>% select(KEY.individ, monitor.consent, dewormed.reported), .id = "data.source") %>% 
    filter(!is.na(monitor.consent), !is.na(dewormed.reported)) %>%  
    group_by(KEY.individ) %>%  
    summarize(monitor.consent = any(!is.na(monitor.consent) & monitor.consent), # Consider as reconsented if at least one acceptance
              dewormed.reported = ifelse(n_distinct(dewormed.reported) == 1, first(dewormed.reported), NA)) %>% # Multiple contradictory responses
    ungroup
}

prepare.analysis.data <- function(.census.data, .takeup.data, .endline.data, .consent.dewormed.reports, .cluster.strat.data,
                                  max.name.match.cost = 1) {
  dewormed.day.data <- .takeup.data %>% 
    filter(!is.na(KEY.individ)) %>% 
    group_by(KEY.individ) %>% 
    summarize(dewormed.day = min(deworming.day)) %>% # If dewormed multiple times, take the first day only
    ungroup
  
  analysis.data <- .census.data %>% 
           filter(!is.na(wave)) %>%  # Remove clusters no longer in study
    left_join(.consent.dewormed.reports, "KEY.individ") %>% 
    mutate(dewormed = KEY.individ %in% .takeup.data$KEY.individ, # TRUE if individual found in take-up data
           dewormed = ifelse(monitored, dewormed, NA),
           monitor.consent = !is.na(monitor.consent) & monitor.consent,
           sms.treated = sms.treatment %in% c("social.info", "reminder.only"),
           sms.treatment.2 = fct_explicit_na(sms.treatment, "sms.control") %>% fct_relevel("sms.control"),
           have.phone.bool = have_phone == "Yes",
           sms.ctrl.subpop = if_else(is.na(sms.ctrl.subpop) & !sms.treated, 
                                     if_else(have.phone.bool, "phone.owner", "non.phone.owner"),
                                     sms.ctrl.subpop)) %>% # NA if not in the monitored group
    left_join(dewormed.day.data, "KEY.individ") %>% 
    filter(!sms.treated | (!is.na(sms.consent) & sms.consent)) # Drop those assigned to SMS treatment but were not consented (hence not treated)
  
  analysis.data %>% 
    filter(is.na(dewormed) | !dewormed) %>% # For anyone in study with with unknown or negative deworming status
    group_by(cluster.id) %>% 
    do(name.match.monitored(., filter(takeup.data, cluster.id %in% .$cluster.id), max.cost = max.name.match.cost)) %>% 
    ungroup %>% 
    # select(KEY.individ, starts_with("dewormed.matched"), contains("min.name.match.dist")) %>% 
    right_join(analysis.data, c("cluster.id", "KEY.individ")) %>% 
    left_join(transmute(takeup.data, KEY.survey.individ, dewormed.day.matched = deworming.day), c("which.min.name.match.dist" = "KEY.survey.individ")) %>% 
    mutate(monitored = !is.na(monitored) & !is.na(wave) & monitored, # Remove those dropped from the study 
           dewormed.any = (!is.na(dewormed) & dewormed) | dewormed.matched,
           # dewormed.any_0 = (!is.na(dewormed) & dewormed) | dewormed.matched_0,
           dewormed.day.any = if_else(!is.na(dewormed.day), as.integer(dewormed.day), dewormed.day.matched), 
           gender = factor(gender, levels = 1:2, labels = c("male", "female"))) %>% 
    left_join(select(.endline.data, KEY.individ, age, school, floor, ethnicity, ethnicity2, any.sms.reported, gift_choice,
                     hh_cal, cal_value, hh_bracelet, number_bracelet), "KEY.individ") %>% 
    mutate_at(vars(age, age.census), funs(squared = (.)^2)) %>% 
    left_join(select(.cluster.strat.data, wave, county, cluster.id, dist.pot.group), c("wave", "county", "cluster.id")) %>% 
    `attr<-`("class", c("takeup_df", class(.))) %>%
    unite(county_dist_stratum, county, dist.pot.group, remove = FALSE) %>% 
    unite(county_dist_mon_stratum, county, dist.pot.group, true.monitored, remove = FALSE) %>% 
    mutate_at(vars(county, county_dist_stratum, county_dist_mon_stratum), factor) %>% 
    mutate(mon_status = factor(true.monitored, levels = c(TRUE, FALSE), labels = c("monitored", "unmonitored")),
           name_matched = !true.monitored,
           phone_owner = sms.treatment.2 != "sms.control" | sms.ctrl.subpop == "phone.owner")
}

multi.factor <- function(.col, labels, levels, ...) {
  if (missing(levels)) {
    levels <- c(seq_along(labels), 97:99)
    labels %<>% c("prefer not say", "DK", "other")
  } 
  
  map(str_split(.col, " "), factor, levels = levels, labels = labels, ...)
}

yes.no.factor <- function(.col, .yes.no = 1:2) .col %>% 
  factor(levels = c(.yes.no, 97:98), labels = c("yes", "no", "prefer not say", "DK"))

base.prepare.baseline.endline.data <- function(.data, .cluster.strat.data) { #, .census.data) {
  .data %>% 
    mutate(who_worms = multi.factor(who_worms, 
                                    labels = c("child", "adult", "sick", "healthy", "pregnant", "old", "everyone")), 
           effect_worms = multi.factor(effect_worms, 
                                       labels = c("stomachache", "malnourishment", "stop child growing", "tired", "diarrhea", 
                                                  "stop child school")),
           how_spread = multi.factor(how_spread, 
                                     labels = c("drinking dirty water", "open def", "swim/bath dirty water", "no shoes", "no washing hands",
                                                "sex", "uncooked food")),
           stop_worms = multi.factor(stop_worms, 
                                     labels = c("wash hands", "wearing shoes", "using toilets", "drink clean water", "medicine", "clean home")),
           when_treat = multi.factor(when_treat, 
                                     labels = c("every week", "every month", "every 2 months", "every 3 months", "every 6 months", 
                                                "every year", "never", "when symptoms", "hw says")),
           have_phone = factor(have_phone, levels = 0:1, labels = c("No", "Yes")),
           school = factor(school,
                           levels = c(1:15, 99, 97),
                           labels = c("Never gone to school", paste("Primary", 1:8), paste("Secondary", 1:4), "College", "University", "Other", "Prefer Not to Say")),
           floor = factor(floor, levels = c(1:3, 99, 97), labels = c("Cement", "Earth", "Tiles", "Other", "Prefer Not to Say")),
           ethnicity = factor(ethnicity, levels = c(1:7, 99, 97),
                              labels = c("Luo", "Luhya", "Kisii", "Kalengin", "Kikukyu", "Kamba", "Iteso", "Other", "Prefer Not to Say"))) %>% 
    mutate_at(vars(worms_affect, neighbours_worms_affect), funs(yes.no.factor(., .yes.no = 1:0))) %>% 
    mutate_at(vars(spread_worms), yes.no.factor) %>% 
    mutate_at(vars(school, floor, ethnicity), funs(recode_factor(., "Other" = NA_character_, "Prefer Not to Say" = NA_character_) %>% 
                                                     fct_drop)) %>% 
    mutate(ethnicity2 = fct_lump(ethnicity, other_level = "Other Ethnicities", prop = 0.05)) %>% 
    left_join(select(.cluster.strat.data, wave, county, cluster.id, assigned.treatment, dist.pot.group), c("wave", "county", "cluster.id")) %>% 
    tidyr::unite(county_dist_stratum, county, dist.pot.group, remove = FALSE) %>% 
    mutate_at(vars(county, county_dist_stratum), factor)
}

prepare.endline.data <- function(.data, .census.data, .cluster.strat.data) {
  .data %>% 
    filter(present, interview, consent) %>% 
    arrange(KEY.individ, SubmissionDate) %>% 
    group_by(KEY.individ) %>% # If more than one entry for an individual, take first one (there are 22 such individuals)
    filter(row_number() == 1) %>% 
    ungroup %>% 
    base.prepare.baseline.endline.data(.cluster.strat.data) %>% 
    mutate_at(vars(know_deworm, chv_visit, flyer, any.sms.reported, hh_bracelet, hh_cal, cal_value), 
              funs(yes.no.factor(., .yes.no = 1:0))) %>% 
    mutate_at(vars(treat_begin, days_available, treat_end), funs(factor(., levels = c(1, 98), c("knows", "DK")))) %>% 
    mutate(treat_begin_date = ymd(sprintf("2016-%d-%d", month_treat_begin, day_treat_begin)),
           treat_end_date = ymd(sprintf("2016-%d-%d", month_treat_end, day_treat_end)),
           #where_offered = labelled(where_offered, c("somewhere else" = 0, "home" = 3, "DK" = 98)) %>% as_factor)
           find_out = multi.factor(find_out, 
                                   levels = c(1:9, 99), 
                                   labels = c("friend", "family", "chv", "elder", "church", "flyer", "poster", "enumerator", "baraza",
                                              "other")),
           gift_choice = factor(gift_choice, levels = 1:3, labels = c("bracelet", "calendar", "neither"))) %>%
    left_join(select(.census.data, KEY.individ, dist.to.pot, sms.ctrl.subpop, true.monitored, monitored), "KEY.individ")  
           # text_content = factor(text_content, levels = c(1:3, 99), labels = c("reminders", "when/where", "social info", "other")))
}

prepare.baseline.data <- function(.data, .cluster.strat.data) {
  .data %>% 
    base.prepare.baseline.endline.data(.cluster.strat.data) %>% 
    mutate(more_less = factor(more_less, levels = 1:4, labels = c("more", "less", "no diff", "DK")),
           treated_when = factor(treated_when, 
                                 levels = c(1:9, 97), 
                                 labels = c("1-2 mon", "3-5 mon", "6-7 mon", "8-9 mon", "10-11 mon",
                                            "1 year", "2 year", "3 year", "4 year or more", 
                                            "prefer no say")),
           who_treated = factor(who_treated, 
                                levels = 1:3, 
                                labels = c("child", "adult", "both")),
           dworm_proportion = factor(dworm_proportion, 
                                      levels = c(1:6, 97:98),
                                      labels = c("few", "nearly half", "half", "more than half", "many", "all", 
                                                 "prefer not say", "DK")),
           ink_more_less = factor(ink_more_less, levels = c(1:3, 97:98), labels = c("more", "less", "same", 
                                                                                    "prefer not say", "DK"))) %>% 
    mutate_at(vars(treated, family_treated), yes.no.factor) %>% 
    mutate_at(vars(few_deworm, many_deworm, matches("(praise|stigma)_(immunize|clothes|deworm|defecate)[A-D]$")), 
              funs(factor(., levels = 0:2, labels = c("no", "yes", "maybe")))) %>% 
    mutate_at(vars(treated_where, where_family_treated), 
              funs(factor(., 
                          levels = c(1:4, 97:99), 
                          labels = c("bought", "school MDA", "non-school MDA", "hosp/clinic", 
                                     "prefer not say", "DK", "other"))))
}


prepare.cluster.takeup.data <- function(.data, 
                                        monitored.only = TRUE, 
                                        consented.only = monitored.only, 
                                        exclude.baseline.sample = TRUE,
                                        add_group_by = NULL) {
 .data %>% 
   filter(monitored | !monitored.only, 
          monitor.consent | !consented.only, 
          !(hh.baseline.sample & exclude.baseline.sample)) %>% 
   transmute(county, dist.pot.group, county_dist_stratum, cluster.id, assigned.treatment, sms.treatment = sms.treatment.2, dewormed.any, mon_status, phone_owner) %>% 
   # unite(stratum, county, dist.pot.group, sep = " ") %>% 
   group_by_(.dots = c("assigned.treatment", "sms.treatment", "county_dist_stratum", "cluster.id", add_group_by)) %>% 
   # group_by(assigned.treatment, sms.treatment, stratum, cluster.id) %>% 
   summarize(takeup.prop = mean(dewormed.any)) %>% 
   # group_by(assigned.treatment, sms.treatment, stratum) %>% 
   group_by_(.dots = c("assigned.treatment", "sms.treatment", "county_dist_stratum", add_group_by)) %>% 
   mutate(outlier = is_outlier(takeup.prop)) %>% 
   ungroup
}

bstrp.uniquify <- function(.data) {
  .data %>% 
    group_by(KEY.individ) %>% 
    mutate(KEY.individ.unique = paste(KEY.individ, seq(1, n()), sep = "-")) %>% 
    ungroup %>% 
    mutate(KEY.individ = KEY.individ.unique)
}

nest_exclude <- function(data, key_col, exclude_nest_cols = character()) {
  dplyr::select_vars(colnames(data), everything()) %>% 
    unname %>% 
    dplyr::setdiff(exclude_nest_cols) %>% 
    nest_(data, key_col = key_col, nest_cols = .)
}

calc.strata.stats <- function(.data, 
                              .treatment = c("assigned.treatment"), 
                              .strat.by = c("county", "dist.pot.group"), # This is what we stratified the experiment by
                              .interact.with = NULL) {
  calc.stratum.ate <- function(.stratum.data) {
    make.cross.assign.tibble <- function(.col) {
      matrix(.col, 
             nrow = length(.col), 
             ncol = length(.col), 
             dimnames = list(.stratum.data$treatment.group)) %>% 
        t %>% as_tibble %>% 
        mutate(lhs.treatment.group = .stratum.data$treatment.group) %>% 
        gather(rhs.treatment.group, val, -lhs.treatment.group) 
    }
   
    .stratum.data %>% {
      generated.cross.tbl <- select(., 
                                    stratum.assign.mean.dewormed, 
                                    stratum.assign.size, 
                                    stratum.assign.sample.var.dewormed) %>% 
        map2(paste0("stratum.", c("ate", "size", "sample.var")), 
             ~ make.cross.assign.tibble(.x) %>% 
               set_names(str_replace(names(.), fixed("val"), fixed(.y)))) %>% 
        reduce(function(.left, .right) left_join(.left, .right, c("lhs.treatment.group", "rhs.treatment.group")))
        
      left_join(., generated.cross.tbl, c("treatment.group" = "lhs.treatment.group")) 
    } %>% 
      filter(treatment.group != rhs.treatment.group) %>% 
      mutate(stratum.ate = stratum.assign.mean.dewormed - stratum.ate,
             stratum.size = stratum.size + stratum.assign.size,
             stratum.sample.var = stratum.sample.var + stratum.assign.sample.var.dewormed)
  }
 
  calc.strata.ate <- function(.data) .data %>% 
    unite_("treatment.group", .treatment) %>% 
    group_by_(.dots = c(.strat.by, "treatment.group")) %>% 
    summarize(stratum.assign.size = n(),
              stratum.assign.mean.dewormed = mean(dewormed.any),
              stratum.assign.sample.var.dewormed = var(dewormed.any)/stratum.assign.size) %>%  
    ungroup %>% 
    mutate(total.size = sum(stratum.assign.size)) %>% 
    group_by_(.dots = .strat.by) %>% 
    do(calc.stratum.ate(.)) %>% 
    group_by_(.dots = c("treatment.group", .interact.with)) %>% 
    do(mutate(., treatment.mean.dewormed = weighted.mean(.$stratum.assign.mean.dewormed, .$stratum.assign.size))) %>% 
    ungroup
  
  strata.data <- .data %>% 
    calc.strata.ate
  
  strata.data %>% {
    left_join(distinct_(., .dots = c(.strat.by, "treatment.group", .interact.with), .keep_all = TRUE) %>% 
                select(-c(rhs.treatment.group, stratum.ate, stratum.size, stratum.sample.var)) %>%
                nest_exclude("strata.data", c("treatment.group", .interact.with)),
              group_by_(., .dots = c("treatment.group", "rhs.treatment.group", .interact.with)) %>%
                summarize(ate = weighted.mean(stratum.ate, stratum.size),
                          sample.var = sum(stratum.sample.var * ((stratum.size/total.size)^2)),
                          treatment.mean.dewormed = first(treatment.mean.dewormed)) %>% 
                ungroup %>% 
                nest_exclude("treatment.data", c("treatment.group", .interact.with)),
              c("treatment.group", .interact.with))
  }
}


analyze.neyman.blk.bs <- function(.data, .reps = 1000, .interact.with = NULL, ...) {
  .data %<>% 
    select(county, dist.pot.group, cluster.id, assigned.treatment, sms.treatment, dewormed.any) %>% {
     original.stats <- calc.strata.stats(., .interact.with = .interact.with, ...)
     
     if (.reps > 0) {
       bs.analysis <- block.bootstrap(., .reps, "cluster.id") %>%
         do(calc.strata.stats(., .interact.with = .interact.with, ...)) %>% 
         ungroup %>% 
         select_(.dots = c("replicate", "treatment.group", .interact.with, "treatment.data")) %>% 
         unnest(treatment.data) %>% 
         group_by_(.dots = c("treatment.group", "rhs.treatment.group", .interact.with)) %>%
         summarize_at(vars(ate, treatment.mean.dewormed), funs(mean(., na.rm = TRUE), sd(., na.rm = TRUE))) %>%
         mutate(ci.lower.treatment.mean.dewormed = treatment.mean.dewormed_mean - 1.64 * treatment.mean.dewormed_sd,
                ci.upper.treatment.mean.dewormed = treatment.mean.dewormed_mean + 1.64 * treatment.mean.dewormed_sd,
                ci.lower.ate = ate_mean - 1.64 * ate_sd,
                ci.upper.ate = ate_mean + 1.64 * ate_sd,
                treatment.mean.dewormed_tstat = treatment.mean.dewormed_mean/treatment.mean.dewormed_sd,
                ate_tstat = ate_mean/ate_sd,
                ate_pvalue = calc.pvalue(ate_tstat))
       
       original.stats %<>% 
         unnest(treatment.data) %>% 
         left_join(bs.analysis, c("treatment.group", "rhs.treatment.group", .interact.with)) %>% 
         nest_exclude("treatment.data", c("treatment.group", .interact.with)) %>% 
         left_join(original.stats %>% select(-treatment.data), ., c("treatment.group", .interact.with))
     } 
     
     return(original.stats)
   }
}

get_treatment_map <- function(analysis_data, analysis_formula) {
  treatment_vars <- analysis_formula %>% 
    terms() %>% 
    delete.response() %>% 
    all.vars() 
  
  analysis_data %>% 
    expand_(treatment_vars) %>% 
    semi_join(analysis_data, treatment_vars) 
}

get_treatment_map_design_matrix <- function(analysis_data, analysis_formula, treatment_map = get_treatment_map(analysis_data, analysis_formula)) {
  rhs_formula <- analysis_formula %>% 
    terms() %>% 
    delete.response() 
  
  old.contrasts <- getOption("contrasts")
  old.contrasts["unordered"] <- "contr.Treatment"
  old.options <- options(contrasts = old.contrasts)
  on.exit(options(old.options), add = TRUE)
 
  treatment_map %>%  
    modelr::model_matrix(rhs_formula) %>% 
    # magrittr::extract(, map_lgl(., ~ n_distinct(.) > 1)) %>% 
    set_colnames(str_replace_all(colnames(.), "([^:\\[]+)\\[T\\.([^\\]]+)]", "\\2"))
}

# Plotting code -----------------------------------------------------------

know.bel.cat.plot <- function(var, .baseline.data = baseline.data, .endline.data = endline.data, na.rm = FALSE) {
  list(Baseline = .baseline.data, Endline = .endline.data) %>% 
    compact %>% 
    map_df(select_, .dots = c(var, "KEY"), .id = "survey.type") %>% 
    unnest_(var) %>% { 
      if(na.rm) filter_(., sprintf("!is.na(%s)", var)) else return(.)
    } %>% 
    group_by(survey.type) %>% 
    do({
      mutate(count_(., c("survey.type", var)), n = n/n_distinct(.$KEY))
    }) %>% 
    ungroup %>% {
      if (n_distinct(.$survey.type) > 1) {
        ggplot(., aes_string(x = var, "n", fill = "survey.type")) +
          scale_fill_discrete("") 
      } else {
        ggplot(., aes_string(x = var, "n")) 
      }
    } %>% {
      . +
        geom_col(position = "dodge", alpha = 0.5) +
        coord_flip() +
        labs(x = "", y = "Proportion")
    }
}

multi.know.bel.cat.plot <- function(question.info, 
                                    .baseline.data = baseline.data, 
                                    .endline.data = endline.data, 
                                    .census.data = NULL,
                                    na.rm = FALSE, 
                                    preprocess.fun = function(.data) .data) {
  stopifnot(is.data.frame(question.info))
  
  list(Census = .census.data, 
       Baseline = .baseline.data %>% { if (!empty(.)) mutate(., KEY.individ = KEY) else return(.) }, 
       Endline = .endline.data) %>% 
    compact %>% 
    map_df(~ select_(.x, .dots = c(intersect(names(.x), question.info$col.name), "KEY.individ")), .id = "survey.type") %>% 
    map_df(question.info$col.name, 
           function(.col.name, .data) {
             .data %>% 
               unnest_(.col.name) %>% {
                 if(na.rm) filter_(., sprintf("!is.na(%s)", .col.name)) else return(.)
               } %>%
               group_by(survey.type) %>% 
               do(mutate(count_(., c("survey.type", .col.name)), n = n/n_distinct(.$KEY.individ))) %>% 
               ungroup %>% 
               rename_(response = .col.name) %>% 
               mutate(question = .col.name)
           },
         .data = .) %>% 
    mutate(question = factor(question, levels = question.info$col.name, labels = question.info$label),
           response = fct_recode(response, "Don't Know" = "DK") %>% 
             fct_relabel(str_to_title) %>% 
             fct_relabel(function(.label) str_replace(.label, "Chv", "CHV")) %>% 
             fct_reorder(n)) %>%
    preprocess.fun %>%
    mutate(survey.type = fct_relevel(survey.type, "Census", "Baseline", "Endline")) %>% { 
      if (n_distinct(.$survey.type) > 1) {
        ggplot(., aes(x = response, n, fill = survey.type)) +
          scale_fill_discrete("") 
      } else {
        ggplot(., aes(x = response, n)) 
      }
    } %>% {
      . +
        geom_col(position = "dodge", alpha = 0.5, na.rm = na.rm) +
        coord_flip() +
        labs(x = "", y = "Proportion") 
    } %>% {
      . +
        facet_grid(question ~ ., scales = "free", space = "free", labeller = label_wrap_gen()) +
        theme(strip.text.y = element_text(angle = 0)) 
    } 
}

plot.pref.unfaceted <- function(.data, .second.group.by = NULL) {
  .data %>% 
    filter(!is.na(gift_choice), monitored, monitor.consent, !hh.baseline.sample, !is.na(sms.treatment)) %>% 
    group_by_(.dots = c("assigned.treatment", .second.group.by)) %>% 
    mutate(arm.size = n()) %>% 
    group_by(gift_choice, add = TRUE) %>%
    summarize(pref.prop = n()/first(arm.size), arm.size = first(arm.size)) %>% 
    ungroup %>% 
    mutate(assigned.treatment = paste(assigned.treatment, "arm")) %>% {
      bind_rows(., 
                group_by_(., .dots = c("gift_choice", .second.group.by)) %>%
                  summarize(assigned.treatment = "(All Arms)", 
                            pref.prop = (pref.prop %*% arm.size)/sum(arm.size)) %>% 
                  ungroup)
    } %>%
    mutate_at(vars(gift_choice, assigned.treatment), funs(fct_relabel(factor(.), str_to_title))) %>% 
    mutate(assigned.treatment = fct_relevel(assigned.treatment, 
                                            "(All Arms)", "Control Arm", "Ink Arm", "Calendar Arm", "Bracelet Arm")) %>% 
    ggplot(aes(gift_choice, pref.prop)) +
    geom_col(alpha = 0.5, color = alpha("black", 0.5)) +
    labs(x = "Preferred Gift", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


prep.sms.ctrl.plot.data <- function(.reg.output, .interact.with = NULL) {
  incentive.treatment.terms <- c("control", "ink", "calendar", "bracelet")
  
  factor.incentive.treatment <- . %>% 
    factor(., levels = incentive.treatment.terms, labels = str_to_title(incentive.treatment.terms))
  
  calendar.ref.terms <- lst(ref.grp = c("(intercept) + calendar", "bracelet - calendar"))
  
  update.formula.with.interaction <- function(.formula, .keep.intercept = FALSE, .add.interact.intercept = TRUE) {
    append.intercept <- .keep.intercept & str_detect(.formula, fixed("(intercept)"))
    
    .formula %>% 
      str_replace(., pattern = "\\(intercept\\)\\s*(\\+\\s*)?", "") %>% 
      str_split("\\s+", simplify = FALSE) %>% 
      map_if(~ length(.x) > 1 | nchar(.x[1]) > 0, ~ c(.x, "+", str_replace(.x, "([^-+]+)", sprintf("\\1:%s", .interact.with)))) %>% 
      map(paste, collapse = " ") %>% 
      map_if(~ .add.interact.intercept & !str_detect(.x, "-") & nchar(.x) > 0, ~ paste(.x, .interact.with, sep = " + ")) %>% 
      map_if(~ .add.interact.intercept & nchar(.x) == 0, ~ .interact.with) %>% 
      map_if(append.intercept, ~ paste("(intercept) +", .x)) %>% 
      unlist
  }
  
  if (!is.null(.interact.with)) {
    calendar.ref.terms %<>%
      update_list(compare.grp = update.formula.with.interaction(.$ref.grp, .keep.intercept = TRUE)) 
  }
  
  
  calendar.ref.plot.data <- map_df(calendar.ref.terms, . %>% 
                                     linear_tester(.reg.output, .) %>% 
                                     magrittr::extract(., c(rep(1, 2), 2), ) %>% 
                                     mutate(incentive.treatment = c("calendar", rep("bracelet", 2)),
                                            ref = str_detect(linear.test, "(intercept)"),
                                            bar.size = estimate + if_else(ref, 0, estimate[1])),
                                   .id = "grp")
 
  add.interactions <- function(.data) {
    if (is.null(.interact.with)) {
      return(.data)
    } else {
      .data$term %>%
        update.formula.with.interaction(.add.interact.intercept = FALSE) %>% 
        map_if(~ nchar(.x) == 0, ~ sprintf("(intercept) + %s", .interact.with)) %>% 
        unlist %>% 
        linear_tester(.reg.output, .) %>% 
        bind_cols(.data %>% select(incentive.treatment)) %>% 
        mutate(ref = row_number() <= 4,
               bar.size = estimate + if_else(ref, 0, estimate[1]),
               estimate = if_else(estimate < 0, 0, estimate)) %>% 
        bind_rows(ref.grp = .data, compare.grp = ., .id = "grp")
    }
  } 
  
  .reg.output %>% 
    tidy %>% {
      if (is.null(.interact.with)) return(.) else filter(., !str_detect(.$term, .interact.with))
    } %>%
    filter(term %in% c("(intercept)", incentive.treatment.terms)) %>% 
    bind_rows(magrittr::extract(., rep(1, 3), ), .) %>% 
    mutate(incentive.treatment = incentive.treatment.terms %>% { c("control", rep(.[-1], 2)) },
           ref = term == "(intercept)",  
           bar.size = estimate + if_else(ref, 0, estimate[1]),
           estimate = if_else(estimate < 0, 0, estimate)) %>% 
    add.interactions %>% 
    bind_rows(control = ., calendar = calendar.ref.plot.data, .id = "ref.treatment") %>% 
    mutate(ci.lb = bar.size - std.error * 1.64,
           ci.ub = bar.size + std.error * 1.64) %>% 
    mutate_at(vars(incentive.treatment, ref.treatment), factor.incentive.treatment)  
}

plot.sms.ctrl.takeup <- function(.preped.data, .facet.formula = . ~ ref.treatment) {
  ggplot(.preped.data, aes(x = incentive.treatment)) +
    geom_col(aes(y = estimate, fill = !ref), 
             position = position_stack(reverse = TRUE), alpha = 0.25, color = "grey50") +
    geom_text_repel(aes(y = bar.size, label = sprintf("%.2f", bar.size)), 
                    nudge_x = 0.2, nudge_y = 0.02, segment.colour = "grey50",
                    data = . %>% filter(!ref | incentive.treatment == ref.treatment)) +
    geom_point(aes(y = bar.size), data = . %>% filter(!ref | incentive.treatment == ref.treatment)) +
    # geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
    #               width = 0.075, position = "dodge", 
    #               data = . %>% filter(!ref)) +
    geom_crossbar(aes(y = bar.size, ymin = ci.lb, ymax = ci.ub), 
                  width = 0.075, fatten = 0, 
                  data = . %>% filter(!ref)) +
    scale_fill_discrete("", labels = c("Reference Group", "Comparison Group")) +
    scale_y_continuous("Proportion of Take-up", breaks = seq(0, 1, 0.05)) +
    scale_x_discrete("Treatment") +
    facet_grid(.facet.formula, scales = "free_x", space = "free_x", 
               labeller = labeller(.cols = as_labeller(. %>% paste0("Reference Group: ", .)),
                                   .rows = as_labeller(. %>% paste0("Distance: ", .)))) +
    labs(title = "Estimated Take-up in Response to Incentive Treatment",
         caption = "Intervals shown identify the 90% confidence intervals, estimated using cluster robust standard errors.\nConfidence intervals test the null hypothesis of no difference in take-up from the reference group.") 
}

prep.sms.treat.plot.data <- function(.reg.output) {
  incentive.treatment.terms <- c("control", "ink", "calendar", "bracelet")
  sms.treatment.interact.terms <- sprintf("%s:social.info", incentive.treatment.terms[-1])
  
  sms.ctrl.linear.restrict <- c("(intercept)", sprintf("(intercept) + %s", incentive.treatment.terms[-1])) 
  social.info.add.effect.restrict <- c("social.info", sprintf("social.info + %s:social.info", incentive.treatment.terms[-1]))
  reminder.only.add.effect.restrict <- c("reminder.only - social.info")
  
  c(sms.ctrl.linear.restrict, social.info.add.effect.restrict, reminder.only.add.effect.restrict) %>% 
    linear_tester(.reg.output, .) %>% 
    mutate(incentive.treatment = incentive.treatment.terms %>% 
             factor(c(rep(., 2), "control"), levels = ., labels = str_to_title(.)),  
           sms.treatment = c(rep(c("sms.control", "social.info"), each = 4), "reminder.only") %>% 
             factor(levels = c("sms.control", "social.info", "reminder.only"), 
                    labels = c("None", "Social Information", "Reminders Only")),
           bar.size = estimate + if_else(sms.treatment == "None", 0, c(rep(estimate[1:4], 2), first(estimate) + estimate[5])),
           ci.lb = bar.size - std.error * 1.64,
           ci.ub = bar.size + std.error * 1.64) 
}

plot.sms.treat.takeup <- function(.preped.data) { 
  ggplot(.preped.data, aes(x = incentive.treatment)) +
    geom_col(aes(y = estimate, fill = sms.treatment), 
             position = position_stack(reverse = TRUE), alpha = 0.25, color = "grey50") +
    geom_crossbar(aes(y = bar.size, ymin = ci.lb, ymax = ci.ub, group = sms.treatment, color = sms.treatment), 
                  width = 0.075, fatten = 0,
                  data = . %>% filter(sms.treatment != "None" | incentive.treatment != "Control")) +
    geom_text_repel(aes(y = bar.size, label = sprintf("%.2f", bar.size)), nudge_x = 0.2, nudge_y = 0.02, segment.colour = "grey50") +
    geom_point(aes(y = bar.size)) +
    scale_y_continuous("Proportion of Take-up", breaks = seq(0, 1, 0.05)) +
    scale_x_discrete("Treatment") +
    scale_fill_discrete("SMS Treatment") +
    scale_color_discrete("SMS Treatment") +
    labs(title = "Estimated Take-up in Response to Incentive and SMS Treatment", 
         caption = "Intervals shown identify the 90% confidence intervals, estimated using cluster robust standard errors.\nRed intervals test the null hypothesis of no difference in take-up from that in the control group (with non SMS treatment).\nGreen and blue intervals test the null hypothesis of no difference from the SMS treatment lower on the same column.") 
}

plot.takeup.dynamics <- function(.data, .aes = aes(x = dewormed.day.any, y = takeup.prop, color = assigned.treatment)) {
  .data %>% 
    ggplot(.aes) +
    geom_line() + 
    geom_point(aes(color = assigned.treatment)) + 
    annotate(geom = "rect", xmin = after.msg.days - 0.5, xmax = after.msg.days + 0.5, ymin = 0, ymax = 1, alpha = 0.2) +
    annotate(geom = "rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 1, alpha = 0.2, fill = "darkred") +
    scale_x_continuous("Deworming Day", breaks = 1:12) +
    scale_color_discrete("Incentive Treatment") +
    coord_cartesian(y = c(0, max(.data[[deparse(.aes$y)]]) + 0.0001)) +
    labs(caption = "Grey vertical bars mark the days after SMS messages are received.\nThe red bar marks the first deworming day, so subjects would have only received a reminder message the day before.")
}

prep.sms.treat.dist.plot.data <- function(.reg.output) {
  interact.with <- "far"
  
  incentive.treatment.terms <- c("control", "ink", "calendar", "bracelet")
  sms.treatment.interact.terms <- sprintf("%s:social.info", incentive.treatment.terms[-1])
  
  sms.ctrl.linear.restrict <- c("(intercept)", sprintf("(intercept) + %s", incentive.treatment.terms[-1])) 
  social.info.add.effect.restrict <- c("social.info", sprintf("social.info + %s:social.info", incentive.treatment.terms[-1]))
  reminder.only.add.effect.restrict <- c("reminder.only - social.info")
  
  col.inference <- c(social.info.add.effect.restrict, reminder.only.add.effect.restrict) %>% 
    c(map(., str_split, "\\s+", simplify = TRUE) %>% 
        map(~ c(.x, "+", str_replace(.x, "([^-+]+)", sprintf("\\1:%s", interact.with)))) %>% 
        map(paste, collapse = " ") %>% 
        # map_if(~ !str_detect(.x, "-"), ~ paste(.x, sprintf("+ %s", interact.with))) %>% 
        unlist) %>% 
    c(sms.ctrl.linear.restrict %>% 
        c(paste(., interact.with, sep = " + ") %>% 
            paste(c("", sprintf("%s:%s", incentive.treatment.terms[-1], interact.with)), sep = " + "))) %>% 
    linear_tester(.reg.output, .) %>% 
    mutate(dist = rep(rep(c("Close", "Far"), 2), c(5, 5, 4, 4)),
           incentive.treatment = incentive.treatment.terms %>% 
             factor(c(rep(c(., "control"), 2), rep(., 2)), levels = ., labels = str_to_title(.))) %>% 
    group_by(dist) %>% 
    mutate(sms.treatment = c(rep("social.info", 4), "reminder.only", rep("sms.control", 4)) %>% 
             factor(levels = c("sms.control", "social.info", "reminder.only"), 
                    labels = c("None", "Social Information", "Reminders Only")),
           bar.size = estimate + if_else(sms.treatment == "None", 0, c(rep(estimate[6:9], 2), 0)),
           bar.size = bar.size + if_else(sms.treatment == "Reminders Only", estimate[1], 0)) %>% 
    ungroup %>% 
    mutate(ci.lb = bar.size - std.error * 1.64,
           ci.ub = bar.size + std.error * 1.64) 
}

# Old plotting code -------------------------------------------------------

prep.ref.group.col <- . %>% {
      intercept.est <- filter(., term == "(intercept)") %$% estimate
      mutate(., 
             ate = estimate,
             estimate = if_else(term != "(intercept)", estimate + intercept.est, estimate))
    } 

incentive.treat.barplot <- function(.data) { 
  .data %<>% 
    filter(ref == "ctrl.ref" | term == "bracelet") %>% 
    mutate(term = factor(term, levels = c("(intercept)", "ink", "calendar", "bracelet"), 
                         labels = c("Control", "Ink", "Calendar", "Bracelet")),
           p.val.label = sprintf("p-value = %.3f", p.value),
           ci.lb = estimate - 1.64 * std.error, 
           ci.ub = estimate + 1.64 * std.error) 
  
  .data %>% 
    ggplot(aes(term, estimate)) +
    geom_col(aes(fill = term %in% c("Control", "Ink")), color = "black", alpha = 0.25, data = . %>% filter(ref == "ctrl.ref")) +
    geom_errorbar(aes(color = ref, ymin = ci.lb, ymax = ci.ub), 
                  width = 0.1, position = "dodge", 
                  data = . %>% filter((term != "Control" & ref == "ctrl.ref") | (ref == "calendar.ref"))) +
    geom_text(aes(y = max(.data$ci.ub) + 0.055, label = p.val.label, color = ref), 
              data = . %>% filter(term != "Control", ref == "ctrl.ref")) +
    geom_text(aes(y = max(.data$ci.ub) + 0.02, label = p.val.label, color = ref), 
              data = . %>% filter(term != "Control", ref == "calendar.ref")) +
    geom_text(aes(label = sprintf("%.3f", estimate)), nudge_x = 0.2, nudge_y = -0.02, data = . %>% filter(ref == "ctrl.ref")) +
    geom_text(aes(y = 0.02, label = sprintf("ATE = %.3f", ate), color = ref), 
              data = . %>% filter(ref == "calendar.ref", term != "Control")) +
    geom_text(aes(y = 0.055, label = sprintf("ATE = %.3f", ate), color = ref), 
              data = . %>% filter(ref == "ctrl.ref", term != "Control")) +
    scale_y_continuous("Proportion of Take-up", breaks = seq(0, 1, 0.05)) +
    scale_x_discrete("Treatment") +
    scale_color_manual("Reference Arm", values = c("red", "black"), labels = c("Calendar", "Control")) +
    scale_fill_manual("", values = c("red", "black"), labels = c("Calendar vs Bracelet", "Control vs Ink")) 
}

# Table Generators --------------------------------------------------------

pval.stars <- . %>% 
  symnum(cutpoints = c(0, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ""))

print.reg.table <- function(.reg.table.data, 
                            caption = "", label = "", font.size = "footnotesize", estimate.buffer = TRUE, landscape = FALSE) {
  all.pt.est <- .reg.table.data %>% 
    select(spec, num.col, pt.est) %>% 
    mutate(spec = forcats::as_factor(spec)) %>% 
    unnest(pt.est) %>%
    mutate(term = fct_relevel(term, 
                              "bracelet", "calendar", "ink", "reminder.only", "social.info", "control",
                              "bracelet:far", "calendar:far", "ink:far", "far",
                              "bracelet:social.info", "calendar:social.info", "ink:social.info"),
           ref.type =  fct_relevel(forcats::as_factor(ref.type),
                                   "control", "ink", "calendar", "bracelet", "reminder.only", "social.info",
                                   "far", "ink:far", "calendar:far", "bracelet:far",
                                   "ink:social.info", "calendar:social.info", "bracelet:social.info")) %>% 
    mutate_at(vars(term, ref.type), 
              funs(fct_relabel(., function(.label) str_replace_all(.label, c(":" = " [x] ", "\\." = " ")) %>% 
                                 str_to_title %>% 
                                 str_replace_all(fixed("[X]"), "$\\times$")))) %>% 
    right_join(modelr::data_grid(., spec, ref.type, term), c("spec", "ref.type", "term")) 
  
  total.num.cols <- sum(.reg.table.data$num.col)
 
  cmidrules <- accumulate(.reg.table.data$num.col, ~ last(.x) + c(1, .y), .init = 1) %>% 
    magrittr::extract(-1) %>% 
    map(~ sprintf("\\cmidrule(r){%d-%d}", .x[1], .x[2])) %>% 
    str_c(collapse = " ") %>% 
    str_c(if_else(estimate.buffer, "\\\\\n", "\n"), collapse = " ")
  
  cat(sprintf("%s\\begin{table}\\centering\\%s%s%s\n", 
              if_else(landscape, "\\begin{landscape}", ""), 
              font.size,
              if_else(nchar(caption) > 0, sprintf("\\caption{%s}", caption), ""),
              if_else(nchar(label) > 0, sprintf("\\label{%s}", label), "")))
  cat("\\begin{threeparttable}\n")
  cat(sprintf("\\begin{tabular}{l*{%d}{c}}\n", total.num.cols))
  cat("\\toprule\n")
  cat(paste(sprintf("& (%d)", seq_len(total.num.cols)), collapse = " "))
  cat("\\\\\n")
  
  cat(paste("\\emph{Treatment Cell} ", 
            .reg.table.data %$% map2(spec, num.col, ~ sprintf("& \\multicolumn{%d}{c}{%s}", .y, .x)) %>% 
              str_c(collapse = " "),
            "\\\\\n", collapse = " "))
 
  cat(cmidrules)
  
  d_ply(all.pt.est, .(term), function(term.rows) {
    term.rows %>% 
      select(ref.pt, estimate, p.value) %>% 
      pmap(function(ref.pt, estimate, p.value) {
        if (is.na(estimate)) {
          return("&")
        } else if (ref.pt) {
          sprintf("& \\cellcolor[gray]{0.9}{\\bf %.4f}", estimate)
        } else {
          sprintf("& $%.4f^{%s}$", estimate, pval.stars(p.value))
        }
      }) %>% 
      str_c(collapse = " ") %>% 
      str_c(term.rows$term[1], ., "\\\\\n", collapse = " ") %>% 
      cat
    
    term.rows %$%
      map2(ref.pt, std.error, function(ref.pt, .se) {
        if (is.na(.se) || ref.pt) {
          return("&")
        } else {
          sprintf("& (%.4f)", .se)
        }
      }) %>%
      str_c(collapse = " ") %>% {
        if (!all(coalesce(term.rows$ref.pt, TRUE)) || !estimate.buffer) {
          str_c(., "\\\\\n", collapse = " ") 
        } else {
          str_c(., "\n", collapse = " ") 
        } 
      } %>% 
      cat
    
    if (estimate.buffer) cat("\\\\\n")
  })
  
  cat(cmidrules)
  
  cat(paste("Observations",  
            .reg.table.data %$% 
              map2(num.obs, num.col, ~ sprintf("& \\multicolumn{%d}{c}{%d}", .y, .x)) %>% 
              str_c(collapse = " "),
            "\\\\\n", collapse = " "))
  
  cat(paste("Mean Obs/Cluster",  
            .reg.table.data %$% 
              map2(average.cluster.num.obs, num.col, ~ sprintf("& \\multicolumn{%d}{c}{%.2f}", .y, .x)) %>% 
              str_c(collapse = " "),
            "\\\\\n", collapse = " "))
  
  cat(paste("Number of Clusters",  
            .reg.table.data %$% 
              map2(num.clusters, num.col, ~ sprintf("& \\multicolumn{%d}{c}{%d}", .y, .x)) %>% 
              str_c(collapse = " "),
            "\\\\\n", collapse = " "))
  
  cat(paste("$R^2$",  
            .reg.table.data %$% 
              map2(r.squared, num.col, ~ sprintf("& \\multicolumn{%d}{c}{%.4f}", .y, .x)) %>% 
              str_c(collapse = " "),
            "\\\\\n", collapse = " "))
  
  cat(paste("$\\text{Adjusted }R^2$",  
            .reg.table.data %$% 
              map2(adj.r.squared, num.col, ~ sprintf("& \\multicolumn{%d}{c}{%.4f}", .y, .x)) %>% 
              str_c(collapse = " "),
            "\\\\\n", collapse = " "))
  
  cat(cmidrules)
 
  cat("\\emph{Joint Tests ($P$-values)} & \\\\\n") 
  # cat("\\emph{($P$-values)} & \\\\\n") 
 
  linear.tests <- .reg.table.data %>%
    select(spec, num.col, joint.tests.res) %>%
    unnest(joint.tests.res) %>%
    d_ply(.(joint.test.type), 
          . %$% 
            str_c(first(joint.test.type), 
                  str_c(map2(p.value, num.col, ~ sprintf("& \\multicolumn{%d}{c}{%.4f}", .y, .x)), collapse = " "),
                  "\\\\\n",
                  collapse = " ") %>% 
            cat)
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\begin{tablenotes}\\footnotesize\n")
 
  c("Reported analysis is from a stratified linear probability model regression. Highlighted cells identify the reference take-up level against which average treatment effects are reported. Each row reports the reference level or average treatment effect in a particular treatment cell. Standard errors computed using cluster robust Huber-White estimators are reported in parentheses under estimates.",
    "Endline survey samples control for ethnicity, gender, age, schooling, material of household floor, and mobile phone ownership.",
    "${}^{***}: P < 0.01, {}^{**}: P < 0.05, {}^*: P < 0.1$") %>% 
    walk(~ cat(sprintf("\\item %s\n", .)))
  
  cat("\\end{tablenotes}\n")
  cat("\\end{threeparttable}\n")
  cat(sprintf("\\end{table}%s\n", if_else(landscape, "\\end{landscape}", "")))
}

print.sms.interact.table <- function(.reg.table.data,
                                     caption = "", label = "", font.size = "footnotesize", estimate.buffer = TRUE, landscape = FALSE) {
  all.pt.est <- .reg.table.data %>% 
    bind_rows(.id = "spec") %>% 
    mutate(index = seq_len(n())) %>% 
    right_join(modelr::data_grid(., spec, linear.test), c("spec", "linear.test")) %>% 
    mutate_at(vars(linear.test), 
              funs(str_replace_all(., c(":" = " <x> ", "\\." = " ")) %>% 
                     str_to_title %>% 
                     str_replace_all(c("<X>" = "$\\\\times$",
                                       "([-\\+=])" = "$\\1$")))) %>% 
    arrange(index) %>% 
    mutate_at(vars(spec, linear.test), funs(. %>% forcats::as_factor())) %>% 
    select(-index)
  
  total.num.cols <- length(.reg.table.data) * 1
  
  cat(sprintf("%s\\begin{table}\\centering\\%s%s%s\n", 
              if_else(landscape, "\\begin{landscape}", ""), 
              font.size,
              if_else(nchar(caption) > 0, sprintf("\\caption{%s}", caption), ""),
              if_else(nchar(label) > 0, sprintf("\\label{%s}", label), "")))
  cat("\\begin{threeparttable}\n")
  cat(sprintf("\\begin{tabular}{l*{%d}{c}}\n", total.num.cols))
  cat("\\toprule\n")
  
  all.pt.est %$% map_chr(unique(spec), ~ sprintf("& \\multicolumn{1}{c}{%s}", .x)) %>% 
    c(rep("\\\\\n", 1)) %>% 
    str_c(collapse = " ") %>% 
    cat
  
  seq(2, total.num.cols + 1) %>% 
    map(~ sprintf("\\cmidrule(r){%d-%d}", .x, .x)) %>% 
    str_c(collapse = " ") %>% 
    str_c(if_else(estimate.buffer, "\\\\\n", "\n"), collapse = " ") %>% 
    cat
 
  all.pt.est %>% 
    d_ply(.(linear.test), function(test.df) {
      sprintf("%s %s \\\\\n", 
              test.df$linear.test[1], 
              alply(test.df, 1, . %$% sprintf("& $%.4f^{%s}$", estimate, pval.stars(p.value))) %>% 
                str_c(collapse = " "), collapse = " ") %>% 
        cat
      
      alply(test.df, 1, . %$% sprintf("& (%.4f)", std.error)) %>% 
        c(rep("\\\\\n", 1 + 1*(estimate.buffer))) %>% 
        str_c(collapse = " ") %>% 
        cat
    })
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\begin{tablenotes}\\footnotesize\n")
 
  c("Reported analysis is from a stratified linear probability model regression. Each row reports the results of linear hypothesis test. Standard errors computed using cluster robust Huber-White estimators are reported in parentheses under estimates.",
    "Endline survey samples control for ethnicity, gender, age, schooling, material of household floor, and mobile phone ownership.",
    "${}^{***}: P < 0.01, {}^{**}: P < 0.05, {}^*: P < 0.1$") %>% 
    walk(~ cat(sprintf("\\item %s\n", .)))
  
  cat("\\end{tablenotes}\n")
  cat("\\end{threeparttable}\n")
  cat(sprintf("\\end{table}%s\n", if_else(landscape, "\\end{landscape}", "")))
}

# Bayesian Analysis -------------------------------------------------------

prepare_bayesian_analysis_data <- function(origin_prepared_analysis_data, 
                                           wtp_data,
                                           ...,
                                           # treatment_formula = ~ assigned.treatment * dist.pot.group * (phone_owner + sms.treatment.2),
                                           # treatment_formula = ~ assigned.treatment * dist.pot.group * sms.treatment.2,
                                           treatment_col = c("assigned.treatment", "dist.pot.group", "sms.treatment.2", "hh.baseline.sample"),
                                           treatment_map = NULL,
                                           treatment_formula = NULL,
                                           subgroup_col = "phone_owner",
                                           all_ate = NULL,
                                           # exclude_from_eval = NULL, #  "name_matched",
                                           endline_covar = c("ethnicity", "floor", "school")) {
  prep_data_arranger <- function(prep_data, ...) prep_data %>% arrange(stratum_id, new_cluster_id, name_matched, dewormed.any, ...)
  
  prepared_analysis_data <- origin_prepared_analysis_data %>% 
    mutate(new_cluster_id = factor(cluster.id) %>% as.integer(),
           age = if_else(!is.na(age), age, age.census),
           age_squared = age^2,
           missing_covar = is.na(floor)) %>% 
    # unite(county_phone_owner_stratum, county, phone_owner, remove = FALSE) %>%
    mutate(#county_phone_owner_stratum = factor(county_phone_owner_stratum),
           stratum = county) %>% 
    mutate_at(vars(county_dist_stratum, county_dist_mon_stratum, county, stratum, gender, school, floor, ethnicity), 
              funs(id = as.integer)) %>% 
    prep_data_arranger() 

  scale_covar <- function (covar_data) {
    counted_map <- "n" %in% names(covar_data)
    
    covar_data %<>% 
      select_if(~ n_distinct(.) > 1) %>% 
      na.omit()
    
    if (counted_map) {
      map_count <- covar_data$n 
      
      covar_data %<>%
        select(-n) 
    } else {
      map_count <- rep(1, nrow(covar_data))
    }
    
    is_factor_col <- map(covar_data, 
                         ~ switch(class(.),
                                  factor = rep(TRUE, nlevels(.) - 1),
                                  logical = TRUE,
                                  FALSE)) %>% 
      unlist()
    
    design_matrix <- model_matrix(covar_data, ~ .) %>% 
      magrittr::extract(, -1) # get rid of intercept column
   
    map_means <- map2_dbl(design_matrix, is_factor_col, ~ if_else(!.y, weighted.mean(.x, map_count), 0))
   
    design_matrix %>%  
      scale(scale = map2_dbl(., map_means, 
                             # function(col, col_mean) (sum(map_count * (col - col_mean))^2)/(sum(map_count) - 1) * 2) %>% 
                             function(col, col_mean) sqrt(sum(map_count * ((col - col_mean)^2))/(sum(map_count) - 1)) * 2) %>%  
              if_else(is_factor_col, 1, .),
      # scale(scale = map2_dbl(., is_factor_col, ~ if_else(!.y, .x -  sd(.x) * 2, 1)),
            center = map_means) # Center and scale to have SD = 0.5
  }
  
  remove_dup_treatment_dm_cols <- TRUE
 
  if (is.null(treatment_map)) { 
    treatment_map <- expand_(prepared_analysis_data, c(treatment_col, subgroup_col)) %>% {
        if ("phone_owner" %in% names(.)) arrange(., phone_owner) else return(.) 
      } 
  } else {
    treatment_col <- names(treatment_map)
    # remove_dup_treatment_dm_cols <- FALSE
  }
  
  treatment_map %<>% 
    mutate(all_treatment_id = seq_len(n())) %>% 
    mutate_if(is.factor, funs(id = as.integer(.)))
  
  census_covar_map <- count_(prepared_analysis_data, c("age", "gender", subgroup_col)) %>% 
    mutate(census_covar_id = seq_len(n())) 
    
  prepared_analysis_data %<>% 
    left_join(treatment_map, c(treatment_col, subgroup_col)) %>% 
    left_join(select(census_covar_map, -n), c("age", "gender", subgroup_col)) %>% 
    prep_data_arranger() %>% 
    mutate(obs_index = seq_len(n()))
 
  # Get rid of treatment cells that don't exist in the data 
  # treatment_map %<>% 
  #   semi_join(prepared_analysis_data, "all_treatment_id") %>% 
  #   arrange(all_treatment_id) %>% 
  #   mutate(new_all_treatment_id = seq_len(n()))
 
  # prepared_analysis_data %<>%
  #   left_join(select(treatment_map, ends_with("all_treatment_id")), "all_treatment_id") %>% 
  #   select(-all_treatment_id) %>% 
  #   rename(all_treatment_id = new_all_treatment_id) %>% 
  #   arrange(obs_index)
  # 
  # treatment_map %<>%
  #   select(-all_treatment_id) %>% 
  #   rename(all_treatment_id = new_all_treatment_id) 
 
  if (is.null(treatment_formula)) { 
    treatment_formula <- str_c("~ ", str_c(c(treatment_col, subgroup_col), collapse = " * ")) %>% {
        if (!is.null(subgroup_col)) str_c(., " - ", str_c(subgroup_col, collapse = " - ")) else return(.)
      } %>% 
      as.formula()
  }
  
  treatment_map_design_matrix <- treatment_map %>%  
    model_matrix(treatment_formula) %>% 
    magrittr::extract(, -1) %>% # get rid of intercept column
    magrittr::extract(, map_lgl(., ~ n_distinct(.) > 1)) %>% {
      if (!remove_dup_treatment_dm_cols) return(.) else magrittr::extract(., , !duplicated(t(.))) # Remove redundant columns
    }
  
  # experiment_coef <- treatment_map_design_matrix %>% 
  #   names() %>% 
  #   str_detect("name_matched") %>% 
  #   not() %>% 
  #   which()
  
  private_value_calendar_coef <- treatment_map_design_matrix %>% 
    names() %>% 
    str_detect("private_valuecalendar") %>% 
    which()
  
  private_value_bracelet_coef <- treatment_map_design_matrix %>% 
    names() %>% 
    str_detect("private_valuebracelet") %>% 
    which()
  
  # name_match_interact_formula <- str_c(c(treat_variables$interact, treat_variables$direct_only), collapse = "*") %>% 
  #   str_c("~", .) %>% 
  #   as.formula()
  # 
  # name_match_interact_map_design_matrix <- treatment_map %>%  
  #   model_matrix(name_match_interact_formula) %>%
  #   magrittr::extract(, -1) %>% # get rid of intercept column
  #   select(-one_of(names(treatment_map_design_matrix)))
  
  census_covar_map_dm <- census_covar_map %>% 
    mutate(age_squared = age ^ 2) %>% 
    select(-census_covar_id) %>% 
    scale_covar()
  
  census_covar_map %<>% select(-n)
  
  endline_covar_dm <- prepared_analysis_data %>% 
    filter(!missing_covar) %>% 
    prep_data_arranger(obs_index) %>% 
    select_(.dots = endline_covar) %>%
    scale_covar()
  
  # missing_treatment <- prepared_analysis_data %>% 
  #   ddply(., .(all_treatment_id), 
  #         function(treatment_group, all_obs) { 
  #           filter(all_obs, 
  #                  all_treatment_id != treatment_group$all_treatment_id[1],
  #                  phone_owner == treatment_group$phone_owner[1]) %>%
  #             select(stratum_id, obs_index)
  #         }, 
  #         all_obs = .)
  # 
  # eval_treatment_prop_id <- treatment_map %>% {
  #     if("name_matched" %in% names(.)) filter(., !name_matched) else return(.) 
  #   } %>% 
  #   pull(all_treatment_id)
  
  name_matched_data <- filter(prepared_analysis_data, name_matched == 1) %>% 
    prep_data_arranger() # Should already be in the right order, but to be on the safe size
  
  monitored_data <- filter(prepared_analysis_data, name_matched == 0) %>% 
    prep_data_arranger() # Should already be in the right order, but to be on the safe size
  
  matching_error_data <- prepared_analysis_data %>% 
    filter(#assigned.treatment == "control", 
           sms.treatment.2 == "sms.control",
           true.monitored | !monitored) # using only true monitored and those not monitored because they weren't in a study sample
    # group_by(stratum_id, new_cluster_id, name_matched) %>% 
    # summarize(dewormed_prop = sum(dewormed.any)/n()) %>% 
    # ungroup()
     
  stratum_map <- distinct(prepared_analysis_data, stratum_id, stratum)
  
  # calendar_treated <- prepared_analysis_data %>% 
  #   filter(assigned.treatment == "calendar") %>% 
  #   arrange(stratum_id, obs_index)
  
  bracelet_treated <- prepared_analysis_data %>%
    filter(assigned.treatment == "bracelet") %>%
    arrange(stratum_id, obs_index)

  incentive_choice_data <- origin_prepared_analysis_data %>%
    filter(!is.na(gift_choice) & gift_choice != "neither",
           assigned.treatment == "control",
           sms.treatment.2 == "sms.control") %>%
    select(county, gift_choice) %>%
    mutate(offer = 0,
           response = "keep")

  incentive_choice_data <- wtp_data %>%
    filter(!is.na(first_choice)) %>%
    transmute(county,
              gift_choice = first_choice,
              offer = price,
              response = second_choice) %>%
    bind_rows(incentive_choice_data) %>%
    left_join(stratum_map, c("county" = "stratum")) %>%
    mutate(gift_choice = 2 * (gift_choice == "calendar") - 1,
           response = 2 * (response == "switch") - 1) %>%
    arrange(stratum_id, response)
 
  if (!is.null(all_ate)) {
    left_right_cells_colnames <- str_subset(names(all_ate), "_(left|right)$")
    both_cells_colnames <- setdiff(names(all_ate), left_right_cells_colnames)
    left_right_cells_colnames %<>% 
      str_replace("_(left|right)$", "") %>% 
      unique()
    
    all_ate %<>% 
      inner_join(treatment_map, 
                c(left_right_cells_colnames %>% setNames(paste0(., "_left")), both_cells_colnames)) %>% 
      inner_join(treatment_map, 
                c(left_right_cells_colnames %>% setNames(paste0(., "_right")), both_cells_colnames),
                suffix = c("_left", "_right"))
    
    non_phone_owner_treatments <- all_ate %>%
      filter(!phone_owner) %>%
      select(all_treatment_id_left, all_treatment_id_right) %>%
      gather(value = id) %>%
      distinct(id) %>%
      arrange(id) %>%
      mutate(rank_id = seq_len(n()))

    phone_owner_treatments <- all_ate %>%
      filter(phone_owner) %>%
      select(all_treatment_id_left, all_treatment_id_right) %>%
      gather(value = id) %>%
      # bind_rows(all_subgroup_treatments) %>% 
      distinct(id) %>%
      arrange(id) %>%
      mutate(rank_id = seq_len(n()))
    
    non_phone_owner_missing_treatment <- non_phone_owner_treatments$id %>% 
      map(~ filter(prepared_analysis_data, !phone_owner, all_treatment_id != .x) %>% pull(obs_index))
    
    non_phone_owner_observed_treatment <- non_phone_owner_treatments$id %>% 
      map(~ filter(prepared_analysis_data, !phone_owner, all_treatment_id == .x) %>% pull(obs_index))
    
    phone_owner_missing_treatment <- phone_owner_treatments$id %>% 
      map(~ filter(prepared_analysis_data, phone_owner, all_treatment_id != .x) %>% pull(obs_index))
    
    phone_owner_observed_treatment <- phone_owner_treatments$id %>% 
      map(~ filter(prepared_analysis_data, phone_owner, all_treatment_id == .x) %>% pull(obs_index))
    
    num_non_phone_owner_ate_pairs <- all_ate %>% filter(!phone_owner) %>% nrow()
    num_phone_owner_ate_pairs <- all_ate %>% filter(phone_owner) %>% nrow()
    
    non_phone_owner_ate_pairs <- all_ate %>% 
      filter(!phone_owner) %>% 
      select(all_treatment_id_left, all_treatment_id_right) %>%
      left_join(non_phone_owner_treatments, c("all_treatment_id_left" = "id")) %>% 
      left_join(non_phone_owner_treatments, c("all_treatment_id_right" = "id"), suffix = c("_left", "_right")) %>% 
      select(rank_id_left, rank_id_right)
    
    phone_owner_ate_pairs <- all_ate %>% 
      filter(phone_owner) %>% 
      select(all_treatment_id_left, all_treatment_id_right) %>% 
      left_join(phone_owner_treatments, c("all_treatment_id_left" = "id")) %>% 
      left_join(phone_owner_treatments, c("all_treatment_id_right" = "id"), suffix = c("_left", "_right")) %>% 
      select(rank_id_left, rank_id_right)
  } else {
    non_phone_owner_treatments <- NULL
    phone_owner_treatments <- NULL
    
    non_phone_owner_missing_treatment <- NULL 
    non_phone_owner_observed_treatment <- NULL
    phone_owner_missing_treatment <- NULL
    phone_owner_observed_treatment <- NULL
    
    num_non_phone_owner_ate_pairs <- NULL
    num_phone_owner_ate_pairs <- NULL
    
    non_phone_owner_ate_pairs <- NULL
    phone_owner_ate_pairs <- NULL
  }
  
  lst(
    # Save meta data 
    
    prepared_analysis_data,
     
    treatment_map,
    
    all_ate,
    
    non_phone_owner_missing_treatment,
    non_phone_owner_observed_treatment,
    phone_owner_missing_treatment,
    phone_owner_observed_treatment,
    
    # missing_treatment,
   
    census_covar_map, 
    stratum_map,
    cluster_map = distinct(prepared_analysis_data, new_cluster_id, cluster.id),
    
    # Data passed to model
    
    treatment_map_design_matrix,
    # name_match_interact_map_design_matrix,
    # experiment_coef,
    num_all_treatments = nrow(treatment_map_design_matrix),
    num_all_treatment_coef = ncol(treatment_map_design_matrix), 
    # num_experiment_coef = length(experiment_coef),
    
    num_private_value_calendar_coef = length(private_value_calendar_coef),
    num_private_value_bracelet_coef = length(private_value_bracelet_coef),
    private_value_calendar_coef,
    private_value_bracelet_coef,
    not_private_value_bracelet_coef = seq_len(num_all_treatment_coef) %>% setdiff(private_value_bracelet_coef),
    # num_not_private_value_bracelet_coef = length(not_private_value_bracelet_coef),
    
    num_bracelet_treated = nrow(bracelet_treated),
    bracelet_treated_id = bracelet_treated$obs_index,
    strata_bracelet_sizes = if (num_bracelet_treated > 0) {
      bracelet_treated %>% count(stratum_id) %>% arrange(stratum_id) %>% pull(n)
    } else {
      rep(0, num_strata)
    },
    
    # num_name_match_interact_coef = ncol(name_match_interact_map_design_matrix),
    name_matched = prepared_analysis_data$name_matched,
    num_name_matched = sum(name_matched),
    name_matched_id = name_matched_data %>% pull(obs_index), 
    monitored_id = monitored_data %>% pull(obs_index), 
    name_matched_strata_sizes = count(name_matched_data, stratum_id) %>% arrange(stratum_id) %$% n,
    name_matched_dewormed_strata_sizes = filter(name_matched_data, dewormed.any) %>% 
      count(stratum_id) %>% 
      arrange(stratum_id) %$% n,
    name_matching_error_ids = matching_error_data %>% arrange(obs_index) %$% obs_index,
    num_name_matching_errors_ids = length(name_matching_error_ids),
    
    phone_owner_indices = prepared_analysis_data$phone_owner + 1,
    
    census_covar_map_dm,
    num_census_covar_coef = ncol(census_covar_map_dm),
    num_distinct_census_covar = nrow(census_covar_map_dm),
    census_covar_id = prepared_analysis_data$census_covar_id,
    
    endline_covar_dm,
    num_endline_covar_coef = ncol(endline_covar_dm),
    
    stratum_covar_id = prepared_analysis_data %>% prep_data_arranger(missing_covar, obs_index) %$% obs_index, # First obs then missing
    stratum_missing_covar_sizes = prepared_analysis_data %>% 
      filter(missing_covar | (select_(., .dots = endline_covar) %>% map(is.na) %>% reduce(or))) %>% 
      count(stratum_id) %>% 
      arrange(stratum_id) %$% 
      n,
    
    obs_treatment = prepared_analysis_data$all_treatment_id,
  
    treatment_id = prepared_analysis_data %>% arrange(all_treatment_id, obs_index) %$% obs_index,
    
    num_obs = length(obs_treatment), 
    
    treatment_sizes = count(prepared_analysis_data, all_treatment_id) %>% arrange(all_treatment_id) %>% pull(n),
    
    dewormed_any = prepared_analysis_data$dewormed.any,
    
    cluster_id = prepared_analysis_data$new_cluster_id,
    stratum_id = prepared_analysis_data$stratum_id,
    num_clusters = n_distinct(cluster_id),
    
    num_strata = n_distinct(stratum_id),
    strata_sizes = count(prepared_analysis_data, stratum_id) %>% arrange(stratum_id) %>% pull(n),
    cluster_sizes = count(prepared_analysis_data, stratum_id, new_cluster_id) %>% arrange(stratum_id, new_cluster_id) %>% pull(n),
    strata_num_clusters = distinct(prepared_analysis_data, stratum_id, new_cluster_id) %>% 
      count(stratum_id) %>% 
      arrange(stratum_id) %>% 
      pull(n),
    
    # ATE
    
    num_non_phone_owner_treatments = nrow(non_phone_owner_treatments),
    num_phone_owner_treatments = nrow(phone_owner_treatments),
    
    non_phone_owner_treatments = non_phone_owner_treatments$id,
    phone_owner_treatments = phone_owner_treatments$id,
    
    missing_non_phone_owner_obs_ids = unlist(non_phone_owner_missing_treatment),
    num_missing_non_phone_owner_obs_ids = length(missing_non_phone_owner_obs_ids),
    missing_non_phone_owner_treatment_sizes = if (!is.null(all_ate)) map_int(non_phone_owner_missing_treatment, length),
    observed_non_phone_owner_obs_ids = unlist(non_phone_owner_observed_treatment),
    num_observed_non_phone_owner_obs_ids = length(observed_non_phone_owner_obs_ids),
    observed_non_phone_owner_treatment_sizes = if (!is.null(all_ate)) map_int(non_phone_owner_observed_treatment, length),
    
    missing_phone_owner_obs_ids = unlist(phone_owner_missing_treatment),
    num_missing_phone_owner_obs_ids = length(missing_phone_owner_obs_ids), 
    missing_phone_owner_treatment_sizes = if (!is.null(all_ate)) map_int(phone_owner_missing_treatment, length),
    observed_phone_owner_obs_ids = unlist(phone_owner_observed_treatment),
    num_observed_phone_owner_obs_ids = length(observed_phone_owner_obs_ids),
    observed_phone_owner_treatment_sizes = if (!is.null(all_ate)) map_int(phone_owner_observed_treatment, length),
    
    num_non_phone_owner_ate_pairs,
    num_phone_owner_ate_pairs,
    
    non_phone_owner_ate_pairs,
    phone_owner_ate_pairs,
    
    # WTP data
    
    num_wtp_obs = nrow(incentive_choice_data),
    wtp_strata_sizes = count(incentive_choice_data, stratum_id) %>% arrange(stratum_id) %>% pull(n),
    gift_choice = incentive_choice_data$gift_choice,
    wtp_offer = incentive_choice_data$offer,
    wtp_response = incentive_choice_data$response,
    
    ...
  )
} 
