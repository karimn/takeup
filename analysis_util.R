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

prepare.analysis.data <- function(.census.data, .takeup.data, .endline.data, .baseline.data, .consent.dewormed.reports, .cluster.strat.data,
                                  max.name.match.cost = 1) {
  dewormed.day.data <- .takeup.data %>% 
    filter(!is.na(KEY.individ)) %>% 
    rename(dewormed.day = deworming.day) %>% 
    group_by(KEY.individ) %>% 
    summarize_at(vars(dewormed.day, dewormed.date), min) %>% # If dewormed multiple times, take the first day only
    ungroup
  
  analysis.data <- .census.data %>% 
    filter(!is.na(wave)) %>%  # Remove clusters no longer in study
    left_join(.consent.dewormed.reports, "KEY.individ") %>% 
    mutate(dewormed = KEY.individ %in% discard(.takeup.data$KEY.individ, is.na), # TRUE if individual found in take-up data
           dewormed = ifelse(monitored, dewormed, NA),
           # hh.baseline.sample = KEY %in% .baseline.data$KEY, # We don't have a clear link between baseline surveys and the census. See the field notebook
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
    right_join(analysis.data, c("cluster.id", "KEY.individ")) %>% 
    left_join(transmute(takeup.data, KEY.survey.individ, dewormed.day.matched = deworming.day, dewormed.date.matched = dewormed.date), 
              c("which.min.name.match.dist" = "KEY.survey.individ")) %>% 
    mutate(monitored = !is.na(monitored) & !is.na(wave) & monitored, # Remove those dropped from the study 
           dewormed.any = (!is.na(dewormed) & dewormed) | dewormed.matched,
           dewormed.day.any = if_else(!is.na(dewormed.day), as.integer(dewormed.day), dewormed.day.matched), 
           dewormed.date.any = if_else(!is.na(dewormed.date), dewormed.date, dewormed.date.matched),
           gender = factor(gender, levels = 1:2, labels = c("male", "female"))) %>% 
    left_join(
        transmute(.endline.data, 
                  KEY.individ, age, school, floor, ethnicity, ethnicity2, any.sms.reported, gift_choice, hh_cal, cal_value, hh_bracelet, number_bracelet, 
                  endline_deworm_rate = dworm_rate), "KEY.individ") %>% 
    mutate_at(vars(age, age.census), funs(squared = (.)^2, group = cut(., breaks = c(seq(18, 58, 10), 120), right = FALSE))) %>% 
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
          !(hh.baseline.sample.pool & exclude.baseline.sample)) %>% 
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
       Baseline = .baseline.data %>% { if (!plyr::empty(.)) mutate(., KEY.individ = KEY) else return(.) }, 
       Endline = .endline.data) %>% 
    compact %>% 
    map_df(~ select_(.x, .dots = c(intersect(names(.x), question.info$col.name), "KEY.individ")), .id = "survey.type") %>% 
    map_df(question.info$col.name, 
           function(.col.name, .data) {
             .data %>% {
                 if (is.list(.[[.col.name]])) { 
                   if (is.factor(.[[.col.name]][[1]])) {
                     .[[.col.name]] <- map(.[[.col.name]], as.character)
                   }
                   
                   unnest_(., .col.name) 
                 } else return(.)
               } %>% {
                 if(na.rm) filter_(., sprintf("!is.na(%s)", .col.name)) else return(.)
               } %>%
               group_by(survey.type) %>% 
               do(mutate(count_(., .col.name), n = n/n_distinct(.$KEY.individ))) %>% 
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
    filter(!is.na(gift_choice), monitored, monitor.consent, !hh.baseline.sample.pool, !is.na(sms.treatment)) %>% 
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

plot_dyn_takeup <- function(takeup_summ_data, takeup_data = NULL, combiner = NULL, data_preparer = function(data) data) {
  inner_data_preparer <- . %>% 
    mutate_at(vars(starts_with("incentive_treatment")), funs(fct_relabel(., str_to_title))) %>% 
    data_preparer() 
  
  takeup_data %<>% 
      inner_data_preparer()
  
  if (!is_null(combiner)) {
    takeup_data %<>% combiner()
  } 
    
  takeup_summ_data %>%
    inner_data_preparer() %>% { 
      plot_obj <- ggplot(., aes(incentive_treatment_static, mean_est)) 
      caption_text <-
        "Points represent mean point estimates and circles represent observed take-up levels.
         The thick and thin vertical lines show the 90% and 95% posterior probability ranges, respectively."
      
      if (!is_null(takeup_data)) {
        plot_obj <- plot_obj + geom_violin(aes(y = wtd_iter_est), 
                                           # draw_quantiles = c(0.25, 0.5, 0.75), color = "white", fill = "darkgrey", data = takeup_data) 
                                           draw_quantiles = c(0.25, 0.5, 0.75), color = "lightgrey", fill = "darkgrey", data = takeup_data)
        
        caption_text %<>% str_c("Vertical lines in the posterior distribution density identify the 25th, 50th, and 75th percentiles.", 
                                sep = "\n")
      } 
      
      plot_obj +
        # ggplot(aes(incentive_treatment_static, mean_est, color = sms.treatment.2_static)) +
        geom_pointrange(aes(ymin = lb_95, ymax = ub_95), size = 0.5, position = position_dodge(width = 0.5)) +
        geom_linerange(aes(ymin = lb_90, ymax = ub_90), size = 1.5, position = position_dodge(width = 0.5)) +
        geom_text_repel(aes(label = sprintf("%.3f", mean_est), color = NULL), nudge_x = -0.25, size = 3, segment.color = NA) +
        geom_point(aes(y = observed_takeup_prop), alpha = 0.5, shape = 21, stroke = 1.5, size = 5, position = position_dodge(width = 0.5)) +
        # scale_color_discrete("SMS Treatment", labels = c("None", "Reminders Only", "Social Information")) +
        coord_flip() +
        labs(title = "Deworming Take-up Levels", 
             # subtitle = "Across all treatments", 
             x = "Incentive/Signal Treatment", 
             y = "Probability of Deworming", 
             caption = caption_text) + 
        theme(legend.position = "right") 
    }
}

plot_dyn_takeup_daily <- function(daily_takeup_summ) {
  daily_takeup_summ %>% 
    filter(day < 13) %>% 
    mutate_at(vars(starts_with("incentive_treatment_static")), funs(fct_relabel(., str_to_title))) %>% 
    ggplot(aes(day, mean_est)) +
    geom_line(aes(group = incentive_treatment_static, color = incentive_treatment_static)) +
    scale_x_continuous("Deworming Day", breaks = 1:12, minor_breaks = FALSE) +
    scale_color_discrete("Incentive/Signal Treatment") +
    labs(title = "Daily Mean Probability of Deworming Take-up", y = "Posterior Mean Probability of Deworming") 
}

plot_dyn_ate <- function(ate_summ_data, ate_data = NULL, combiner = NULL, data_preparer = function(data) data) {
  inner_data_preparer <- . %>% 
    mutate_at(vars(starts_with("incentive_treatment")), funs(fct_relabel(., str_to_title))) %>% 
    data_preparer() 
  
  ate_data %<>% 
      inner_data_preparer()
  
  if (!is_null(combiner)) {
    ate_data %<>% combiner()
  } 
  
  # mutate(ref = if_else(incentive_treatment_right_static == "control" & sms.treatment.2_right == "sms.control", 
  #                      "control-sms.control", 
  #                      if_else(sms.treatment.2_right == "sms.control", "own sms.control", "control-social.info")),
  #        ref = fct_relevel(ref, "control-sms.control", "own sms.control", "control-social.info")) %>%
  ate_summ_data %>% 
    inner_data_preparer() %>% { 
      plot_obj <- ggplot(., aes(incentive_treatment_static_left, mean_est)) 
        # ggplot(aes(incentive_treatment_static, mean_est, color = sms.treatment.2_static)) +
      caption_text <-
        "Points represent mean point estimates.
         The thick and thin vertical lines show the 90% and 95% posterior probability ranges, respectively."
        
      if (!is_null(ate_data)) {
        plot_obj <- plot_obj + geom_violin(aes(y = wtd_iter_est), 
                                           draw_quantiles = c(0.25, 0.5, 0.75), color = "lightgrey", fill = "darkgrey", data = ate_data)
        
        caption_text %<>% str_c("Vertical lines in the posterior distribution density identify the 25th, 50th, and 75th percentiles.", 
                                sep = "\n")
      }
      
      plot_obj +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_pointrange(aes(ymin = lb_95, ymax = ub_95), size = 0.5, position = position_dodge(width = 0.5)) +
        geom_linerange(aes(ymin = lb_90, ymax = ub_90), size = 1.5, position = position_dodge(width = 0.5)) +
        geom_text_repel(aes(label = sprintf("%.3f", mean_est), color = NULL), nudge_x = -0.25, size = 3, segment.color = NA) +
        scale_y_continuous(breaks = seq(-1, 1, 0.05)) +
        coord_flip() +
        labs(title = "Deworming Take-up Average Treatment Effect", 
             # subtitle = "Across all treatments", 
             x = "Incentive/Signal Treatment", 
             y = "Posterior Probability of Deworming Average Treatment Effect", 
             caption = caption_text) + 
        theme(legend.position = "right") 
    }
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
  
  plyr::d_ply(all.pt.est, .(term), function(term.rows) {
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
    plyr::d_ply(.(joint.test.type), 
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
    plyr::d_ply(.(linear.test), function(test.df) {
      sprintf("%s %s \\\\\n", 
              test.df$linear.test[1], 
              plyr::alply(test.df, 1, . %$% sprintf("& $%.4f^{%s}$", estimate, pval.stars(p.value))) %>% 
                str_c(collapse = " "), collapse = " ") %>% 
        cat
      
      plyr::alply(test.df, 1, . %$% sprintf("& (%.4f)", std.error)) %>% 
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

prepare_bayes_wtp_data <- function(
  origin_prepared_analysis_data, 
  wtp_data,
  stratum_map = origin_prepared_analysis_data %>% 
    mutate(stratum_id = as.integer(stratum)) %>% 
    distinct(stratum_id, stratum),
  ...) {
  
  incentive_choice_data <- origin_prepared_analysis_data %>%
    filter(!is.na(gift_choice) & gift_choice != "neither",
           assigned.treatment == "control",
           sms.treatment.2 == "sms.control") %>%
    select(county, gift_choice, phone_owner) %>%
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
 
  lst(
    num_strata = nrow(stratum_map),
    
    num_wtp_obs = nrow(incentive_choice_data),
    wtp_strata_sizes = count(incentive_choice_data, stratum_id) %>% arrange(stratum_id) %>% pull(n),
    gift_choice = incentive_choice_data$gift_choice,
    wtp_offer = incentive_choice_data$offer,
    wtp_response = incentive_choice_data$response,
    
    ...
  ) 
}

identify_treatment_id <- function(ate_pairs, treatment_map) {
  left_right_cells_colnames <- str_subset(names(ate_pairs), "_(left|right)$")
  both_cells_colnames <- setdiff(names(ate_pairs), left_right_cells_colnames)
  left_right_cells_colnames %<>% 
    str_replace("_(left|right)$", "") %>% 
    unique()
  
  ate_pairs %>% 
    inner_join(treatment_map, 
              c(left_right_cells_colnames %>% setNames(paste0(., "_left")), 
                both_cells_colnames)) %>% 
    inner_join(treatment_map, 
              c(left_right_cells_colnames %>% setNames(paste0(., "_right")),
                both_cells_colnames),
              suffix = c("_left", "_right")) 
}

get_unique_treatments <- function(ate_pairs) {
  ate_pairs %>% 
    select(matches("(all|dynamic)_treatment_id_(left|right)$")) %>%
    mutate(row_id = seq_len(n())) %>%
    gather(key = treatment_key, value = id, matches("_(left|right)$")) %>% 
    separate(treatment_key, c("treatment_key", "left_right"), "_(?=left|right)") %>% 
    unite(row_id, row_id, left_right) %>% 
    spread(treatment_key, id) %>% 
    distinct_(.dots = str_subset(names(.), "(all|dynamic)_treatment_id$")) %>% 
    arrange_(.dots = str_subset(names(.), "(all|dynamic)_treatment_id$")) %>% 
    mutate(rank_id = seq_len(n()))
}

prepare_dynamic_treatment_maps <- function(dynamic_treatment_map_config, prepared_analysis_data, drop_intercept = TRUE, day_interactions_only = TRUE) {
  dyn_var <- names(dynamic_treatment_map_config$trends)
  dyn_var_re <- str_c(dyn_var, collapse = "|")
  dyn_formula_var <- all.vars(dynamic_treatment_map_config$formula)
  original_dyn_var <- setdiff(dyn_formula_var, dyn_var)
  
  dynamic_treatment_mask_map <- prepared_analysis_data %>% 
    expand(crossing_(original_dyn_var), dynamic_treatment_map_config$trends) 
 
  dynamic_treatment_map <- dynamic_treatment_mask_map %>% 
    model_matrix(dynamic_treatment_map_config$formula) %>% {
      if (drop_intercept) {
        select(., matches(dyn_var_re), -contains("control")) 
      } else {
        select(., matches(dyn_var_re), -contains("control"), "(Intercept)") 
      }
    } %>% {
      if (day_interactions_only) select(., matches(str_c(original_dyn_var, collapse = "|"))) else return(.) 
    } %>%
    select(one_of(dyn_var), which(str_detect(names(.), ":(?!phone_owner)"))) %>% {
      if (dynamic_treatment_map_config$scale %||% FALSE) scale(., center = FALSE, scale = rep(max(.), ncol(.))) else return(.)
      # if (dynamic_treatment_map_config$scale %||% FALSE) map_dfc(., ~ scale(.x, center = FALSE, scale = max(.x))) else return(.)
    } %>% 
    as_tibble()
  
  trend_length <- nrow(dynamic_treatment_map_config$trends)
  
  dynamic_treatment_mask_map %<>% 
    mutate(dynamic_treatment_id = rep(seq_len(nrow(.) %/% trend_length), each = trend_length))
  #   {
  #   mask <- select(., which(!str_detect(names(.), fixed("phone_owner")) & str_detect(names(.), fixed(original_dyn_var)))) 
  #   
  #   if (!is_empty(mask)) {
  #     mask %<>% equals(0) %>% not() %>% multiply_by(1)
  #     
  #     cbind(select(., -one_of(colnames(mask))), mask) %>%
  #       left_join(distinct_(., .dots = original_dyn_var) %>%
  #                   mutate(dynamic_treatment_id = seq_len(n())), original_dyn_var) %>%
  #       arrange(dynamic_treatment_id)
  #   } else {
  #     return(NULL)
  #   }
  # }
  
  dynamic_treatment_map %<>% {
    if (is_null(dynamic_treatment_mask_map)) {
      stop("Dynamic IDs are wrong.")
      mutate(., dynamic_treatment_id = seq_len(n()))
    } else {
      bind_cols(., dynamic_treatment_mask_map %>% select(dynamic_treatment_id))
    }
  }
    # mutate(dynamic_treatment_id = seq_len(n()))
    
  lst(map = dynamic_treatment_map, 
      mask = dynamic_treatment_mask_map,
      original_dyn_var)
}

prepare_treatment_map <- function(static_treatment_map, dynamic_treatment_maps) { #, ate) {
  prepared <- static_treatment_map %>% 
    mutate(all_treatment_id = seq_len(n())) %>% 
    mutate_if(is.factor, funs(id = as.integer(.))) 
  
  if (!missing(dynamic_treatment_maps)) {
    prepared %<>% {
      if (!is_empty(dynamic_treatment_maps$original_dyn_var)) {
        left_join(., select(dynamic_treatment_maps$mask, dynamic_treatment_maps$original_dyn_var, dynamic_treatment_id) %>% distinct(), 
                  by = dynamic_treatment_maps$original_dyn_var) 
      } else { # This means there is only one dynamic treatment
        mutate(., dynamic_treatment_id = 1) 
      }
    }
  }
  
  return(prepared)
}

to_new_dyn_ate <- function(all_ate, treatment_map) {
  ate_names <- str_replace(names(all_ate), "_(left|right)", "")
  static_ate_names <- c("private_value", "social_value", "dist.pot.group", "sms.treatment.2")
  dyn_ate_names <- c("incentive_shift", "signal_observed", "dyn_dist_pot", "reminder_info_stock")
  
  static_ate <- all_ate %>% 
    select(-one_of(str_c(rep(dyn_ate_names, each = 2), c("_left", "_right")))) %>% 
    identify_treatment_id(treatment_map) 
  
  dyn_ate <- all_ate %>% 
    select(-one_of(c(str_c(rep(c("private_value", "social_value", "sms.treatment.2"), each = 2), c("_left", "_right")), "dist.pot.group"))) %>%
    rename_(.dots = setNames(str_c(rep(intersect(dyn_ate_names, ate_names), each = 2), c("_left", "_right")), 
                             str_c(rep(intersect(static_ate_names, ate_names), each = 2), c("_left", "_right")))) %>% 
    identify_treatment_id(treatment_map) 
  
  stopifnot(nrow(static_ate) == nrow(all_ate) & nrow(static_ate) == nrow(dyn_ate))
  
  imap_dfc(lst(static = static_ate, dynamic = dyn_ate), 
           ~ select(.x, starts_with("all_treatment_id")) %>% set_colnames(str_c(.y, names(.), sep = "_")))
}

prepare_bayesian_analysis_data <- function(origin_prepared_analysis_data, 
                                           wtp_data,
                                           ...,
                                           treatment_col = c("assigned.treatment", "dist.pot.group", "sms.treatment.2", "hh.baseline.sample.pool"),
                                           prepared_treatment_maps = FALSE, 
                                           treatment_map = NULL,
                                           dynamic_treatment_map = NULL,
                                           treatment_formula = NULL,
                                           subgroup_col = c("phone_owner", "name_matched"),
                                           drop_intercept_from_dm = TRUE,
                                           nonparam_dynamics = FALSE,
                                           param_poly_order = 1,
                                           all_ate = NULL,
                                           endline_covar = c("ethnicity", "floor", "school")) {
  prep_data_arranger <- function(prep_data, ...) prep_data %>% arrange(stratum_id, new_cluster_id, name_matched, dewormed.any, ...)
  
  prepared_analysis_data <- origin_prepared_analysis_data %>% 
    mutate(new_cluster_id = factor(cluster.id) %>% as.integer(),
           age = if_else(!is.na(age), age, age.census),
           age_group = if_else(!is.na(age_group), age_group, age.census_group),
           age_squared = age^2,
           missing_covar = is.na(floor)) %>% 
    mutate(stratum = county) %>% 
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
                             function(col, col_mean) sqrt(sum(map_count * ((col - col_mean)^2))/(sum(map_count) - 1)) * 2) %>%  
              if_else(is_factor_col, 1, .),
            center = map_means) # Center and scale to have SD = 0.5
  }
  
  remove_dup_treatment_dm_cols <- TRUE
  
  is_dynamic_model <- !is.null(dynamic_treatment_map)
  
  if (is_dynamic_model) {
    if (prepared_treatment_maps) {
      dyn_treat_maps <- dynamic_treatment_map
    } else {
      dyn_treat_maps <- prepare_dynamic_treatment_maps(dynamic_treatment_map, prepared_analysis_data)
    }
    
    dynamic_treatment_map <- dyn_treat_maps$map
    dynamic_treatment_mask_map <- dyn_treat_maps$mask
    dyn_formula_var <- dyn_treat_maps$dyn_formula_var
  } else {
    dyn_formula_var <- NULL
    dynamic_treatment_mask_map <- NULL
    dynamic_treatment_dm <- NULL
  }
 
  if (is.null(treatment_map)) { 
    treatment_map <- expand_(prepared_analysis_data, c(treatment_col, subgroup_col)) %>% {
        if ("phone_owner" %in% names(.)) arrange(., phone_owner) else return(.) 
      } 
  } else {
    treatment_col <- intersect(names(treatment_map), names(prepared_analysis_data))
  }
  
  if (!prepared_treatment_maps) {
    if (is_dynamic_model) treatment_map %<>% prepare_treatment_map(dyn_treat_maps) else treatment_map %<>% prepare_treatment_map()
  }
  
  census_covar_map <- count_(prepared_analysis_data, c("age", "gender")) %>% 
    mutate(census_covar_id = seq_len(n())) 
    
  prepared_analysis_data %<>% 
    left_join(treatment_map, unique(c(treatment_col, dyn_formula_var, subgroup_col))) %>% 
    left_join(select(census_covar_map, -n), c("age", "gender")) %T>% {
      if (any(is.na(.$all_treatment_id))) warning(sprintf("%d observations with NA all treatment ID", sum(is.na(.$all_treatment_id))))
    } %>%
    # filter(!is.na(all_treatment_id)) %>% 
    prep_data_arranger() %>% 
    mutate(obs_index = seq_len(n()))
  
 
  if (is.null(treatment_formula)) { 
    treatment_formula <- str_c("~ ", str_c(c(treatment_col, subgroup_col), collapse = " * ")) %>% {
        if (!is.null(subgroup_col)) str_c(., " - ", str_c(subgroup_col, collapse = " - ")) else return(.)
      } %>% 
      as.formula()
  }
  
  # Find redundant columns, excluding the intercept if we are required to keep it
  detect_redund_col <- . %>% { map_lgl(., ~ n_distinct(.) > 1) | (str_detect(names(.), fixed("intercept")) & !drop_intercept_from_dm) }
  
  treatment_map_design_matrix <- treatment_map %>%  
    model_matrix(treatment_formula) %>% 
    rename(intercept = `(Intercept)`) %>% 
    magrittr::extract(, detect_redund_col(.)) %>% {
      if (!remove_dup_treatment_dm_cols) return(.) else magrittr::extract(., , !duplicated(t(.))) # Remove redundant columns (rows?)
    }
  
  if (!is_null(subgroup_col)) {
    subgroup_treatment_map_dm <- treatment_map %>% 
      distinct_(.dots = subgroup_col) %>% 
      mutate(subgroup_id = seq_len(n())) %>% 
      left_join(treatment_map, ., subgroup_col) %>% 
      select(subgroup_id) %>% 
      bind_cols(treatment_map_design_matrix, .) %>% 
      group_by(subgroup_id) %>% 
      summarize_all(~ 1 * (sum(.) > 1)) %>% 
      ungroup() 
    
    omitted_subgroup_treatment_col <- subgroup_treatment_map_dm %>% 
      select(-subgroup_id) %>% 
      select_if(~ sum(.) > 1) %>% 
      names()
    
    omitted_subgroup_id <- subgroup_treatment_map_dm %>% 
      select(-one_of(omitted_subgroup_treatment_col)) %>% 
      filter_at(vars(-subgroup_id), all_vars(. == 0)) %>% 
      pull(subgroup_id)
   
    omitted_subgroup <- tibble(subgroup_treatment_col = omitted_subgroup_treatment_col) %>% 
      mutate(subgroup_id = omitted_subgroup_id)
    
    subgroup_treatment_col_map <- subgroup_treatment_map_dm %>% 
      select(-one_of(omitted_subgroup_treatment_col)) %>% 
      anti_join(omitted_subgroup, "subgroup_id") %>%
      gather(subgroup_treatment_col, col_mask, -subgroup_id) %>% 
      group_by(subgroup_id) %>% 
      filter(col_mask == 1) %>% 
      ungroup() %>% 
      select(-col_mask) %>% 
      bind_rows(omitted_subgroup) %>% 
      arrange(subgroup_id)
    
    treatment_map_design_matrix %<>% select(subgroup_treatment_col_map$subgroup_treatment_col)
    
    subgroup_treatment_col_sizes <- subgroup_treatment_col_map %>% 
      count(subgroup_id) %>% 
      pull(n)
  } else {
    subgroup_treatment_col_sizes <- ncol(treatment_map_design_matrix)
  }
  
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
  
  name_matched_data <- filter(prepared_analysis_data, name_matched == 1) %>% 
    prep_data_arranger() # Should already be in the right order, but to be on the safe size
  
  monitored_data <- filter(prepared_analysis_data, name_matched == 0) %>% 
    prep_data_arranger() # Should already be in the right order, but to be on the safe size
  
  matching_error_data <- prepared_analysis_data %>% 
    filter(sms.treatment.2 == "sms.control",
           true.monitored | !monitored) # using only true monitored and those not monitored because they weren't in a study sample
     
  stratum_map <- distinct(prepared_analysis_data, stratum_id, stratum)
 
  if (!is.null(all_ate)) {
    num_ate_pairs <- nrow(all_ate)
    
    observed_takeup_total <- prepared_analysis_data %>% 
      group_by(all_treatment_id) %>% 
      summarize(takeup_total = sum(dewormed.any)) %>% 
      ungroup() 
    
    if (nonparam_dynamics) {
      combined_ate <- to_new_dyn_ate(all_ate, treatment_map)
      
      ate_treatments <- combined_ate %>% 
        get_unique_treatments()
      
      if (!is_empty(subgroup_col)) {
        ate_treatments %<>% 
          left_join(select(treatment_map, all_treatment_id, subgroup_col), c("static_all_treatment_id" = "all_treatment_id")) 
      }
      
      missing_treatment <- ate_treatments %>%
        plyr::alply(1, function(treat_row) {
          (if (!is_empty(subgroup_col)) {
            semi_join(prepared_analysis_data, select(treat_row, subgroup_col), subgroup_col) 
          } else { prepared_analysis_data }) %>% 
            filter(all_treatment_id != treat_row$static_all_treatment_id | all_treatment_id != treat_row$dynamic_all_treatment_id) %>% 
            pull(obs_index)
        })
      
      observed_stratum_treatment <- ate_treatments %>%
        plyr::adply(1, function(treat_row) {
          (if (!is_empty(subgroup_col)) {
            semi_join(prepared_analysis_data, select(treat_row, subgroup_col), subgroup_col) 
          } else { prepared_analysis_data }) %>% 
            filter(all_treatment_id == treat_row$static_all_treatment_id & all_treatment_id == treat_row$dynamic_all_treatment_id) %>% 
            select(stratum, obs_index)
        }) 
      
      observed_treatment <- observed_stratum_treatment %>% 
        group_by_at(vars(ends_with("all_treatment_id"))) %>% 
        do(treatment_ids = .$obs_index) %>% 
        ungroup() %>%
        right_join(select(ate_treatments, static_all_treatment_id, dynamic_all_treatment_id), 
                   str_c(c("static", "dynamic"), "_all_treatment_id")) %>% 
        pull(treatment_ids)
      
      ate_pairs <- combined_ate %>% 
        select(matches("(static|dynamic)_all_treatment_id_(left|right)")) %>% 
        left_join(ate_treatments, c("static_all_treatment_id_left" = "static_all_treatment_id",
                                    "dynamic_all_treatment_id_left" = "dynamic_all_treatment_id")) %>% 
        left_join(ate_treatments, c("static_all_treatment_id_right" = "static_all_treatment_id",
                                    "dynamic_all_treatment_id_right" = "dynamic_all_treatment_id"), suffix = c("_left", "_right")) %>% 
        select(rank_id_left, rank_id_right)
      
      ate_treatments %<>% 
        left_join(observed_takeup_total, c("static_all_treatment_id" = "all_treatment_id")) %>% 
        mutate(takeup_total = coalesce(takeup_total, 0L))
      
      observed_stratum_treatment %<>% 
        rename(all_treatment_id = static_all_treatment_id)
    } else {
      all_ate %<>% identify_treatment_id(treatment_map) 
      
      stopifnot(nrow(all_ate) == num_ate_pairs)
      
      ate_treatments <- all_ate %>% get_unique_treatments()
      
      if (!is_empty(subgroup_col)) {
        ate_treatments %<>% left_join(select(treatment_map, all_treatment_id, subgroup_col), "all_treatment_id")
      }
      
      missing_treatment <- ate_treatments %>%
        plyr::alply(1, function(treat_row) {
          (if (!is_empty(subgroup_col)) {
            semi_join(prepared_analysis_data, select(treat_row, subgroup_col), subgroup_col) 
          } else { prepared_analysis_data }) %>% 
            filter(all_treatment_id != treat_row$all_treatment_id) %>% 
            pull(obs_index)
        })
      
      observed_stratum_treatment <- ate_treatments %>%
        plyr::adply(1, function(treat_row) {
          (if (!is_empty(subgroup_col)) {
            semi_join(prepared_analysis_data, select(treat_row, subgroup_col), subgroup_col) 
          } else { prepared_analysis_data }) %>% 
            filter(all_treatment_id == treat_row$all_treatment_id) %>% 
            select(stratum, obs_index)
        })
      
      observed_treatment <- observed_stratum_treatment %>% 
        group_by(all_treatment_id) %>% 
        do(treatment_ids = .$obs_index) %>% 
        ungroup() %>%
        right_join(select(ate_treatments, all_treatment_id), "all_treatment_id") %>% 
        pull(treatment_ids)
      
      ate_pairs <- all_ate %>% 
        select(all_treatment_id_left, all_treatment_id_right) %>%
        left_join(ate_treatments, c("all_treatment_id_left" = "all_treatment_id")) %>% 
        left_join(ate_treatments, c("all_treatment_id_right" = "all_treatment_id"), suffix = c("_left", "_right")) %>% 
        select(rank_id_left, rank_id_right)
      
      ate_treatments %<>% 
        left_join(observed_takeup_total, "all_treatment_id") %>% 
        mutate(takeup_total = coalesce(takeup_total, 0L))
    }
    
    observed_stratum_takeup_total <- prepared_analysis_data %>% 
      group_by(stratum, all_treatment_id) %>% 
      summarize(takeup_total = sum(dewormed.any)) %>% 
      ungroup() 
  } else {
    ate_treatments <- NULL
    
    missing_treatment <- NULL 
    observed_treatment <- NULL
    
    num_ate_pairs <- 0
    
    ate_pairs <- NULL
    
    observed_stratum_takeup_total <- NULL
    observed_stratum_treatment <- NULL
  }
 
  num_deworming_days <- 12L 
  
  hazard_day_map <- diag(c(rep(1, num_deworming_days), 0))[, 1:num_deworming_days]
  hazard_day_triangle_map <- lower.tri(diag(num_deworming_days + 1L)[, seq_len(num_deworming_days)]) * 1
  relevant_latent_var_map <- hazard_day_map + hazard_day_triangle_map
  
  dewormed_day_any <- prepared_analysis_data %>% 
    pull(dewormed.day.any) %>% 
    coalesce(num_deworming_days + 1L)
  
  num_relevant_obs_days <- sum(relevant_latent_var_map[dewormed_day_any, ])
  
  obs_relevant_latent_var <- relevant_latent_var_map[dewormed_day_any, ]
  obs_relevant_days <- rowSums(obs_relevant_latent_var)
  
  census_covar_id <- prepared_analysis_data$census_covar_id
  census_covar_dm <- census_covar_map_dm[census_covar_id, ]
 
  census_covar_dm_long <- census_covar_dm %>% magrittr::extract(rep(seq_len(nrow(.)), obs_relevant_days), ) 
  
  # treatment_map_design_matrix %<>%
  #   magrittr::extract(., , magrittr::extract(., c(unique(prepared_analysis_data$all_treatment_id), ate_treatments$all_treatment_id), ) %>% 
  #                       detect_redund_col())
  
  num_all_treatments <- nrow(treatment_map_design_matrix)
  num_all_treatment_coef <- ncol(treatment_map_design_matrix) 

  obs_treatment <- prepared_analysis_data$all_treatment_id 
  treatment_design_matrix <- treatment_map_design_matrix[obs_treatment, ]
   
  treatment_design_matrix_long <- treatment_map_design_matrix %>% magrittr::extract(rep(obs_treatment, obs_relevant_days), )
  treatment_design_matrix_qr <- qr(treatment_design_matrix_long)
  Q_treatment_design_matrix_long <- qr.Q(treatment_design_matrix_qr) * sqrt(num_relevant_obs_days - 1)
  R_treatment_design_matrix_long <- qr.R(treatment_design_matrix_qr) / sqrt(num_relevant_obs_days - 1)
  R_inv_treatment_design_matrix_long <- solve(R_treatment_design_matrix_long)
  
  cluster_id_long <- rep(prepared_analysis_data$new_cluster_id, obs_relevant_days)
  
  relevant_day_dewormed <- obs_relevant_days %>% # pmin(dewormed_day_any, num_deworming_days) %>% 
    accumulate(add) %>% 
    magrittr::extract(dewormed_day_any <= num_deworming_days)
  
  relevant_dewormed_any_daily <- rep.int(0, num_relevant_obs_days) %>% 
    magrittr::inset(relevant_day_dewormed, 1)
  
  dewormed_day_long <- map(obs_relevant_days, seq_len) %>% unlist()
  
  num_param_dyn_coef <- (num_all_treatment_coef - 1) * param_poly_order
  
  days_poly_trend <- ((1:num_deworming_days) - 1) %>% 
    scale() %>% 
    matrix(num_deworming_days, param_poly_order) %>% 
    plyr::aaply(1, accumulate, multiply_by) %>% 
    magrittr::extract(, rep(seq_len(ncol(.)), each = num_all_treatment_coef - 1))  
  
  param_dyn_treatment_map <- treatment_map_design_matrix[, rep(2:num_all_treatment_coef, param_poly_order)] %>% 
    plyr::alply(1, function(mask_row) matrix(as.numeric(mask_row), num_deworming_days, length(mask_row), byrow = TRUE) * days_poly_trend)
 
  param_dyn_treatment_design_matrix_long <- obs_relevant_days %>% 
    plyr::llply(function(curr_relevant_days) days_poly_trend[1:curr_relevant_days, ]) %>% 
    do.call(rbind, .) %>% 
    multiply_by(treatment_design_matrix_long[, rep(2:num_all_treatment_coef, param_poly_order)])
  
  param_dyn_qr <- qr(param_dyn_treatment_design_matrix_long)
  
  Q_param_dyn_treatment_design_matrix_long <- qr.Q(param_dyn_qr) * sqrt(num_relevant_obs_days - 1)
  R_param_dyn_treatment_design_matrix_long <- qr.R(param_dyn_qr) / sqrt(num_relevant_obs_days - 1)
  R_inv_param_dyn_treatment_design_matrix_long <- solve(R_param_dyn_treatment_design_matrix_long)
  
  private_value_calendar_coef <- treatment_map_design_matrix %>% 
    names() %>% 
    str_detect("private_valuecalendar") %>% 
    which()
  
  private_value_bracelet_coef <- treatment_map_design_matrix %>% 
    names() %>% 
    str_detect("private_valuebracelet") %>% 
    which()
 
  stan_data_list <- lst(
    # Save meta data 
    
    prepared_analysis_data,
     
    treatment_map,
    
    all_ate,
    
    missing_treatment,
    observed_treatment,
    
    census_covar_map, 
    stratum_map,
    cluster_map = distinct(prepared_analysis_data, new_cluster_id, cluster.id),
    
    observed_stratum_takeup_total,
    observed_stratum_takeup_prop = if (!is_null(observed_stratum_takeup_total)) observed_stratum_takeup_total %>% 
      left_join(count(observed_stratum_treatment, all_treatment_id, stratum), c("all_treatment_id", "stratum")) %>% 
      transmute(all_treatment_id, stratum, takeup_prop = takeup_total / n),
    observed_stratum_treatment,
    
    # Data passed to model
    
    param_poly_order,
    
    num_deworming_days, 
    
    treatment_map_design_matrix,
    treatment_design_matrix,
    treatment_design_matrix_long,
    Q_treatment_design_matrix_long,
    R_treatment_design_matrix_long,
    R_inv_treatment_design_matrix_long,
    num_all_treatments,
    num_all_treatment_coef,
    num_subgroups = length(subgroup_treatment_col_sizes),
    subgroup_treatment_col_sizes,
    
    num_dynamic_treatments = if (is_empty(dynamic_treatment_map)) 0 else n_distinct(dynamic_treatment_map$dynamic_treatment_id),
    num_dynamic_treatment_col = if (is_empty(dynamic_treatment_map)) 0 else dynamic_treatment_map %>% ncol() %>% subtract(1),
    all_treatment_dyn_id = treatment_map$dynamic_treatment_id,
    signal_observed_coef = dynamic_treatment_map %>% colnames() %>% str_which("signal_observed"),
    reminder_info_coef = dynamic_treatment_map %>% colnames() %>% str_which("reminder_info"),
    num_signal_observed_coef = length(signal_observed_coef),
    num_reminder_info_coef = length(reminder_info_coef), 
    
    dynamic_treatment_map = if (!is_empty(dynamic_treatment_map)) {
      dynamic_treatment_map %>% dlply(.(dynamic_treatment_id), . %>% select(-dynamic_treatment_id) %>% as.matrix()) 
    },
    
    param_dyn_treatment_map,
    param_dyn_treatment_design_matrix_long,
    Q_param_dyn_treatment_design_matrix_long,
    R_param_dyn_treatment_design_matrix_long,
    R_inv_param_dyn_treatment_design_matrix_long,
    
    num_private_value_calendar_coef = length(private_value_calendar_coef),
    num_private_value_bracelet_coef = length(private_value_bracelet_coef),
    private_value_calendar_coef,
    private_value_bracelet_coef,
    private_value_bracelet_indicator = colnames(treatment_map_design_matrix) %>% str_detect("^private_valuebracelet$") %>% which,
    not_private_value_bracelet_coef = seq_len(num_all_treatment_coef) %>% setdiff(private_value_bracelet_coef),
    
    census_covar_map_dm,
    num_census_covar_coef = ncol(census_covar_map_dm),
    num_distinct_census_covar = nrow(census_covar_map_dm),
    census_covar_id,
    
    census_covar_dm,
    
    census_covar_dm_long,
    
    endline_covar_dm,
    num_endline_covar_coef = ncol(endline_covar_dm),
    
    stratum_covar_id = prepared_analysis_data %>% prep_data_arranger(missing_covar, obs_index) %$% obs_index, # First obs then missing
    stratum_missing_covar_sizes = prepared_analysis_data %>% 
      filter(missing_covar | (select_(., .dots = endline_covar) %>% map(is.na) %>% reduce(or))) %>% 
      count(stratum_id) %>% 
      arrange(stratum_id) %$% 
      n,
    
    obs_treatment,
  
    treatment_id = prepared_analysis_data %>% arrange(all_treatment_id, obs_index) %$% obs_index,
    dynamic_treatment_id = if (is_empty(dynamic_treatment_mask_map)) rep(0, nrow(prepared_analysis_data)) else prepared_analysis_data$dynamic_treatment_id,
    
    num_obs = length(obs_treatment), 
    
    treatment_sizes = count(prepared_analysis_data, all_treatment_id) %>% arrange(all_treatment_id) %>% pull(n),
    
    dewormed_any = prepared_analysis_data$dewormed.any,
    
    dewormed_day_any, # Set dewormed day to 13 for the non-dewormed
    
    relevant_dewormed_any_daily,
    
    dewormed_day_long,
    
    num_dewormed = sum(dewormed_any),
    dewormed_ids = prepared_analysis_data %>% filter(dewormed.any) %>% arrange(stratum_id) %>% pull(obs_index),
    
    cluster_id = prepared_analysis_data$new_cluster_id,
    cluster_id_long,
    stratum_id = prepared_analysis_data$stratum_id,
    num_clusters = n_distinct(cluster_id),
    
    num_strata = n_distinct(stratum_id),
    strata_sizes = count(prepared_analysis_data, stratum_id) %>% arrange(stratum_id) %>% pull(n),
    cluster_sizes = count(prepared_analysis_data, new_cluster_id) %>% arrange(new_cluster_id) %>% pull(n),
    strata_num_clusters = distinct(prepared_analysis_data, stratum_id, new_cluster_id) %>% 
      count(stratum_id) %>% 
      arrange(stratum_id) %>% 
      pull(n),
    strata_cluster_ids = distinct(prepared_analysis_data, stratum_id, new_cluster_id) %>% 
      arrange(stratum_id) %>% 
      pull(new_cluster_id),
    cluster_obs_ids = prepared_analysis_data %>% 
      arrange(new_cluster_id) %>% 
      pull(obs_index),
    cluster_stratum_ids = prepared_analysis_data %>% 
      distinct(new_cluster_id, stratum_id) %>% 
      arrange(new_cluster_id) %>% 
      pull(stratum_id),
    
    strata_dewormed_sizes = prepared_analysis_data %>% 
      filter(dewormed.any) %>% 
      count(stratum_id) %>% 
      arrange(stratum_id) %>% 
      pull(n),
    
    stratum_dewormed_index = strata_sizes[-num_strata] %>% 
      accumulate(add, .init = 0) %>% 
      map2(strata_dewormed_sizes, ~ rep(.x, each = .y)) %>% 
      unlist() %>% 
      subtract(dewormed_ids, .),
    
    hazard_day_map,
    hazard_day_triangle_map,
    relevant_latent_var_map,
    
    num_relevant_obs_days,
    num_param_dyn_coef,
    
    # ATE
    
    num_ate_treatments = if (!is.null(ate_treatments)) nrow(ate_treatments) else 0,
    
    observed_takeup_total = ate_treatments$takeup_total,
    ate_treatments = if (!is.null(ate_treatments)) { 
        if (nonparam_dynamics) {
          ate_treatments %>% select(static_all_treatment_id, dynamic_all_treatment_id)
        } else {
          ate_treatments$all_treatment_id 
        } 
      }else 1,
    
    missing_obs_ids = if (num_ate_treatments > 0) unlist(missing_treatment) else 1,
    num_missing_obs_ids = if (num_ate_treatments > 0) length(missing_obs_ids) else 0,
    missing_treatment_sizes = if (!is.null(all_ate)) map_int(missing_treatment, length) else 1,
    observed_obs_ids = if (num_ate_treatments > 0) unlist(observed_treatment) else 1,
    num_observed_obs_ids = if (num_ate_treatments > 0) length(observed_obs_ids) else 0,
    observed_treatment_sizes = if (!is.null(all_ate)) map_int(observed_treatment, length) else 1,
    
    missing_treatment_stratum_id = stratum_id[missing_obs_ids], 
    missing_treatment_cluster_id = cluster_id[missing_obs_ids],
    
    missing_census_covar_dm = census_covar_dm[missing_obs_ids, ],
    
    num_ate_pairs = if (num_ate_treatments > 0) num_ate_pairs else 0,
    ate_pairs = if (num_ate_pairs > 0) ate_pairs else c(1, 1),
    
    observed_dewormed_day = dewormed_day_any[observed_obs_ids],
    ...
  ) %>% 
    modifyList(prepare_bayes_wtp_data(origin_prepared_analysis_data, wtp_data, stratum_map))
  
  # stan_data_list$dynamic_initializer <- function() { 
  #   num_not_private_value_bracelet_coef <- with(stan_data_list, num_all_treatment_coef - num_private_value_bracelet_coef)
  #   
  #   lst(
  #     hyper_beta_raw = rep(0, num_not_private_value_bracelet_coef),
  #     stratum_beta_raw = array(rep(0, num_not_private_value_bracelet_coef * stan_data_list$num_strata),
  #                              dim = c(stan_data_list$num_strata, num_not_private_value_bracelet_coef)),
  #     hyper_dynamic_treatment_coef_raw = rep(0, stan_data_list$num_dynamic_treatment_col),
  #     stratum_dynamic_treatment_coef_raw = with(stan_data_list, array(rep(0, num_dynamic_treatment_col * num_strata),
  #                                                                     dim = c(num_strata, num_dynamic_treatment_col)))
  #   )
  # }
  
  return(stan_data_list)
} 

estimate_treatment_deworm_prob_dm <- function(applicable_obs,
                                              static_dm_row, 
                                              dyn_dm_row,
                                              sim_stratum_hazards,
                                              kappa_census_covar_data = NULL,
                                              treatment_summarizer = function(d) d,
                                              days_treated = 0:11,
                                              by_id = "obs_index") {
  beta_params <- applicable_obs %>% 
    select(matches("^stratum_beta_\\d+")) %>% 
    magrittr::set_names(static_dm_row %>% select(-matches("private_valuebracelet")) %>% names())
  
  lower_mask_mat <- matrix(nrow = 12, ncol = 12) %>% 
    upper.tri(diag = TRUE) %>%
    magrittr::extract(, days_treated + 1) %>% 
    matrix(nrow = 12 * nrow(applicable_obs), ncol = NCOL(.), byrow = TRUE)
  
  upper_mask_mat <- 2 * (1 - lower_mask_mat) / 2 
  
  num_days_treated <- length(days_treated)
  
  calculate_days_treated_matrix <- function(log_kappa) {
    upper <- matrix(log_kappa, nrow = length(log_kappa), ncol = num_days_treated) * upper_mask_mat
    lower <- (t(plyr::aaply(matrix(rep(log_kappa, num_days_treated), nrow = length(log_kappa)/num_days_treated, ncol = num_days_treated, byrow = TRUE), 
                      2, rep, each = num_days_treated, .drop = FALSE)) * lower_mask_mat) 
    
    return(upper + lower)
  }
 
  log_prod_alpha_to_prob_deworm <- . %>%  
    exp() %>%
    subtract(cbind(1, .[, -12]), .) %>% 
    set_colnames(str_c("prob_deworm_deworming_day_", 1:12))
  
  applicable_obs %>% 
    mutate(log_kappa_static_treatment = 
             c((as.matrix(beta_params) %*% unlist(select(static_dm_row, -matches("private_valuebracelet")))) + 
               (as.matrix(select(beta_params, matches("private_valuecalendar"))) %*% unlist(select(static_dm_row, matches("private_valuebracelet")))))) %>% 
    left_join(sim_stratum_hazards, c("iter_id", "stratum_id")) %>% 
    mutate(cluster_hazard = stratum_hazard * cluster_hazard_effect,
           log_kappa_dyn_treat = c(rowSums(dyn_dm_row[deworming_day, ] * select(., matches("^stratum_dynamic_treatment_coef_\\d+"))))) %>%
    select(-c(cluster_hazard_effect, stratum_hazard),
           -starts_with("stratum_beta"),
           -starts_with("stratum_dynamic_treatment_coef_"),
           -starts_with("dyn_treatment_"),
           -one_of(names(beta_params)),
           -contains("private_valuebracelet")) %>% 
    arrange(subgroup_id, deworming_day) %>% 
    bind_cols(calculate_days_treated_matrix(.$log_kappa_dyn_treat) %>% 
                as_tibble() %>% 
                magrittr::set_names(str_c("log_kappa_dyn_treat_days_", days_treated))) %>% 
    mutate_at(vars(num_range("log_kappa_dyn_treat_days_", days_treated)), 
              funs(log_alpha = - exp(log_kappa_census_covar + log_kappa_static_treatment + .) * cluster_hazard)) %>% 
    {
      if (num_days_treated > 1) {
        gather(., log_key, log_val, matches("^log_kappa_dyn_treat_days_\\d+(_log_alpha)?$")) %>% 
          mutate(dyn_treat_days = as.integer(str_extract(log_key, "\\d+")),
                 log_key = str_replace_all(log_key, c(".+(log_alpha)$" = "\\1", "(.+)(_days_\\d+)$" = "\\1"))) %>% 
          select(-log_kappa_dyn_treat) %>% 
          spread(log_key, log_val) 
      } else mutate(., dyn_treat_days = days_treated)
    } %>% 
    select_at(vars(-cluster_hazard, 
                   -starts_with("log_kappa_static"), 
                   -starts_with("log_kappa_census"), 
                   -starts_with("log_kappa_dyn_treat_days"),
                   -starts_with("stratum_census_covar_coef"))) %>% 
    { 
      inner_join(select(., -log_alpha) %>% 
                   mutate(deworming_day = sprintf("log_kappa_dyn_treat_deworming_day_%02d", deworming_day)) %>% 
                   spread(deworming_day, log_kappa_dyn_treat), 
                 select(., -log_kappa_dyn_treat) %>% 
                   mutate(deworming_day = sprintf("log_alpha_deworming_day_%02d", deworming_day)) %>% 
                   spread(deworming_day, log_alpha),
                 by = setdiff(names(.), c("deworming_day", "log_alpha", "log_kappa_dyn_treat"))) 
    } %>% 
    cbind(select(., log_alpha_deworming_day_01:log_alpha_deworming_day_12) %>% 
            accumulate(`+`) %>% 
            set_names(str_c("log_prod_alpha_deworming_day_", 1:12))) %>% 
    { 
      if (!is.null(kappa_census_covar_data)) {
        inner_join(., kappa_census_covar_data, "cluster_id") %>% 
          cbind(multiply_by(select(., starts_with("log_prod_alpha_deworming_day")), .$kappa_census_covar) %>% 
                  log_prod_alpha_to_prob_deworm()) %>% 
          select(-kappa_census_covar)
      } else {
        cbind(., 
              select(., starts_with("log_prod_alpha_deworming_day")) %>% 
                log_prod_alpha_to_prob_deworm())
      }
    } %>% 
    select(-starts_with("log_prod_alpha")) %>% 
    gather(deworming_day_key, deworming_day_val, 
           log_alpha_deworming_day_01:log_alpha_deworming_day_12, 
           log_kappa_dyn_treat_deworming_day_01:log_kappa_dyn_treat_deworming_day_12, 
           prob_deworm_deworming_day_1:prob_deworm_deworming_day_12) %>% 
    separate(deworming_day_key, c("deworming_day_key", "deworming_day"), sep = "_(?=\\d+)") %>%  
    mutate(deworming_day = as.integer(deworming_day)) %>% 
    spread(deworming_day_key, deworming_day_val) %>% 
    magrittr::set_names(str_replace(names(.), "_deworming_day", "")) %>% 
    treatment_summarizer() 
}

estimate_treatment_deworm_prob <- function(dm_datarow, 
                                           iter_est, 
                                           sim_stratum_hazards,
                                           full_dyn_treatment_dm, 
                                           kappa_census_covar_data = NULL,
                                           days_treated = 0:11,
                                           treatment_summarizer = function(obs_estimates) {
                                             obs_estimates %>% 
                                               group_by(deworming_day, dyn_treat_days) %>% 
                                               summarize_at(vars(starts_with("prob_deworm"), starts_with("log_alpha")), mean) %>% 
                                               ungroup() %>% 
                                               left_join(dplyr::count(obs_estimates, deworming_day, dyn_treat_days), by = c("deworming_day", "dyn_treat_days")) %>% 
                                               dplyr::rename(day_size = n) %>% 
                                               bind_rows(group_by(., dyn_treat_days) %>% 
                                                           summarize(deworming_day = 13, prob_deworm = 1 - sum(prob_deworm))) 
                                           },
                                           subgroup_col = c("phone_owner", "name_matched"),
                                           by_id = "obs_index") {
  static_dm_row <- select(dm_datarow, -all_treatment_id, -one_of(subgroup_col), -starts_with("dyn_treatment_"))
  dyn_dm_row <- select(dm_datarow, starts_with("dyn_treatment_")) %>% 
    unlist() %>% 
    matrix(nrow = nrow(full_dyn_treatment_dm), ncol = ncol(full_dyn_treatment_dm), byrow = TRUE) %>% 
    magrittr::multiply_by(full_dyn_treatment_dm)
  
  left_join(dm_datarow, iter_est, subgroup_col) %>% 
    estimate_treatment_deworm_prob_dm(static_dm_row, dyn_dm_row, sim_stratum_hazards, 
                                      days_treated = days_treated, 
                                      treatment_summarizer = treatment_summarizer,
                                      kappa_census_covar_data = kappa_census_covar_data,
                                      by_id = by_id) 
}

fast_estimate_deworm_prob <- function(iter_cluster_parameters,
                                      sim_stratum_hazards,
                                      treatments_info,
                                      treatment_map,
                                      treatment_map_dm,
                                      all_treat_dyn_id,
                                      dyn_treatment_mask_map,
                                      full_dyn_treatment_map_dm,
                                      average_over = NULL,
                                      days_treated = c(0, 11),
                                      subgroup_col = "phone_owner") {
  dyn_treatment_mask_map %<>% set_names(stringr::str_c("dyn_treatment_", names(.)))
  
  fast_est_group_prob <- function(treatment_group) {
    treatment_group %>%
      distinct_(.dots = c("all_treatment_id", "dynamic_treatment_id"), .keep_all = TRUE) %>% 
      select(all_treatment_id, dynamic_treatment_id, one_of(subgroup_col)) %>% 
      bind_cols(treatment_map_dm[.$all_treatment_id, ], slice(dyn_treatment_mask_map, .$dynamic_treatment_id)) %>% 
      select(-dynamic_treatment_id) %>% 
      plyr::ddply(.(all_treatment_id), estimate_treatment_deworm_prob, 
                  iter_est = iter_cluster_parameters %>% 
                    mutate(subgroup_id = seq_len(nrow(.))), # For some reason n() doesn't seem to work in running in parallel (tried dplyr::n() as well)
                  sim_stratum_hazards = sim_stratum_hazards,
                  full_dyn_treatment_dm = full_dyn_treatment_map_dm,
                  days_treated = days_treated,
                  treatment_summarizer = function(d) d,
                  by_id = "subgroup_id",
                  subgroup_col = subgroup_col) %>% 
      mutate(alpha = exp(log_alpha), 
             kappa_dyn_treat = exp(log_kappa_dyn_treat)) %>%
      dplyr::group_by(deworming_day, dyn_treat_days) %>% 
      mutate(treatment_group_size = sum(subgroup_size)) %>% 
      dplyr::group_by(deworming_day, dyn_treat_days, treatment_group_size) %>% 
      summarize_at(vars(prob_deworm, alpha, kappa_dyn_treat), funs(weighted.mean(., subgroup_size))) %>% 
      ungroup() 
  }
  
  treatments_to_eval <- identify_treatment_id(treatments_info, treatment_map) 
  
  treatment_col <- c("phone_owner", "dist.pot.group", stringr::str_subset(names(treatments_to_eval), "(?<!_id)$")) %>% 
    setdiff(average_over)
  
  plyr::ddply(treatments_to_eval, treatment_col, fast_est_group_prob) 
}

estimate_deworm_prob <- function(iter_cluster_parameters,
                                 obs_meta_data, 
                                 sim_stratum_hazards,
                                 census_covar_dm, 
                                 treatments_info,
                                 treatment_map,
                                 treatment_map_dm,
                                 dyn_treatment_mask_map,
                                 full_dyn_treatment_map_dm,
                                 days_treated = 0:11) {
  tryCatch({
  kappa_census_covar_data <- iter_cluster_parameters %>% 
    left_join(obs_meta_data, c("stratum_id", "cluster_id", "phone_owner", "name_matched")) %>% 
    transmute(obs_index, cluster_id,
              kappa_census_covar = select(., starts_with("stratum_census_covar_coef")) %>% 
                magrittr::multiply_by(census_covar_dm) %>% 
                rowSums() %>% 
                exp()) 
  
  dyn_treatment_mask_map %<>% set_names(stringr::str_c("dyn_treatment_", names(.)))
  
  est_group_prob <- . %>% 
      distinct_(.dots = c("all_treatment_id", "dynamic_treatment_id"), .keep_all = TRUE) %>% 
      select(all_treatment_id, dynamic_treatment_id, phone_owner, name_matched) %>% 
      bind_cols(treatment_map_dm[.$all_treatment_id, ], slice(dyn_treatment_mask_map, .$dynamic_treatment_id)) %>% 
      select(-dynamic_treatment_id) %>% 
      plyr::ddply(.(all_treatment_id), estimate_treatment_deworm_prob, 
                  iter_est = iter_cluster_parameters, 
                  sim_stratum_hazards = sim_stratum_hazards,
                  kappa_census_covar_data = kappa_census_covar_data,
                  full_dyn_treatment_dm = full_dyn_treatment_map_dm,
                  days_treated = days_treated,
                  # treatment_summarizer = function(d) d,
                  by_id = "obs_index") %>% 
      group_by(deworming_day, dyn_treat_days) %>% 
      summarize_at(vars(prob_deworm, log_alpha), funs(weighted.mean(., day_size))) %>%
      ungroup()
  
  identify_treatment_id(treatments_info, treatment_map) %>% 
    plyr::ddply(setdiff(c("phone_owner", "dist.pot.group", stringr::str_subset(names(.), "(?<!_id)$")), "name_matched"), 
                est_group_prob,
                .progress = "text") 
  }, error = function(err) sprintf("Failed in iteration %d: %s", iter_cluster_parameters$iter_id[1], err))
}

estimate_deworm_prob_ate <- function(iter_parameters,
                                     obs_meta_data, 
                                     sim_stratum_hazards,
                                     census_covar_dm, 
                                     ate_info,
                                     treatment_map,
                                     treatment_map_dm,
                                     all_treat_dyn_id,
                                     dyn_treatment_mask_map,
                                     full_dyn_treatment_map_dm,
                                     inner_parallel = FALSE) {
  est_group_ate <- function(ate_group, iter_est, sim_stratum_hazards, treatment_map_dm, dyn_treatment_mask_map, full_dyn_treatment_map_dm,
                            group_suffix = "right", other_group_suffix = "left") {
    group_id_colnames <- stringr::str_c(c("all_treatment_id_", "dynamic_treatment_id_"), group_suffix)
    obs_est <- ate_group %>%
      distinct_(.dots = group_id_colnames, .keep_all = TRUE) %>% 
      rename_(.dots = group_id_colnames %>% setNames(stringr::str_replace(., stringr::str_interp("_${group_suffix}$"), ""))) %>% 
      select(all_treatment_id, dynamic_treatment_id, phone_owner, name_matched) %>% 
      bind_cols(treatment_map_dm[.$all_treatment_id, ], slice(dyn_treatment_mask_map, .$dynamic_treatment_id)) %>% 
      select(-dynamic_treatment_id) %>% 
      plyr::ddply(.(all_treatment_id), estimate_treatment_deworm_prob, 
                  iter_est = iter_est, 
                  sim_stratum_hazards = sim_stratum_hazards,
                  full_dyn_treatment_dm = full_dyn_treatment_map_dm) %>% 
      group_by(deworming_day, dyn_treat_days) %>% 
      summarize_at(vars(prob_deworm, alpha), funs(weighted.mean(., day_size))) %>% 
      ungroup()
    
    if (!is.null(other_group_suffix)) {
      other_group_fun <- partial(est_group_ate, 
                                 iter_est = iter_est, sim_stratum_hazards = sim_stratum_hazards, 
                                 treatment_map_dm = treatment_map_dm, dyn_treatment_mask_map = dyn_treatment_mask_map, 
                                 full_dyn_treatment_map_dm = full_dyn_treatment_map_dm,
                                 group_suffix = "left", other_group_suffix = NULL)
      ate_group %>% 
        plyr::ddply(c("phone_owner", "dist.pot.group", stringr::str_subset(names(.), stringr::str_interp("(?<!_id)_${other_group_suffix}$"))), 
                    other_group_fun) %>% 
        select(deworming_day, dyn_treat_days, alpha, prob_deworm) %>% #, all_treatment_id) %>%
        inner_join(obs_est, c("deworming_day", "dyn_treat_days"), suffix = c("_other", "_this")) %>% 
        mutate(prob_deworm_diff = prob_deworm_other - prob_deworm_this) %>% 
        rename_(.dots = stringr::str_subset(names(.), "_(this|other)$") %>% 
                          setNames(., stringr::str_replace_all(., c("this$" = group_suffix, "other$" = other_group_suffix)))) 
    } else {
      return(obs_est)
    }
  }
  
  iter_est <- iter_parameters %>% 
    left_join(obs_meta_data, c("stratum_id", "cluster_id")) %>% 
    mutate(log_kappa_census_covar = select(., starts_with("stratum_census_covar_coef")) %>% 
             magrittr::multiply_by(census_covar_dm) %>% 
             rowSums()) %>% 
    select(-starts_with("stratum_census_covar_coef"))
  
  dyn_treatment_mask_map %<>% set_names(stringr::str_c("dyn_treatment_", names(.)))
  
  identify_treatment_id(ate_info, treatment_map) %>% 
    plyr::ddply(c("phone_owner", "dist.pot.group", stringr::str_subset(names(.), "(?<!_id)_right$")), 
                est_group_ate, 
   .parallel = inner_parallel, 
   .paropts = lst(.packages = c("tidyverse"),
                  .export = c("estimate_treatment_deworm_prob", "estimate_treatment_deworm_prob_dm")),
    iter_est = iter_est,
    sim_stratum_hazards = sim_stratum_hazards,
    treatment_map_dm = treatment_map_dm,
    dyn_treatment_mask_map = dyn_treatment_mask_map,
    full_dyn_treatment_map_dm = full_dyn_treatment_map_dm)
}

subgroup_combiner <- function(iter_data, outcome_var, subgroups = NULL, 
                              group_by_regex = "(incentive_treatment|(private|social)_value|sms.treatment.2)") {
  iter_data %>% 
    group_by_at(c(vars(matches(group_by_regex), iter_id), subgroups)) %>% 
    summarize_at(vars(outcome_var), funs(wtd_iter_est = weighted.mean(., treatment_size))) %>% 
    ungroup() 
}


# Dynamic ATE -------------------------------------------------------------

get_dyn_ate <- function() {
  dyn_sms_control_ate <- tribble(
      ~ private_value_left, ~ social_value_left, ~ private_value_right, ~ social_value_right,
      "control",            "ink",               "control",             "control",
      "calendar",           "control",           "control",             "control",
      "calendar",           "bracelet",          "control",             "control",
      "calendar",           "bracelet",          "calendar",            "control",
      "control",            "bracelet",          "control",             "control" 
    ) %>% 
    mutate(signal_observed_left = social_value_left,
           signal_observed_right = social_value_right) %>% 
    bind_rows(
      tribble(
        ~ private_value_left, ~ social_value_left, ~ signal_observed_left, ~ private_value_right, ~ social_value_right, ~ signal_observed_right,
        "control",            "bracelet",          "control",              "control",             "control",            "control",
        "control",            "ink",               "control",              "control",             "control",            "control",
        "control",            "bracelet",          "bracelet",             "control",             "bracelet",           "control",
        "control",            "ink",               "ink",                  "control",             "ink",                "control"
    )) %>%
    mutate(incentive_shift_left = private_value_left,
           incentive_shift_right = private_value_right) %>% 
    bind_rows(
      tribble(
        ~ private_value_left, ~ incentive_shift_left,  ~ private_value_right, ~ incentive_shift_right,
        "calendar",           "calendar",              "calendar",            "control",
        "calendar",           "control",               "control",             "control" 
      ) %>% 
      mutate(social_value_left = "control",
             social_value_right = "control",
             signal_observed_left = "control",
             signal_observed_right = "control")) %>%
    mutate(dyn_dist_pot_left = NA,
           dyn_dist_pot_right = NA) %>% 
    bind_rows(
      tribble(
        ~private_value_left, ~ social_value_left, ~ dyn_dist_pot_left, ~ private_value_right, ~ social_value_right, ~ dyn_dist_pot_right,
        "control",           "control",              "close",          "control",             "control",            "far",
        "control",           "control",              "far",            "control",             "control",            "close",
        "calendar",          "control",              "close",          "calendar",            "control",            "far",
        "calendar",          "control",              "far",            "calendar",            "control",            "close",
        "control",           "bracelet",             "close",          "control",             "bracelet",           "far",
        "control",           "bracelet",             "far",            "control",             "bracelet",           "close",
        "control",           "ink",                  "close",          "control",             "ink",                "far",
        "control",           "ink",                  "far",          "control",             "ink",                  "close"
      ) %>% 
        mutate(signal_observed_left = social_value_left,
               signal_observed_right = social_value_right,
               incentive_shift_left = private_value_left,
               incentive_shift_right = private_value_right) 
    ) %>% 
    bind_rows("close" = ., "far" = ., .id = "dist.pot.group") %>% 
    bind_rows(`TRUE` = ., `FALSE` = ., .id = "phone_owner") %>%
    bind_rows(`TRUE` = ., `FALSE` = ., .id = "name_matched") %>%
    mutate(sms.treatment.2_left = "control", 
           sms.treatment.2_right = "control",
           # name_matched = FALSE,
           # hh.baseline.sample = FALSE,
           reminder_info_stock_left =  sms.treatment.2_left,
           reminder_info_stock_right = sms.treatment.2_right) %>% 
    mutate_at(vars(phone_owner, name_matched), as.logical) %>%
    mutate_at(vars(starts_with("dyn_dist_pot")), funs(coalesce(., dist.pot.group))) 
    # filter(!name_matched) %>% select(-name_matched)
  
  dyn_phone_owners_ate <- tribble(
    ~ private_value_left, ~ social_value_left, ~ sms.treatment.2_left, ~ private_value_right, ~ social_value_right, ~ sms.treatment.2_right,
    "control",            "control",           "reminder.only",        "control",             "control",            "control",
    "control",            "control",           "social.info",          "control",             "control",            "control",
    "control",            "control",           "social.info",          "control",             "control",            "reminder.only",
  
    "control",            "ink",               "social.info",          "control",             "control",            "control",
    "control",            "ink",               "social.info",          "control",             "control",            "social.info",
    "control",            "ink",               "social.info",          "control",             "ink",                "control",
  
    "calendar",           "control",           "social.info",          "control",             "control",            "control",
    "calendar",           "control",           "social.info",          "control",             "control",            "social.info",
    "calendar",           "control",           "social.info",          "calendar",            "control",            "control",
  
    "calendar",           "bracelet",          "social.info",          "control",             "control",            "control",
    "calendar",           "bracelet",          "social.info",          "control",             "control",            "social.info",
    "calendar",           "bracelet",          "social.info",          "calendar",            "control",            "social.info",
    "calendar",           "bracelet",          "social.info",          "calendar",            "bracelet",           "control",
    
    "control",           "bracelet",           "social.info",          "control",             "control",            "control",
    "control",           "bracelet",           "social.info",          "control",             "control",            "social.info",
    "control",           "bracelet",           "social.info",          "calendar",            "control",            "social.info",
    "control",           "bracelet",           "social.info",          "control",             "bracelet",           "control"
  ) %>%
    mutate(#incentive_shift_left = private_value_left,
           #incentive_shift_right = private_value_right,
           signal_observed_left = social_value_left,
           signal_observed_right = social_value_right,
           reminder_info_stock_left = sms.treatment.2_left,
           reminder_info_stock_right = sms.treatment.2_right) %>% 
    bind_rows(
      tribble(
        ~ social_value_left, ~ social_value_right, ~ sms.treatment.2_right,
        "bracelet",          "bracelet",           "control",
        "bracelet",          "control",            "social.info",
        "ink",               "ink",                "control",    
        "ink",               "control",            "social.info"
      ) %>% 
      mutate(private_value_left = "control",
             private_value_right = "control",
             signal_observed_left = "control",
             signal_observed_right = "control",
             reminder_info_stock_left = "control",
             reminder_info_stock_right = "control",
             sms.treatment.2_left = "social.info")) %>% 
    bind_rows(
      tribble(
        ~ private_value_left, ~ private_value_right, ~ sms.treatment.2_right,
        "calendar",           "calendar",            "control",
        "calendar",           "control",             "social.info"
      ) %>% 
      mutate(social_value_left = "control",
             social_value_right = "control",
             signal_observed_left = "control",
             signal_observed_right = "control",
             reminder_info_stock_left = "control",
             reminder_info_stock_right = "control",
             sms.treatment.2_left = "social.info")) %>% 
    bind_rows(
      tribble(
        ~ social_value_left, ~ reminder_info_stock_left, ~ social_value_right,  ~ sms.treatment.2_right, 
        "bracelet",          "social.info",              "bracelet",            "social.info",   
        "bracelet",          "control",                  "bracelet",            "control",       
        "ink",               "social.info",              "ink",                 "social.info",   
        "ink",               "control",                  "ink",                 "control"       
    ) %>% 
      mutate(private_value_left = "control",
             private_value_right = "control",
             signal_observed_left = social_value_left,
             signal_observed_right = social_value_right,
             sms.treatment.2_left = "social.info",
             reminder_info_stock_right = "control")) %>% 
    bind_rows(
      tribble(
         ~ social_value_left, ~ private_value_left, ~ dyn_dist_pot_left, ~ dyn_dist_pot_right,
         "control",           "control",            "close",             "far",
         "control",           "calendar",           "close",             "far",
         "ink",               "control",            "close",             "far",
         "bracelet",          "control",            "close",             "far" 
      ) %>%
      mutate(
        private_value_right = private_value_left,
        social_value_right = social_value_left,
        signal_observed_left = "control",
        signal_observed_right = "control",
        sms.treatment.2_left = "social.info",
        sms.treatment.2_right = "social.info",
        reminder_info_stock_left = sms.treatment.2_left,
        reminder_info_stock_right = sms.treatment.2_right,
        # incentive_shift_left = private_value_left,
        # incentive_shift_right = private_value_right
      )
    ) %>% 
    bind_rows("close" = ., "far" = ., .id = "dist.pot.group") %>%
    mutate(
      phone_owner = TRUE,
      incentive_shift_left = private_value_left,
      incentive_shift_right = private_value_right,
      # incentive_shift_left = "control",
      # incentive_shift_right = "control",
      # hh.baseline.sample = FALSE, 
      name_matched = FALSE
    ) %>%
    mutate_at(vars(starts_with("dyn_dist_pot")), funs(coalesce(., dist.pot.group))) %>% 
    # select(-name_matched)
    bind_rows(filter(., sms.treatment.2_left == "control" & sms.treatment.2_right == "control") %>% mutate(name_matched = TRUE))
  
  dyn_all_ate <- bind_rows(dyn_sms_control_ate, dyn_phone_owners_ate) 
  
  assertthat::assert_that(all(map_lgl(dyn_all_ate, ~ !any(is.na(.)))))
  
  return(dyn_all_ate)
}
