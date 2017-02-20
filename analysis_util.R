# ---- name.match.monitored
# Name matching function. Given the census data and take-up from a particular attempt to 
# find individuals whose names were recorded at the PoT.
name.match.monitored <- function(census.cluster.data, 
                                takeup.cluster.data, 
                                max.cost = 1) { # This is the maximum number of "edits" or difference allowed between names
  dist.mat <- adist(census.cluster.data$name1st, takeup.cluster.data$name1st, ignore.case = TRUE) %>%
    magrittr::add(adist(census.cluster.data$last_name, takeup.cluster.data$last_name, ignore.case = TRUE)) 
  
  census.cluster.data %>% 
    mutate(min.name.match.dist = aaply(dist.mat, 1, . %>% min(na.rm = TRUE)) %>% na_if(Inf),
           which.min.name.match.dist = ifelse(!is.na(min.name.match.dist), 
                                              aaply(dist.mat, 1, . %>% 
                                                      which.min %>% 
                                                      magrittr::extract(takeup.cluster.data$KEY.survey.individ, .)),
                                              NA),
           dewormed.matched = !is.na(which.min.name.match.dist) & min.name.match.dist <= max.cost,
           which.min.name.match.dist = ifelse(dewormed.matched, which.min.name.match.dist, NA))
}
# ---- end
# Misc Functions and Constants ----

reg.covar <- c("school", "floor", "ethnicity", "sms.ctrl.subpop", "age", "gender")

prepare.consent.dewormed.data <- function(.all.endline.data, .reconsent.data) {
  list(endline.survey = .all.endline.data, 
       reconsent = .reconsent.data) %>% 
    map_df(. %>% select(KEY.individ, monitor.consent, dewormed.reported), .id = "data.source") %>% 
    filter(!is.na(monitor.consent), !is.na(dewormed.reported)) %>%  
    group_by(KEY.individ) %>%  
    summarize(monitor.consent = any(monitor.consent), # Consider as reconsented if at least one acceptance
              dewormed.reported = ifelse(n_distinct(dewormed.reported) == 1, first(dewormed.reported), NA)) %>% # Multiple contradictory responses
    ungroup
}

prepare.analysis.data <- function(.census.data, .takeup.data, .endline.data, .consent.dewormed.reports, .cluster.strat.data) {
  dewormed.day.data <- .takeup.data %>% 
    filter(!is.na(KEY.individ)) %>% 
    group_by(KEY.individ) %>% 
    summarize(dewormed.day = min(deworming.day)) %>% # If dewormed multiple times, take the first day only
    ungroup
  
  analysis.data <- .census.data %>% 
           filter(!is.na(wave)) %>%  # Remove clusters no longer in study
    left_join(.consent.dewormed.reports, "KEY.individ") %>% 
    mutate(dewormed = KEY.individ %in% .takeup.data$KEY.individ, # TRUE if individual found in take-up data
           dewormed = ifelse(monitored, dewormed, NA)) %>% # NA if not in the monitored group
    left_join(dewormed.day.data, "KEY.individ")
  
  analysis.data %>% 
    filter(is.na(dewormed) | !dewormed) %>% # For anyone in study with with unknown or negative deworming status
    group_by(cluster.id) %>% 
    do(name.match.monitored(., filter(takeup.data, cluster.id %in% unique(.$cluster.id)))) %>% 
    ungroup %>% 
    select(KEY.individ, dewormed.matched, ends_with("min.name.match.dist")) %>% 
    right_join(analysis.data, "KEY.individ") %>% 
    left_join(transmute(takeup.data, KEY.survey.individ, dewormed.day.matched = deworming.day), c("which.min.name.match.dist" = "KEY.survey.individ")) %>% 
    mutate(monitored = !is.na(wave) & monitored, # Remove those dropped from the study 
           dewormed.any = (!is.na(dewormed) & dewormed) | dewormed.matched,
           dewormed.day.any = if_else(!is.na(dewormed.day), dewormed.day, dewormed.day.matched), 
           baseline.sample = !is.na(baseline.sample.wave),
           gender = factor(gender, levels = 1:2, labels = c("male", "female"))) %>% 
    left_join(select(.endline.data, KEY.individ, age, school, floor, ethnicity, any.sms.reported), "KEY.individ") %>% 
    left_join(select(.cluster.strat.data, wave, county, cluster.id, dist.pot.group), c("wave", "county", "cluster.id")) %>% 
    `attr<-`("class", c("takeup_df", class(.)))
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

base.prepare.baseline.endline.data <- function(.data) { #, .census.data, .cluster.strat.data) {
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
                                                "every year", "never", "when symptoms", "hw says"))) %>% 
    mutate_at(vars(worms_affect, neighbours_worms_affect), funs(yes.no.factor(., .yes.no = 1:0))) %>% 
    mutate_at(vars(spread_worms), yes.no.factor) 
}

prepare.endline.data <- function(.data, .census.data, .cluster.strat.data) {
  .data %>% 
    filter(present, interview, consent) %>% 
    arrange(KEY.individ, SubmissionDate) %>% 
    group_by(KEY.individ) %>% # If more than one entry for an individual, take first one (there are 22 such individuals)
    filter(row_number() == 1) %>% 
    ungroup %>% 
    base.prepare.baseline.endline.data %>% 
    mutate_at(vars(know_deworm, chv_visit, flyer, any.sms.reported), funs(yes.no.factor(., .yes.no = 1:0))) %>% 
    mutate_at(vars(treat_begin, days_available, treat_end), funs(factor(., levels = c(1, 98), c("knows", "DK")))) %>% 
    mutate(treat_begin_date = ymd(sprintf("2016-%d-%d", month_treat_begin, day_treat_begin)),
           treat_end_date = ymd(sprintf("2016-%d-%d", month_treat_end, day_treat_end)),
           #where_offered = labelled(where_offered, c("somewhere else" = 0, "home" = 3, "DK" = 98)) %>% as_factor)
           find_out = multi.factor(find_out, 
                                   levels = c(1:9, 99), 
                                   labels = c("friend", "family", "chv", "elder", "church", "flyer", "poster", "enumerator", "baraza",
                                              "other"))) %>%
    left_join(select(.cluster.strat.data, cluster.id, assigned.treatment, dist.pot.group), c("cluster.id")) %>% 
    left_join(select(.census.data, KEY.individ, dist.to.pot, sms.ctrl.subpop), "KEY.individ")  
           # text_content = factor(text_content, levels = c(1:3, 99), labels = c("reminders", "when/where", "social info", "other")))
}

prepare.baseline.data <- function(.data) {
  .data %>% 
    base.prepare.baseline.endline.data %>% 
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


prepare.cluster.takeup.data <- function(.data) {
 .data %>% 
   filter(monitored, monitor.consent) %>% 
   select(county, dist.pot.group, cluster.id, assigned.treatment, sms.treatment, dewormed.any) %>% 
   unite(stratum, county, dist.pot.group, sep = " ") %>% 
   group_by(assigned.treatment, sms.treatment, stratum, cluster.id) %>% 
   summarize(takeup.prop = mean(dewormed.any)) %>% 
   group_by(assigned.treatment, sms.treatment, stratum) %>% 
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


# Plotting code -----------------------------------------------------------

know.bel.cat.plot <- function(var, .baseline.data = baseline.data, .endline.data = endline.data, na.rm = FALSE) {
  list(baseline = .baseline.data, endline = .endline.data) %>% 
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
        geom_col(position = "dodge") +
        coord_flip() +
        theme(legend.position = "bottom") +
        labs(x = "", y = "Proportion")
    }
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
    facet_grid(.facet.formula, scales = "free_x", space = "free_x") +
    labs(title = "Estimated Take-up in Response to Incentive Treatment",
         caption = "Intervals shown identify the 90% confidence intervals, estimated using cluster robust standard errors.\nIntervals test the null hypothesis of no difference in take-up from that in the reference group.") +
    theme(legend.position = "bottom")
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
         caption = "Intervals shown identify the 90% confidence intervals, estimated using cluster robust standard errors.\nRed intervals test the null hypothesis of no difference in take-up from that in the control group (with non SMS treatment).\nGreen and blue intervals test the null hypothesis of no difference from the SMS treatment lower on the same column.") +
    theme(legend.position = "bottom")
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
  interact.with <- "close"
  
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
    mutate(dist = rep(rep(c("Far", "Close"), 2), c(5, 5, 4, 4)),
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
    scale_fill_manual("", values = c("red", "black"), labels = c("Calendar vs Bracelet", "Control vs Ink")) +
    theme(legend.position = "bottom")
}
