library(foreach)

dummy_generate_strat_reg_data <- function(.data, .formula, .strat.by, .strata.contrasts, .covariates) {
  rhs.vars <- terms(.formula) %>%
    delete.response() %>%
    all.vars()
  
  stratum.contrasts.list <- setNames(list(.strata.contrasts), .strat.by)
  
  if (.strat.by %in% names(.data)) {
    design.mat <- update.formula(.formula, as.formula(paste("~ . *", .strat.by))) %>%
      model.matrix(.data, contrasts.arg = stratum.contrasts.list) %>%
      magrittr::extract(, stringr::str_detect(colnames(.), .strat.by)) #%>%
    # set_colnames(stringr::str_replace_all(colnames(.), sprintf("(%s)\\[T\\.([^\\]]+)\\]", paste(rhs.vars, collapse = "|")), "\\2"))
  } else {
    design.mat <- model.matrix(.formula, .data, contrasts.arg = stratum.contrasts.list)
  }
  
  design.mat %<>%
    set_colnames(stringr::str_replace_all(colnames(.), sprintf("(%s)\\[T\\.([^\\]]+)\\]", paste(rhs.vars, collapse = "|")), "\\2"))
  
  # Handling the covariates on their own because we want to demean them for each stratum. That way the intercept has a clearer
  # interpretation
  if (!purrr::is_empty(.covariates)) {
    design.mat <- .data %>%
      transmute_(.dots = c(.strat.by, setNames(.covariates, paste0("covar_", .covariates)))) %>%
      group_by_(.strat.by) %>%
      do(strat.covar.mat = model.matrix(as.formula(sprintf("~ (%s) * %s", paste(paste0("covar_", .covariates), collapse = " + "), .strat.by)),
                                        data = .,
                                        contrasts.arg = stratum.contrasts.list) %>%
           magrittr::extract(, stringr::str_detect(colnames(.), paste0("covar_.*", .strat.by))) %>%
           t %>%
           magrittr::subtract(rowMeans(.))) %>%
      ungroup %>% {
        do.call(cbind, .$strat.covar.mat)
      } %>%
      t %>%
      cbind(design.mat, .)
    
    
    # covar.row.means <- covar.design.mat %>%
    #   group_by_(.strat.by) %>%
    #   do(strat.covar.means = rowMeans(.$strat.covar.mat[[1]])) %>%
    #   ungroup
    #
    # design.mat <- covar.design.mat %>% {
    #   do.call(cbind, .$strat.covar.mat)
    #   } %>%
    #   t %>%
    #   cbind(design.mat, .)
    
  }
  
  return(design.mat)
  
}

foo <- function(.data,
                                  .formula,
                                  .strat.by,
                                  .cluster,
                                  .covariates = NULL, ...) {
  stopifnot(length(.strat.by) == 1)
  stopifnot(is.factor(magrittr::extract2(.data, .strat.by)))
  
  clean.data <- .data %>%
    select_(.dots = c(all.vars(.formula), .strat.by, .cluster, .covariates)) %>%
    na.omit()
  
  # Let's make the naming of contrasts (factors) a bit easier to parse. Changing back to default on function exit
  old.contrasts <- getOption("contrasts")
  old.contrasts["unordered"] <- "contr.Treatment"
  old.options <- options(contrasts = old.contrasts)
  on.exit(options(old.options), add = TRUE)
  
  rhs.vars <- terms(.formula) %>%
    delete.response() %>%
    all.vars()

  strata.contrasts <- clean.data %>%
    select_(.dots = c(rhs.vars, .strat.by, .covariates)) %>% {
      strata.sizes <- model.matrix(as.formula(paste("~ ", .strat.by)), .) %>% colSums

      contr.Treatment(levels(.[[.strat.by]]), contrasts = FALSE) %>%
        magrittr::inset(1, , strata.sizes)
    }

  strata <- colnames(strata.contrasts) %>%
    stringr::str_replace("\\[(.+)\\]", "\\1")
  colnames(strata.contrasts)[1] <- ""

  strata.contrasts[1, ] %<>% magrittr::divide_by(.[1] - sum(.[-1]))
  strata.contrasts[1, -1] %<>% magrittr::multiply_by(-1)
   
   design.mat <- dummy_generate_strat_reg_data(clean.data, .formula, .strat.by, strata.contrasts, .covariates) %>%
     set_colnames(stringr::str_replace_all(colnames(.), setNames(c("", "(intercept)"), c(stringr::str_interp(":?${.strat.by}$"), "^$"))))

  y <- model.frame(.formula, clean.data) %>% model.response() #reg.data$response
  lm.fit(design.mat[1:5000, ], y[1:5000])
  # fm <- stats::lm.fit(design.mat, y)
  # 
  # na.coef <- is.na(fm$coefficients)
  # 
  # fm$coefficients %<>% magrittr::extract(!na.coef)
  # 
  # fm$strat.by <- .strat.by
  # fm$strata <- strata
  # fm$strata.constrasts <- strata.contrasts
  # fm$cluster <- unname(unlist(clean.data[, .cluster])) #reg.data$cluster
  # fm$cluster.var.name <- .cluster
  # fm$model <- cbind(y, design.mat[, !na.coef])
  # fm$formula <- .formula
  # fm$covariates <- .covariates
  # 
  # class(fm) <- "lm_strat"
  # 
  # return(fm)
}

# Test --------------------------------------------------------------------

xx <- analysis.data %>%  
  filter(monitored, sms.treatment.2 == "sms.control", !hh.baseline.sample) %>% 
  anti_join(outlier.cells, c("assigned.treatment","sms.treatment", "mon_status", "cluster.id", "phone_owner")) %>% 
  strat_mht(dewormed.any ~ assigned.treatment, strat_by = "county_dist_stratum", cluster = "cluster.id", covar = census.reg.covar, 
            hypotheses = c("ink", "calendar", "bracelet", "bracelet - calendar"), 
            num_resample = 5)

# All model ---------------------------------------------------------------

hypo_data <- tribble(
  ~ hypothesis,
  
  # close, non_phone
  "ink",                                                                      
  "calendar",                                                                
  "bracelet",                                                                
  "bracelet - calendar",                                                     
  
  # far, non_phone
  "ink + ink:far",                                                           
  "calendar + calendar:far",                                                 
  "bracelet + bracelet:far",                                                
  "bracelet - calendar + bracelet:far - calendar:far",                      
  
  # close, phone
  "ink + ink:phone_owner",                                                  
  "calendar + calendar:phone_owner",                                        
  "bracelet + bracelet:phone_owner",                                        
  "bracelet - calendar + bracelet:phone_owner - calendar:phone_owner",    
  
  "reminder.only",
  "social.info",
  "social.info - reminder.only",
  
  # "ink:social.info",
  # "calendar:social.info",
  # "bracelet:social.info",
  # "bracelet:social.info - calendar:social.info",
  # "bracelet:social.info - ink:social.info",
  # "calendar:social.info - ink:social.info",
  
  # far, phone
  "ink + ink:phone_owner + ink:far + ink:far:phone_owner",                  
  "calendar + calendar:phone_owner + calendar:far + calendar:far:phone_owner",
  "bracelet + bracelet:phone_owner + bracelet:far + bracelet:far:phone_owner",
  "bracelet - calendar + bracelet:phone_owner - calendar:phone_owner + bracelet:far - calendar:far + bracelet:far:phone_owner - calendar:far:phone_owner",
  
  "reminder.only + reminder.only:far",
  "social.info + social.info:far",
  "social.info - reminder.only + social.info:far - reminder.only:far"
  
  # "ink:social.info",
  # "calendar:social.info",
  # "bracelet:social.info",
  # "bracelet:social.info - calendar:social.info",
  # "bracelet:social.info - ink:social.info",
  # "calendar:social.info - ink:social.info",
)

analysis.data %>% 
  filter(!hh.baseline.sample) %>% 
  anti_join(outlier.cells, c("assigned.treatment","sms.treatment", "mon_status", "cluster.id", "phone_owner")) %>% 
  mutate(sms.treatment = sms.treatment.2,
         phone_owner = factor(phone_owner, levels = c(FALSE, TRUE), labels = c("non_phone_owner", "phone_owner"))) %>% 
  # run_strat_reg(dewormed.any ~ assigned.treatment * sms.treatment * mon_status * dist.pot.group * phone_owner, 
  #               .strat.by = "county", .cluster = "cluster.id", .covariates = census.reg.covar)
  strat_mht(dewormed.any ~ assigned.treatment * sms.treatment * mon_status * dist.pot.group * phone_owner, 
            strat_by = "county", cluster = "cluster.id", covar = census.reg.covar,
            hypotheses = hypo_data$hypothesis, num_resample = 100)

# Test treatment map ------------------------------------------------------

mht_analysis_formula <- dewormed.any ~ assigned.treatment * sms.treatment * mon_status * dist.pot.group * phone_owner

mht_analysis_data <- analysis.data %>% 
  filter(!hh.baseline.sample) %>% 
  anti_join(outlier.cells, c("assigned.treatment","sms.treatment", "mon_status", "cluster.id", "phone_owner")) %>% 
  mutate(sms.treatment = sms.treatment.2,
         phone_owner = factor(phone_owner, levels = c(TRUE, FALSE), labels = c("phone_owner", "non_phone_owner"))) 

mht_selected_treatment <- get_treatment_map(mht_analysis_data, mht_analysis_formula) %>% 
  filter(mon_status == "monitored", sms.treatment == "sms.control") %>% 
  arrange(phone_owner, dist.pot.group, assigned.treatment, sms.treatment)  

mht_treatment_map_dm <- get_treatment_map_design_matrix(mht_analysis_data, mht_analysis_formula, mht_selected_treatment) %>% 
  select_if(~ n_distinct(.x) > 1)

ate_hypo_data <- tribble(~ left, ~ right,
                         2,      1, 
                         3,      1, 
                         4,      1,
                         4,      3) 

test_mat <- ate_hypo_data %$% subtract(mht_treatment_map_dm[left, ], mht_treatment_map_dm[right, ]) %>% as.matrix()

mht_analysis_data %>%
  run_strat_reg(mht_analysis_formula,
                .strat.by = "county", .cluster = "cluster.id", .covariates = census.reg.covar)

mht_analysis_data %>%
  strat_mht(mht_analysis_formula,
            strat_by = "county", cluster = "cluster.id", covar = census.reg.covar,
            hypotheses = test_mat, num_resample = 100)

mht_analysis_data %>% 
  strat_mht(mht_analysis_formula,
            strat_by = "county", cluster = "cluster.id", covar = census.reg.covar,
            hypotheses = c("ink", "calendar", "bracelet", "bracelet - calendar"), num_resample = 100)

# mht examples ------------------------------------------------------------

read_csv("data.csv") %>% 
  transmute(amount, treatment = factor(ratio, levels = 0:3), dummy_stratum = factor(1), cluster = 1:n()) %>% 
  # run_strat_reg(amount ~ treatment, .cluster = "cluster", .strat.by = NULL) %>% 
  # linear_tester(c("treatment[T.1]", "treatment[T.2]", "treatment[T.3]"))
  # tidy()
  strat_mht(amount ~ treatment, cluster = "cluster", strat_by = NULL, covar = NULL,
            hypotheses = c("treatment[T.1]", "treatment[T.2]", "treatment[T.3]"),
            num_resample = 30) 
