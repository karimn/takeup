library(tidyverse)
library(fixest)
library(broom)
library(marginaleffects)

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# stick to monitored sms.treatment group
# remove sms.treatment.2
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_sms_data <- analysis.data %>% 
  filter(mon_status == "monitored") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()



nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


analysis_data <- monitored_nosms_data %>% 
    mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group)

social_perception_data = baseline.data %>% 
  select(matches("^(praise|stigma)_[^_]+$")) %>% 
  gather(key = key, value = response) %>% 
  separate(key, c("praise.stigma", "topic"), "_") %>% 
  separate(topic, c("topic", "question.group"), -2) %>% 
  filter(!is.na(response))  %>%
  count(praise.stigma, topic, response) %>% 
  group_by(praise.stigma, topic) %>% 
  mutate(n = n/sum(n)) %>%
  mutate_at(vars(praise.stigma, response), ~ fct_relabel(factor(.), str_to_title)) %>% 
  mutate(topic = fct_recode(factor(topic), 
                            "Wearing/not wearing nice clothes to church" = "clothe",
                            "Use Latrine/open defecation" = "defecat",
                            "Deworming/not deworming during MDA" = "dewor",
                            "Immunize/not immunize children" = "immuniz")) 


externality_tab_df = baseline.data %>%
    select(worms_affect, neighbours_worms_affect) %>%
    group_by(worms_affect, neighbours_worms_affect) %>%
    summarise(
        n = n()
    )  %>%
    filter(!is.na(worms_affect)) %>%
    spread(
        worms_affect, 
        n
    ) %>%
    rename(
        "(neighbour)/worms_affect" = neighbours_worms_affect
    )  


baseline.data = baseline.data %>%
    mutate(
        any_externality_knowledge = worms_affect == "yes" | neighbours_worms_affect == "yes"
    )


baseline.data %>%
    group_by(
        cluster.id
    ) %>%
    summarise(
        frac_externality_knowledge = mean(any_externality_knowledge, na.rm = TRUE), 
        treatment = unique(assigned.treatment)
    ) %>%
    ggplot(aes(
        x = frac_externality_knowledge
    )) +
    geom_histogram(
        colour = "black"
    ) +
    theme_minimal() + 
    labs(
        title = "Fraction Externality Knowledge Per Community", 
        x = "Fraction Aware of Worm Externality"
    )
ggsave(
    "temp-plots/hist-frac-aware-externalities.pdf", 
    width = 8, 
    height = 6
)


#### Externalities ####
cluster_clean_externality_df  = baseline.data %>%
    group_by(
        cluster.id
    ) %>%
    summarise(
        frac_externality_knowledge = mean(any_externality_knowledge, na.rm = TRUE)
    ) %>%
    inner_join(
        analysis_data %>%
            group_by(cluster.id) %>%
            summarise(
                takeup = mean(dewormed), 
                assigned_treatment = unique(assigned_treatment), 
                assigned_dist_group = unique(assigned_dist_group)
            ) %>%
            select(
                cluster.id, 
                assigned_treatment, 
                assigned_dist_group, 
                takeup
            ), 
    by = "cluster.id"
    )


clean_externality_df = analysis_data %>%
    left_join(
        baseline.data %>%
            group_by(
                cluster.id
            ) %>%
            summarise(
                frac_externality_knowledge = mean(any_externality_knowledge, na.rm = TRUE)
            ) %>% 
            ungroup() %>%
            mutate(
                above_med_externality_knowledge = frac_externality_knowledge > median(frac_externality_knowledge, na.rm = TRUE)
            ),
    by = "cluster.id"
    )

#### Perceptions ####
clean_perception_data = baseline.data %>% 
  select(cluster.id, matches("^(praise|stigma)_[^_]+$")) %>% 
  gather(key = key, value = response, -cluster.id) %>% 
  separate(key, c("praise.stigma", "topic"), "_") %>% 
  separate(topic, c("topic", "question.group"), -2)  %>%
  filter(!is.na(response))  


long_topic_judgement_score_df = clean_perception_data %>%
        count(cluster.id, praise.stigma, topic, response)  %>%
        group_by(cluster.id, praise.stigma, topic) %>% 
        mutate(n = n/sum(n))   %>%
        group_by(topic, cluster.id) %>%
        summarise(
            judge_score = n[praise.stigma == "praise" & response == "yes"] * n[praise.stigma == "stigma" & response == "yes"],
            stigma_score = n[praise.stigma == "stigma" & response == "yes"],
            praise_score = n[praise.stigma == "praise" & response == "yes"]
        )  %>%
        ungroup() %>%
        group_by(topic) %>%
        mutate(
            median_judge_score = median(judge_score, na.rm = TRUE),
            median_stigma_score = median(stigma_score, na.rm = TRUE)
        )


overall_judgement_score_df = clean_perception_data %>%
  count(cluster.id, praise.stigma, topic, response)  %>%
  group_by(cluster.id, praise.stigma, topic) %>% 
  mutate(n = n/sum(n))   %>%
  group_by(cluster.id) %>%
  filter(response == "yes") %>%
  summarise(
    judge_score = prod(n), 
  ) %>% 
  ungroup() %>%
  mutate(
    median_judge_score = median(judge_score, na.rm = TRUE)
  ) %>%
  mutate(topic = "overall")

judgement_score_df = bind_rows(
    long_topic_judgement_score_df, 
    overall_judgement_score_df
)

judgement_score_df %>%
    filter(topic == "dewor") %>%
    ggplot(aes(
        x = judge_score, 
    )) +
    geom_histogram(colour = "black") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        title = "Deworming 'Judgemental' Score Distribution"
    )

ggsave(
    "temp-plots/deworm-judgemental-score-histogram.pdf", 
    width = 10, 
    height = 10
)



wide_judgement_score_df = judgement_score_df %>%
    pivot_wider(
        names_from = topic, 
        values_from = c(judge_score, stigma_score, praise_score, median_judge_score, median_stigma_score)
    )


judge_analysis_data = left_join(
    analysis_data, 
    wide_judgement_score_df,
    by = "cluster.id"
) 

judge_analysis_data = judge_analysis_data %>%
    mutate(
        judgemental_dewor = judge_score_dewor >= median_judge_score_dewor,
        judgemental_immuniz = judge_score_immuniz >= median_judge_score_immuniz,
        judgemental_defecat = judge_score_defecat >= median_judge_score_defecat,
        judgmental_overall = judge_score_overall >= median_judge_score_overall,
    )

#### Community Centrality ####

analysis_data %<>% 
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
        knows_other_dewormed_no = sum(fct_match(.x$dewormed, "no"), na.rm = TRUE),
        thinks_other_knows = sum(fct_match(.x$second.order, c("yes", "no")), na.rm = TRUE),
        thinks_other_knows_yes = sum(fct_match(.x$second.order, "yes"), na.rm = TRUE),
        thinks_other_knows_no = sum(fct_match(.x$second.order, "no"), na.rm = TRUE),
      )
    }
  )) %>%
  mutate(
    in_knowledge_survey = map_lgl(knowledge_data, ~nrow(.x) > 0), 
    obs_know_person = if_else(in_knowledge_survey == FALSE, NA_real_, obs_know_person),
    obs_know_person_prop = if_else(in_knowledge_survey == FALSE, NA_real_, obs_know_person_prop)
  )

analysis_data = analysis_data %>%
    group_by(cluster.id) %>%
    mutate(
        cluster_obs_know_person_prop = mean(obs_know_person_prop, na.rm = TRUE),
        cluster_obs_know_person = mean(obs_know_person, na.rm = TRUE)
    ) %>%
    ungroup()



cluster_census_pop_df = census.data %>%
    filter(!is.na(cluster.id)) %>%
    group_by(
        cluster.id
    ) %>%
    summarise(
        census_cluster_pop = sum(num.individuals)
    )

analysis_data = analysis_data %>%
    left_join(
        cluster_census_pop_df, 
        by = "cluster.id"
    ) 

#### Ethnic Fractionalisation ####


analysis_data %>%
    filter(!is.na(ethnicity)) %>%
    group_by(
        cluster.id, ethnicity
    ) %>%
    mutate(
        n_ethnicity = n()
    ) %>%
    group_by(cluster.id) %>%
    mutate(
        ethnicity_share = n_ethnicity/sum(n_ethnicity)
    ) %>%
    summarise(
        fractionalisation = 1 - sum(ethnicity_share^2)
    ) %>%
    ggplot(aes(
        x = fractionalisation
    )) +
    geom_histogram(bins = 60, colour = "grey") +
    theme_minimal() +
    labs(
        x = "Fractionalisation", 
        y = "Count", 
        title = "Distribution of Fractionalisation Score Across Counties"
    )
ggsave(
    "temp-plots/frac-hist.pdf",
    width = 16,
    height = 8
)

ethnicity_share_df = analysis_data %>%
    filter(!is.na(ethnicity)) %>%
    group_by(
        cluster.id, ethnicity
    ) %>%
    mutate(
        n_ethnicity = n()
    ) %>%
    group_by(cluster.id) %>%
    mutate(
        ethnicity_share = n_ethnicity/sum(n_ethnicity)
    ) %>%
    summarise(
        fractionalisation = 1 - sum(ethnicity_share^2, na.rm = TRUE)
    )  %>%
    ungroup()

analysis_data %>%
    group_by(assigned_treatment, assigned_dist_group) %>%
    count(ethnicity) %>%
    filter(!is.na(ethnicity)) %>%
    ggplot(aes(
        y = assigned_treatment, 
        x = n, 
        fill = ethnicity
    )) +
    geom_col(position = position_dodge(0.9)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~assigned_dist_group, ncol = 1) +
    labs(
        x = "Count", 
        y = "Treatment", 
        title = "Ethnicity Counts By Condition"
    )

ggsave(
    "temp-plots/ethnicity-counts.pdf", 
    width = 10, 
    height = 10
    )



analysis_data = analysis_data %>%
    left_join(
        ethnicity_share_df, 
        by = "cluster.id"
    )



#### Regressions ####
probit_fit = analysis_data %>%
    feglm(
        dewormed ~ 0 + assigned_treatment:assigned_dist_group, 
        data = ., 
        family = binomial(link = "probit"), 
        cluster = ~county
    ) 

externality_fit = clean_externality_df %>%
    feglm(
        dewormed ~  0 + frac_externality_knowledge + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

judgement_fit = judge_analysis_data %>%
    feglm(
        dewormed ~  0 + judge_score_dewor + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

indiv_knowledge_fit = analysis_data %>%
    feglm(
        dewormed ~ 
            0 + 
            log(census_cluster_pop) + 
            obs_know_person  + 
            assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

cluster_knowledge_fit = feglm(
        data = analysis_data,
        dewormed ~ 
            0 + 
            log(census_cluster_pop) + 
            cluster_obs_know_person  + 
            assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )


ethnicity_fit =  analysis_data %>%
    feglm(
        dewormed ~ 
        0 + fractionalisation +
        assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"), 
        cluster = ~county
    )
    


term_dict = c(
        judge_score_dewor = "Judgemental score", 
        frac_externality_knowledge = "Externality knowledge",
        obs_know_person = "Number of people recognised in community",
        cluster_obs_know_person = "Number of people recognised in community, village mean",
        fractionalisation = "Fractionalisation", 
        deworming = "Dewormed"
    )

etable(
    judgement_fit,
    externality_fit, 
    indiv_knowledge_fit, 
    cluster_knowledge_fit,
    ethnicity_fit,
    order = "!assigned", 
    keep = c("!assigned"),
    drop = c("Constant", "census_cluster_pop"),
    cluster = ~county, 
    replace = TRUE, 
    depvar = FALSE,
    dict = term_dict,
    # style.tex = style.tex(var.title = "", fixef.title = "", stats.title = " ") 
    tex = TRUE, 
    file = "temp-data/rf-mechanism-full-regression-table.tex",
    title = "Takeup Heterogeneity: Full Sample", 
    notes = "All regressions include the saturated interaction of incentive and distance 
    conditions. Additional regressors are included without interactions i.e. as level shifts of 
    overall deworming. All regressions cluster standard errors at the county level and use frequentist 
    probit."
)

## Het TEs by social connectedness

indiv_knowledge_het_fit = analysis_data %>%
    feglm(
        dewormed ~ 
            log(census_cluster_pop) + 
            obs_know_person*assigned_treatment,
        family = binomial(link = "probit"),
        cluster = ~county
    )

cluster_knowledge_het_fit = feglm(
        data = analysis_data,
        dewormed ~ 
            log(census_cluster_pop) + 
            cluster_obs_know_person*assigned_treatment,
        family = binomial(link = "probit"),
        cluster = ~county
    )

etable(
    indiv_knowledge_het_fit, 
    cluster_knowledge_het_fit, 
    keep = "recognised",
    cluster = ~county, 
    replace = TRUE, 
    depvar = FALSE,
    dict = c(
        term_dict, 
        "assigned_treatmentink" = "Treatment: Ink", 
        "assigned_treatmentcalendar" = "Treatment: Calendar", 
        "assigned_treatmentbracelet" = "Treatment: Bracelet"
        ),
    tex = TRUE, 
    file = "temp-data/rf-mechanism-het-knowledge-regression-table.tex", 
    title = "Takeup Treatment Effect Heterogeneity: Social Connectedness", 
    notes = "Standard errors clustered at the county level, estimated using frequentist probit. 
    All regressions control for the log of village population and the interaction of `Number of people recognised' with 
    treatment. The first column only uses individuals included in the knowledge survey. The second column averages over 
    surveyed individuals in a village to produce a village level measure of social connectedness. 
    "
)
## Any variation in control group

control_externality_fit = clean_externality_df %>%
    filter(assigned_treatment == "control") %>%
    feglm(
        dewormed ~  0 + frac_externality_knowledge + assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

control_judgement_fit = judge_analysis_data %>%
    filter(assigned_treatment == "control") %>%
    feglm(
        dewormed ~  0 + judge_score_dewor + assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

control_indiv_knowledge_fit = analysis_data %>%
    filter(assigned_treatment == "control") %>%
    feglm(
        dewormed ~ 
            0 + 
            log(census_cluster_pop) + 
            obs_know_person  + 
            assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

control_cluster_knowledge_fit = feglm(
        data = analysis_data %>% filter(assigned_treatment == "control"),
        dewormed ~ 
            0 + 
            log(census_cluster_pop) + 
            cluster_obs_know_person  + 
            assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )


control_ethnicity_fit =  analysis_data %>%
    filter(assigned_treatment == "control") %>%
    feglm(
        dewormed ~ 
        0 + fractionalisation +
        assigned_dist_group,
        family = binomial(link = "probit"), 
        cluster = ~county
    )


etable(
    control_judgement_fit,
    control_externality_fit, 
    control_indiv_knowledge_fit, 
    control_cluster_knowledge_fit,
    control_ethnicity_fit,
    order = "!assigned", 
    keep = c("!assigned"),
    drop = c("Constant", "census_cluster_pop"),
    cluster = ~county,
    replace = TRUE, 
    depvar = FALSE,
    dict = term_dict,
    # style.tex = style.tex(var.title = "", fixef.title = "", stats.title = " ") 
    tex = TRUE, 
    file = "temp-data/rf-mechanism-control-regression-table.tex", 
    title = "Takeup Heterogeneity: Control Condition", 
    notes = "Standard errors clustered at the county level, models estimated using frequentist probit.
    Sample only includes those in the control condition. All regressions control for the covariate of interest 
    and the assigned distance group. Social connectedness regressions also control for the log of village population.
    "
)


## Regressions w/ Interaction for Any Treatment ##
clean_externality_df = clean_externality_df %>%
    mutate(
        any_incentive = assigned_treatment != "control"
    )
judge_analysis_data = judge_analysis_data %>%
    mutate(
        any_incentive = assigned_treatment != "control"
    )
analysis_data = analysis_data %>%
    mutate(
        any_incentive = assigned_treatment != "control"
    )

externality_any_fit = clean_externality_df %>%
    feglm(
        dewormed ~  0 + frac_externality_knowledge:any_incentive + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )


judgement_any_fit = judge_analysis_data %>%
    feglm(
        dewormed ~  0 + judge_score_dewor:any_incentive + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

indiv_knowledge_any_fit = analysis_data %>%
    feglm(
        dewormed ~ 
            0 + 
            obs_know_person:any_incentive  + 
            log(census_cluster_pop) + 
            assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

cluster_knowledge_any_fit = feglm(
        data = analysis_data,
        dewormed ~ 
            0 + 
            cluster_obs_know_person:any_incentive  + 
            log(census_cluster_pop) + 
            assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

fractionalisation_any_fit = analysis_data %>%
    feglm(
        dewormed ~  0 + fractionalisation:any_incentive + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

het_any_tests = list(
    car::lht(
        judgement_any_fit, 
        c(1, -1, rep(0, 8))
    ),
    car::lht(
        externality_any_fit, 
        c(1, -1, rep(0, 8))
    ),
    car::lht(
        indiv_knowledge_any_fit,
        c(0, 1, -1, rep(0, 8))

    ),
    car::lht(
        cluster_knowledge_any_fit,
        c(0, 1, -1, rep(0, 8))
    ),
    car::lht(
        fractionalisation_any_fit, 
        c(1, -1, rep(0, 8))
    )
)

het_pvals = map(het_any_tests, "Pr(>Chisq)") %>% map_dbl(., ~.x[[2]])

etable(
    judgement_any_fit,
    externality_any_fit, 
    indiv_knowledge_any_fit, 
    cluster_knowledge_any_fit,
    fractionalisation_any_fit,
    order = "!assigned", 
    keep = c("!assigned"),
    drop = c("Constant", "census_cluster_pop"),
    cluster = ~county,
    replace = TRUE, 
    depvar = FALSE,
    extralines = list("^_Homogeneous effects across `Any incentive', $p$-value" = het_pvals),
    dict = c(
        term_dict, 
        "any_incentiveTRUE" = "(Any incentive = True)", 
        "any_incentiveFALSE" = "(Any incentive = False)"
        ),
    tex = TRUE, 
    file = "temp-data/rf-mechanism-incentive-regression-table.tex", 
    title = "Takeup Heterogeneity: Any Incentive vs Control",
    notes = "Standard errors clustered at the county level, models estimated using frequentist probit.
     All regressions control for the covariate of interest interacted with an `any incentive' dummy as well as 
     the saturated interaction of treatment assignment with distance group. 
     Social connectedness regressions also control for the log of village population.
    "
)





#### SOB Reason Interlude ####

endline.know.table.data = endline.know.table.data %>%
    mutate(second.order.reason = str_to_lower(second.order.reason) %>% str_remove_all(., "[[:punct:]]")) %>%
    mutate(
        category_sob_reason = case_when(
            str_detect(second.order.reason, "brace") ~ "bracelet",
            str_detect(second.order.reason, "\\bink") ~ "ink",
            str_detect(second.order.reason, "calen") ~ "calendar",
            str_detect(second.order.reason, "saw me|point of treatment|pot|together|(?=.*met)(?=.*treatment|pot)|not seen|(?=.*didnt|did not|dont)(?=.*see)(?=.*treatment|pot)|didnt meet|(?=.*saw me)(?=.*treatment|pot)") ~ "campaign",
            # str_detect(second.order.reason, "health|responsible") | (str_detect(second.order.reason, "(?=.*know)(?=.*me)") & !str_detect(second.order.reason, "didnt|did not|dont|doesnt")) ~ "type",
            str_detect(second.order.reason, "example|good|health|responsible|clean")  ~ "type",
            str_detect(
                second.order.reason, 
                "never  seen|not very|hardly|dk|dont|dont|far|doesnt know|not close|rare|husband|mother|wife|son|daughter|neighbor|neighbour|friend") ~ "relationship",
            str_detect(second.order.reason, "ask|told|tell|talk|(?=.*met)(?!.*treatment|pot)(?=.*met)") ~ "communication",
            str_detect(second.order.reason, "school|hospital|not around|sick|preg|busy") ~ "circumstances",
            str_detect(second.order.reason, "tradi|erbal") ~ "traditional",
            TRUE ~ "other"
        ) 
    )   


knowledge_fit = function(data, dep_var) {
    data %>%
        filter(second.order %in% c("yes", "no")) %>%
        filter(!is.na(second.order)) %>%
        filter(!is.na(second.order.reason)) %>%
        mutate(
            lhs = category_sob_reason_short == dep_var
        ) %>%
        feols(
            data = ., 
            fml = lhs ~ assigned.treatment,
            split = ~second.order
        ) 
}


dep_vars = c(
    "signal",
    "campaign", 
    "communication", 
    "relationship",
    "type", 
    "circumstances",
    "other"
)

endline.know.table.data = endline.know.table.data %>%
    mutate(
    category_sob_reason_short =  if_else(
        category_sob_reason %in% c("bracelet", "ink", "calendar"), 
        "signal", 
        category_sob_reason
    )
    )


endline.know.table.data %>%
    write_csv(
        "temp-data/second-order-reason-raw-data.csv"
    )



endline.know.table.data %>%
    filter(!is.na(second.order)) %>%
    filter(!is.na(second.order.reason)) %>%
    filter(second.order != "don't know") %>%
    filter(category_sob_reason == "other") %>%
    count(
        second.order.reason
    ) 

endline.know.table.data %>%
    filter(!is.na(second.order)) %>%
    filter(!is.na(second.order.reason)) %>%
    filter(second.order != "don't know") %>%
    filter(category_sob_reason == "other") %>%
    count(
        second.order.reason
    ) %>%
    write_csv(
        "temp-data/sob-other-reasons.csv"
    )


gpt_df = read_csv("temp-data/sob-other-reasons-feedback-gpt-3.5-turbo-single-pass.csv")


gpt_df = gpt_df %>%
    mutate(
        gpt_category = str_remove_all(`gpt-3.5-turbo`, "[:punct:]") %>%
            str_to_lower() %>%
            str_remove_all("category|categorization|classification") %>%
            str_trim()
    ) %>%
    mutate(
        gpt_category_clean = case_when(
            gpt_category == "other" ~ 'other', 
            gpt_category == "campaign" ~ 'campaign',
            gpt_category == "type" ~ 'type',
            gpt_category == "relationship" ~ 'relationship',
            gpt_category == "circumstances" ~ 'circumstances',
            gpt_category == "communication" ~ 'communication',
            gpt_category == "signal" ~ 'signal',
            str_detect(gpt_category, "unsure") ~ "unsure",
            TRUE ~ "misc"
        ) 
    ) 


gpt_df

#### SOB Regressions ####


gpt_endline_know_table_data = bind_rows(
    endline.know.table.data %>%
        filter(!is.na(second.order)) %>%
        filter(!is.na(second.order.reason)) %>%
        filter(second.order != "don't know") %>%
        filter(category_sob_reason != "other"),
    left_join(
        endline.know.table.data %>%
            filter(!is.na(second.order)) %>%
            filter(!is.na(second.order.reason)) %>%
            filter(second.order != "don't know") %>%
            filter(category_sob_reason == "other"), 
        gpt_df %>%
            select(
                second.order.reason, 
                gpt_category_clean
            ), 
        by = "second.order.reason"
        )
    ) %>%
    mutate(
        category_sob_reason_short =  if_else(
            category_sob_reason_short == "other", 
            gpt_category_clean,
            category_sob_reason_short
        )
    ) %>%
    # If we can't clean even with GPT just filter 
    # (by default GPT will put in other but even then some it cannot tell at all) 
    filter(
        category_sob_reason_short %in% c(
            "signal",
            "campaign", 
            "communication", 
            "relationship",
            "type", 
            "circumstances",
            "other"
        )
    ) 



sob_reason_fits = map(
        dep_vars, 
        ~knowledge_fit(data = endline.know.table.data, dep_var = .x)
    )

sob_gpt_reason_fits = map(
        dep_vars, 
        ~knowledge_fit(data = gpt_endline_know_table_data, dep_var = .x)
    )

sob_fit_df = imap_dfr(sob_reason_fits, ~map_dfr(., ~tidy(.x, conf.int = TRUE) %>% mutate(N = nobs(.x)), .id = "sample") %>% mutate(var = dep_vars[[.y]]))
sob_gpt_fit_df = imap_dfr(sob_gpt_reason_fits, ~map_dfr(., ~tidy(.x, conf.int = TRUE) %>% mutate(N = nobs(.x)), .id = "sample") %>% mutate(var = dep_vars[[.y]]))



sob_fit_df %>%
    write_csv(
        "temp-data/second-order-reason-distribution.csv"
    )

sob_gpt_fit_df %>%
    write_csv(
        "temp-data/second-order-gpt-reason-distribution.csv"
    )
