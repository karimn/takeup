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
        cluster_obs_know_person_prop = mean(obs_know_person_prop, na.rm = TRUE)
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
        dewormed ~  frac_externality_knowledge + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

judgement_geo_fit = judge_analysis_data %>%
    feglm(
        dewormed ~  judge_score_dewor + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

indiv_knowledge_fit = analysis_data %>%
    feglm(
        dewormed ~ log(census_cluster_pop) + obs_know_person_prop  + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

cluster_knowledge_fit = analysis_data %>%
    feglm(
        dewormed ~ log(census_cluster_pop) + cluster_obs_know_person_prop  + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit"),
        cluster = ~county
    )

etable(
    indiv_knowledge_fit, 
    cluster_knowledge_fit, 
    order = "!assigned"
)


analysis_data %>%
    group_by(cluster.id) %>%
    summarise(
        cluster_obs_know_person_prop = unique(cluster_obs_know_person_prop), 
        takeup = mean(dewormed), 
        treatment = unique(assigned_treatment), 
        dist_group = unique(assigned_dist_group)
    ) %>%
    ggplot(aes(
        x = cluster_obs_know_person_prop, 
        y = takeup, 
        colour = interaction(treatment, dist_group)
    )) +
    geom_point() +
    geom_smooth(
        method = "lm", 
        colour = "black"
    )

etable(
    judgement_geo_fit,
    externality_fit, 
    indiv_knowledge_fit, 
    cluster_knowledge_fit,
    order = "!assigned", 
    keep = c("!assigned"),
    drop = "Constant",
    cluster = ~county
)


#### SOB Reason Interlude ####

endline.know.table.data = endline.know.table.data %>%
    mutate(second.order.reason = str_to_lower(second.order.reason) %>% str_remove_all(., "[[:punct:]]")) %>%
    mutate(
        category_sob_reason = case_when(
            str_detect(second.order.reason, "saw me|point of treatment|pot|together|(?=.*met)(?=.*treatment|pot)|not seen|(?=.*didnt|did not|dont)(?=.*see)(?=.*treatment|pot)|didnt meet|(?=.*saw me)(?=.*treatment|pot)") ~ "campaign",
            # str_detect(second.order.reason, "health|responsible") | (str_detect(second.order.reason, "(?=.*know)(?=.*me)") & !str_detect(second.order.reason, "didnt|did not|dont|doesnt")) ~ "type",
            str_detect(second.order.reason, "example|good|health|responsible|clean")  ~ "type",
            str_detect(
                second.order.reason, 
                "never  seen|not very|hardly|dk|dont|dont|far|doesnt know|not close|rare|husband|mother|wife|son|daughter|neighbor|neighbour|friend") ~ "relationship",
            str_detect(second.order.reason, "brace") ~ "bracelet",
            str_detect(second.order.reason, "(?<!th)ink") ~ "ink",
            str_detect(second.order.reason, "calen") ~ "calendar",
            str_detect(second.order.reason, "ask|told|tell|talk|(?=.*met)(?!.*treatment|pot)(?=.*met)") ~ "communication",
            str_detect(second.order.reason, "school|hospital|not around|sick|preg|busy") ~ "circumstances",
            str_detect(second.order.reason, "tradi|erbal") ~ "traditional",
            TRUE ~ "other"
        ) 
    )   

endline.know.table.data %>%
    filter(!is.na(second.order.reason)) %>%
    filter(str_detect(second.order.reason, "wise")) %>%
    count(second.order.reason) %>%
    arrange(-n) %>%
    pull(second.order.reason) %>%
    head(30)

endline.know.table.data %>%
    group_by(
        assigned.treatment, 
        dist.pot.group,
        category_sob_reason
    ) %>%
    summarise(
        n = n()
    ) %>%
    mutate(frac = 100*n / sum(n)) %>%
    filter(category_sob_reason != "other") %>%
    filter(
        !(category_sob_reason %in% c("traditional", "ink", "bracelet", "calendar"))) %>%
    ggplot(aes(
        y = category_sob_reason, 
        x = frac, 
        fill = assigned.treatment
    )) +
    geom_col(position = position_dodge(0.8)) +
    facet_wrap(~dist.pot.group)  +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        x = "Percent of Respondents", 
        title = "SOB Reason Distribution by Treatment and Distance Group"
    )

ggsave(
    "temp-plots/second-order-reason-distribution.pdf",
    width = 10, 
    height = 10
)

endline.know.table.data %>%
    unique() %>%
    group_by(category_sob_reason) %>%
    filter(category_sob_reason == "type") %>%
    select(second.order.reason) %>%
    sample_n(20, replace = TRUE) %>%
    unique() %>%
    pull()

tmp_endline_reasons %>%
    unique() %>%
    group_by(category_sob_reason) %>%
    sample_n(20, replace = TRUE) %>%
    unique() %>%
    write_csv(
        "temp-data/second-order-reason-check.csv"
    )
tmp_endline_reasons %>%
    unique() %>%
    group_by(category_sob_reason) %>%
    sample_n(20, replace = TRUE) %>%
    unique() %>%
    write_csv(
        "temp-data/second-order-reason-check.csv"
    )
