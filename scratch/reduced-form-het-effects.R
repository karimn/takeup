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


analysis_data %>%
    count(num.individuals)


village.centers

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
        family = binomial(link = "probit")
    )

judgement_geo_fit = judge_analysis_data %>%
    feglm(
        dewormed ~  judge_score_dewor + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit")
    )

indiv_knowledge_fit = analysis_data %>%
    feglm(
        dewormed ~ obs_know_person_prop + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit")
    )

cluster_knowledge_fit = analysis_data %>%
    feglm(
        dewormed ~ cluster_obs_know_person_prop + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit")
    )

endline.know.table.data %>%
    colnames()

endline.know.table.data %>%
    select(rel.size)


endline.know.table.data %>%
    group_by(second.order) %>%
    count(second.order.reason) %>%
    arrange(-n) %>%
    print(n = 200)


endline.know.table.data %>%
    colnames()


endline.know.table.data %>%
select(second.order)

split_reason_table = endline.know.table.data %>%
    filter(!is.na(second.order.reason)) %>%
    filter(second.order != "prefer not say") %>%
    mutate(
        cleanish_second_order_reason = case_when(
            str_detect(second.order.reason, "don't|dont") ~ "dont know",
            str_detect(second.order.reason, "brace") ~ "bracelet",
            str_detect(second.order.reason, "ink") ~ "ink",
            str_detect(second.order.reason, "calen") ~ "calendar",
            str_detect(second.order.reason, "rare") ~ "rarely see",
            str_detect(second.order.reason, "told") ~ "told them",
            TRUE ~ "other"
        ) 
    ) %>%
    group_by(relationship, second.order) %>%
    count(cleanish_second_order_reason)  %>%
    arrange(-n) %>%
    group_by(cleanish_second_order_reason) %>%
    mutate(
        total = sum(n)
    ) %>%
    spread(relationship, n)  %>%
    arrange(-total) 

split_reason_table %>%
    filter(second.order == "yes") %>%
    print(n = 21)



split_reason_table %>%
    gather(variable, value, -second.order, -cleanish_second_order_reason, -total)  %>%
    mutate(variable = fct_collapse(
        variable,
        "extended family" = "extended family",
        "neighbor" = "neighbor", 
        "other" = c(
            "other", 
            "church", 
            "friend", 
            "hh member", "village member" )
    ) %>% fct_rev) %>%
    ggplot(aes(
        x = value,
        fill = second.order, 
        y = cleanish_second_order_reason
    )) +
    geom_col(position = position_dodge(0.8)) +
    facet_wrap(~variable, ncol = 1) +
    theme_bw() +
    labs(
        x = "Count",
        y = "Reason for SOB Response"
    )

ggsave(
    "temp-plots/second-order-reason-count.pdf", 
    width = 10, height = 10
)

reason_table = endline.know.table.data %>%
    filter(!is.na(second.order.reason)) %>%
    mutate(
        cleanish_second_order_reason = case_when(
            str_detect(second.order.reason, "don't|dont") ~ "dont know",
            str_detect(second.order.reason, "brace") ~ "bracelet",
            str_detect(second.order.reason, "ink") ~ "ink",
            str_detect(second.order.reason, "calen") ~ "calendar",
            str_detect(second.order.reason, "rare") ~ "rarely see",
            str_detect(second.order.reason, "told") ~ "told them",
            TRUE ~ "other"
        ) 
    ) %>%
    group_by(relationship) %>%
    count(cleanish_second_order_reason)  %>%
    arrange(-n) %>%
    group_by(cleanish_second_order_reason) %>%
    mutate(
        total = sum(n)
    ) %>%
    spread(relationship, n)  %>%
    arrange(-total)

reason_table

pr_reason_table = reason_table %>%
    ungroup() %>%
    mutate(
        across(where(is.numeric), ~100*.x/sum(.x, na.rm = TRUE))
    )

pr_reason_table

    mutate(

    )
    arrange(cleanish_second_order_reason) %>%






analysis_data %>%
    colnames()

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

