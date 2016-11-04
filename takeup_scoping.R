library(dplyr)
library(plyr)
library(magrittr)
library(tidyr)
library(foreign)
library(ggplot2)
library(stringr)
library(tm)
library(wordcloud)
library(parallel)
library(foreach)
library(Hmisc)

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores) # This is necessary for "stemming" in the text mining part
}, error=function(err) {
  registerDoSEQ()
})

theme_set(theme_bw() + 
            theme(panel.border=element_rect(color=NA), 
                  axis.ticks=element_blank(), 
                  strip.background=element_rect(color=NA, size=2)))

transform.txt.to.factor <- function(.col, .levels) {
  .col %>% 
    sub(sprintf(".*(%s).*", paste(.levels, collapse="|")), "\\L\\1", ., ignore.case=TRUE, perl=TRUE) %>% 
    factor(levels=.levels, labels=capitalize(.levels))
}

transform.yes.no <- . %>% 
  transform.txt.to.factor(c("yes", "no"))

individ.data <- read.dta("~/Data/TakeUp/NP_TAKEUP_INDIVIDUAL_27112015.dta") %>% 
  transmute(have.bednet=ind5_have_bednet,
            children.dewormed=ind7_children_dewormed,
            know.anyone.no.child.deworm=transform.yes.no(ind9_doesnt_allow), 
            adults.can.get.worms=ind11_adults_get_worms,
            problems.deworming=ind14_deworming_pills_problems,
            purchased.drugs=ind15_medicine,
            travel.time.purchase.drugs=ind20_time_to_purchase,
            discuss.deworming=ind31_discuss_deworming,
            come.adult.deworming=ind34_come_deworming,
            peer.influence.deworming=ind38_decision_affected,
            want.deworming.signal=ind41_like_observe,
            want.to.signal.deworming=ind44_show_others,
            know.voting=ind47_know_ppl_voted,
            know.voting.ink.1=grepl("ink", ind48_know_how1, ignore.case=TRUE),
            know.voting.ink.2=grepl("ink", ind48_know_how2, ignore.case=TRUE),
            know.voting.ink=factor(know.voting.ink.1 | know.voting.ink.2, levels=c(TRUE, FALSE), labels=c("Yes", "No")),
            know.friend.voted=transform.yes.no(ind50_friend_voted), 
            care.others.know.voting=ind52_care_ppl_voted %>% 
              tolower %>% 
              sub("don't\\s+care", "no", .) %>% 
              sub("care|happy", "yes", .) %>% 
              transform.yes.no,
            want.to.signal.voting=transform.yes.no(ind53_show_participation), 
            want.voting.signal=transform.yes.no(ind54_know_other_participated),
            selected.incentive=sub(".*(bracelet|ribbon|ink).*", "\\L\\1", ind60_item_chose, ignore.case=TRUE, perl=TRUE) %>% 
              factor(levels=c("bracelet", "ribbon", "ink")))

individ.fgd.data <- read.dta("~/Data/TakeUp/NP_TAKEUP_FOCUS_27112015.dta") %>% 
  select(matches(sprintf("^fgd(%s)", paste(c(11, 13:14, 17, 22) * 100, collapse="|")))) %>% 
  set_names(names(.) %>% gsub("^fgd\\d+_", "", .) %>% gsub("_", ".", .)) %>% 
  dplyr::rename(adults.buy.deworming.meds=deworming.medicines,
                adults.free.treatment.schools=come.treatment,
                ok.not.contribute.community.work=doesnot.contribute)

teacher.fgd.data <- read.dta("~/Data/TakeUp/NP_TAKEUP_TEACHERS_27112015.dta") %>% 
  select(matches(sprintf("^tgd(%s)", paste(c(5, 14:15, 19) * 100, collapse="|")))) %>% 
  set_names(names(.) %>% gsub("^tgd\\d+_", "", .) %>% gsub("_", ".", .)) %>% 
  dplyr::rename(possible.deworm.adults.school=encourage.adults)

healthworker.data <-
  read.dta("~/Data/TakeUp/NP_TAKEUP_HEALTHWORKER_27112015.dta") %>% 
  select(matches(sprintf("^np(%s)", paste(c(3, 5, 6, 7, 10, 14, 19, 23, 28) + 100, collapse="|")))) %>% 
  set_names(names(.) %>% gsub("^np\\d+_", "", .) %>% gsub("_", ".", .)) %>% 
  dplyr::rename(deworm.children.affect.others=effect.dwm,
                dworm.healthworkers=pharm.dwm.pill,
                deworm.parents=parent.dwm,
                how.many.buy.meds=buyers,
                how.long.to.deworm.adults=day.dwm,
                use.incentives.other.programs=use.incentives) %>% 
  mutate(how.many.buy.meds=transform.txt.to.factor(how.many.buy.meds, c("almost all", "very few people", "many")),
         how.long.to.deworm.adults=transform.txt.to.factor(how.long.to.deworm.adults, c("one week", "one deworming day")))

individ.data %>% 
  gather %>% 
  filter(value %in% c("Yes", "No"), !key %in% c("know.voting", "come.adult.deworming")) %>%
  bind_rows(individ.fgd.data %>% 
              gather %>% 
              filter(value %in% c("Yes", "No"))) %>% 
  ggplot() +
  geom_bar(aes(x=key, fill=value), position="dodge") +
  scale_fill_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

individ.data %>% 
  filter(travel.time.purchase.drugs >= 0) %>% 
  ggplot() +
  geom_histogram(aes(travel.time.purchase.drugs), binwidth=5)

individ.data %>% 
  ggplot() +
  geom_bar(aes(x=selected.incentive))

individ.fgd.data %>% 
  gather %>% 
  filter(value %in% c("Yes", "No")) %>%
  ggplot() +
  geom_bar(aes(x=key, fill=value), position="dodge") +
  scale_fill_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

teacher.fgd.data %>% 
  gather %>% 
  filter(value %in% c("Yes", "No")) %>%
  ggplot() +
  geom_bar(aes(x=key, fill=value), position="dodge") +
  scale_fill_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

healthworker.data %>% 
  gather %>% 
  filter(value %in% c("Yes", "No")) %>%
  ggplot() +
  geom_bar(aes(x=key, fill=value), position="dodge") +
  scale_fill_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

healthworker.data %>% 
  select(how.many.buy.meds) %>% 
  gather %>% 
  ggplot() +
  geom_bar(aes(x=key, fill=value), position="dodge") +
  scale_fill_discrete("") +
  scale_x_discrete("") 

healthworker.data %>% 
  select(how.long.to.deworm.adults) %>% 
  gather %>% 
  ggplot() +
  geom_bar(aes(x=key, fill=value), position="dodge") +
  scale_fill_discrete("") +
  scale_x_discrete("")

## Text mining
get.response.term.docs <- function(.data, prefix, question.num, to.merge.question.num) {
  .data %>% 
    select(matches(sprintf("^%s(%s)", prefix, paste(union(question.num, to.merge.question.num), collapse="|")))) %>% 
    gather(question, value, convert=FALSE) %>%
    mutate(question=factor(ifelse(grepl(sprintf("^%s(%s)", prefix, paste(to.merge.question.num, collapse="|")), question), sub("\\d$", "", question), as.character(question)))) %>% 
    group_by(question) %>% 
    dplyr::summarize(text=paste(value, collapse=" ")) %>% 
    ungroup %>% 
    mutate(question=sub(sprintf("%s(\\d+)[a-b]?_(.+)", prefix), "\\2 (\\1)", question) %>% gsub("_", ".", .)) %>% 
    { VCorpus(DataframeSource(.), readerControl=list(reader=readTabular(mapping=list(content="text", id="question")))) } %>% 
    tm_map(stripWhitespace) %>% 
    tm_map(removePunctuation, preserve_intra_word_dashes=TRUE) %>% 
    tm_map(content_transformer(tolower)) %>% 
    tm_map(removeWords, stopwords("english")) %>% 
    # tm_map(stemDocument) %>% 
    TermDocumentMatrix %>% 
    as.matrix %>% 
    adply(.margins=2, function(column) {
      column[column > 1] %>% 
        data.frame(word=names(.), count=., stringsAsFactors=FALSE) %>% 
        arrange(desc(count)) %>%
        mutate(word=factor(word, levels=word)) 
      }) %>% 
    dplyr::rename(question=Docs) %>% 
    mutate(question.num=str_extract(as.character(question), "\\d+(?=\\))") %>% as.integer)
}
  
tdm.individ.data <- read.dta("~/Data/TakeUp/NP_TAKEUP_INDIVIDUAL_27112015.dta") %>% 
  get.response.term.docs("ind", c(21, 26:27, 36, 42, 43), c(28, 29, 39, 45))

read.dta("~/Data/TakeUp/NP_TAKEUP_FOCUS_27112015.dta") %>% 
  get.response.term.docs("fgd", NULL, c(c(5, 10, 15, 18:19, 21, 23:24) * 100, "2000a")) %>% 
  filter(question.num != 1000) %>% 
  mutate(question=ifelse(.$question.num == 1500, sprintf("why.not.buy.meds (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 1800, sprintf("why.not.treat.adults.schools (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 1900, sprintf("learn.about.important.events (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 2000, sprintf("health.promoters (%d)", question.num), as.character(question))) %>% 
  ggplot(aes(word, count)) +
  geom_bar(stat="identity") +
  scale_y_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(~ question, scales="free", ncol=3)

read.dta("~/Data/TakeUp/NP_TAKEUP_TEACHERS_27112015.dta") %>% 
  get.response.term.docs("tgd", NULL, c(6, 11, 16, 21) * 100) %>% 
  mutate(question=ifelse(.$question.num == 600, sprintf("why.deworm.adults (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 1100, sprintf("reasons.parents.refuse.deworm.children (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 1600, sprintf("reasons.parents.dont.come.deworm (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 2100, sprintf("suggestion.deworm.adults (%d)", question.num), as.character(question))) %>% 
  ggplot(aes(word, count)) +
  geom_bar(stat="identity") +
  scale_y_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(~ question, scales="free", ncol=3)

read.dta("~/Data/TakeUp/NP_TAKEUP_HEALTHWORKER_27112015.dta") %>% 
  get.response.term.docs("np", NULL, c(c(11, 20:22, 24, 25) + 100, "126b", "129b")) %>% 
  filter(!is.na(question.num), question.num != 124) %>% 
  mutate(question=ifelse(.$question.num == 111, sprintf("why.deworm.parents (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 126, sprintf("how.healthworker.encourage.adults (%d)", question.num), as.character(question))) %>% 
  mutate(question=ifelse(.$question.num == 129, sprintf("who.health.promoter (%d)", question.num), as.character(question))) %>% 
  ggplot(aes(word, count)) +
  geom_bar(stat="identity") +
  scale_y_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(~ question, scales="free", ncol=3)

tdm.individ.data %>% 
  filter(grepl("convenient.location", Docs)) %>% 
  ggplot(aes(word, count)) +
  geom_bar(stat="identity") +
  scale_y_discrete("") +
  scale_x_discrete("") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(~ question, scales="free", ncol=2)