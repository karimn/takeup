---
title: "Analysis steps and outputs"
---

* [ ] Submit full questionnaires for census, baseline, SMS recruitment, endline survey.  
* PoT and cluster selection: 
  - [X] Describe selection of 150 PoT and selected clusters (criteria used, algorithm). 
  - Output: 
        - [X] Map that shows 3 countries, 150 PoTs, their catchment areas; 
        - [X] Actual distance distribution - slightly different since could not work with schools? PoTs were moved? 
        - [X] Karim can you please give me a sentence stating why the distance from PoT is random i.e. we overcome traditional selection problem that people far from e.g. a clinic/school are different from those close to a clinic/school?
        - [ ] Explain: people were monitored at PoT came from set of selected villages that have a buffer so that people are further away from the next PoT than from their own. closest villages is 2.5km, median 4.2km, average 4.2km, max 7.6 distance from next PoT. whereas distance from your own PoT median/mean is 1.1km, 44% are below 1.25km, max to your 2.9km, min to other 2.5km. this allows us to limit spillovers i.e. people switching their assignment (SUTVA).
        - [ ] Explain how we avoided areas where the pilot was
* [X] Explain stratified assignment mechanism
* Census: 
  - [ ] Describe random sampling of 
    + baseline (purely random sample), 
    + SMS (phone owners), 
    + PoT monitored and endline sample (stratified phone and non-phone owners). 
    + collected and will use GPS coordinates for each house.
  - [ ] Output: use GPS data as controls and for heterogeneity analysis?
    + [ ] Model density into take-up specification, create density measure at cluster level i.e. how closely people live. use that as (a) covariate and (b) to explore heterogeneity in take-up (Karim please confirm we will do that?) since density is likely correlated with people’s likelihood of meeting others i.e. the potential of observing the signal and therefore the effectiveness of social signals; 
    + [ ] model density into social knowledge and aggregate beliefs specification, create social proximity measure for each individual link (i.e. respondent A and 10 other individuals given to her in social knowledge table) —> control for social proximity or explicitly estimate heterogeneity (Karim ?) in knowledge and accuracy of beliefs given differential impact of signaling with varying social proximity. B&T model assumes more accurate knowledge about other individuals’ actions, or at least that people believe they have more accurate knowledge. —> look at individuals who are more or less “central” and people that live more or less “far/close” to me.

    ANNE: indicate these estimations for endline data. 

* [ ] Baselines survey: use to test for balance in relevant characteristics and report their means and std, prior to sensitization. 
  - [ ] Balance for all baseline characteristics: 
    + general characteristics: 
      - level of education 
      - gender 
      - floor 
      - age
      - ethnicity 
      - religion 
      - phone ownership 
    + specific ones, we added to baseline to ensure balance (likely correlated with outcome)
    + knowledge about deworming and externalities of deworming: know_deworm who_worms who_worms_other effect_worms effect_worms_other how_spread how_spread_other effect_worms effect_worms_other stop_worms stop_worms_other when_treat when_treat_other worms_affect w_affect_how neighbours_worms_affect n_w_affect_how
    + beliefs about other people’s deworming choice: a) dworm_rate, b) ink_dworm_rate
    + stigma and praise: yes/no, between 0 and 10 - compare deworming stigma/praise to immunization. comparison across different outcomes and show that similar concerns. reasons as well.  
* SMS treatment: we document the following for this sub-intervention (KEY: intervention was at HH level)
  - From baseline survey: Show balance table for phone and non-phone owners. Make clear difference between 2 population groups. 
  - Recruitment: i) percentage of adults that took up reminder SMS, ii) percentage of adults that took up social info SMS —write one sentence on that, no table or so necessary.  
  - Received SMS: percentage of people for whom i) failed, ii) sent, iii) successful. 
* Information take-up: if we find no treatment effect of text messages, we need to know whether it is b/c people didn’t read messages or reminder/social info treatment truly had no effect. 2 measures to capture if people read the messages. i) information rewards: describe treatment, percentage of people responding to message, ii) endline data: sms (indicator sms.control/reminder.only/social.info),  receive_text (yes/no, among those that were part of SMS treatment), number_text, text_content (text_content_other look at that later since needs coding)
  - Output: How different/similar were phone owners to non-phone owners? Among phone-owners: what was take up of SMS treatment? how many SMS were we able to deliver successfully? did people read the messages? 
* Sensitization on deworming and externalities (Anne: need to pull info on that for presentation) ignore monitoring data from research enumerators. focus on endline data i.e. self-reported answers. 
  Relevant variables to verify: 
  i) were informed about treatment and how: find_out (this was coded but multiple options given)
  ii) people were sensitized by CHV: cv_visit —> test for balance to make sure all had same sensitization
  iii) received info about incentives: flyer 
  iv) informed about dates of deworming: day_treat_begin, day_treat_end,  days_available (knowledge about no of days available), treat_days (sum it — people knew! :))
  v) where treatment was offered (only for wave 2, so do later): where_offered
  vi) treatment for free: pay
  - Output: all numbers/narrative. don’t think we need plots or tables for this stuff. 

6. PoT treatment take-up data: 
Specifications to run: 
0. We are not pulling entire data but I’d like to report # of people we dewormed in total and by arm. 
1. Is there any effect of incentives on deworming take-up conditional on sensitization, among all monitored that did NOT get SMS? Output: plot with confidence intervals. —> compare control to incentives. switch to actual research question: can social signaling affect deworming decisions? Compare ink/control and calendar/bracelet in same plot. 
2. Is there any effect of social incentives on deworming take-up, among phone-owner sample but no SMS? Output: redo plot from 1. but only for phone-owners this time. 1-4-6-8
(Note: can also say sth about difference btw phone and non-phone owners with regards to take-up)
3. Social signaling identification: control for reminder and social learning: 3-5-7-9. 
4. 1vs.2 and 2vs.3 as described further down. 

Heterogeneity: by distance i.e. “close” vs. “far”. Re-do take-up graph with 8 bars that show for each intervention arm take up for close and far communities. 

Robustness: test for spillover effects, people switching assignment. 

2nd measure of take-up behavior: distribution of take-up over time among people that did not receive SMS treatment. apply Kolmogorov-Smirnov test to check for differences in take-up over time. assumption  (optional)


7. Endline data:
1. Show balance on same characteristics (minus belief, stigma/punishment questions).
2. Verify assumptions and important descriptives: 
- Calendar is valid control for bracelet/social signaling: private value of calendar > private value of bracelet. Output: choice test, by treatment (treat_code for control/ink/bracelet/calendar) show demand for calendar/bracelet, split by people that already have a calendar and people who don’t.
- Comparison of two signals and private incentive:
i) visibility: seen_bracelet seen_ink seen_cal
ii) retention: 
— bracelet: got_bracelet (received when came to PoT), have_bracelet and wear_bracelet, not_wear (reason why not), seen_bracelet (enumerator verifies that seen bracelet), bad_bracelet (complaint by respondent)
— ink: got_ink, ink_visible (is ink still visible), hh_ink, number_ink, bad_ink (complaints)
— calendar: got_cal (received when coming for treatment), have_cal, confirm_cal (enumerator checks if has cal), bad_cal (complaint) 
iii) meaning: bracelet_meaning ink_meaning cal_meaning
iv) hh demand “complementarities/substitues”: hh_bracelet and number_bracelet, hh_cal and number_cal (who else has cal)

—> Calendar as control for bracelet based on fact: calendar is privately more valued, particularly look among those that did not come for treatment: how many had seen calendar? 
—> Ink less useful signal: less visible, fewer people seen it. 

3. Outcomes in addition to take-up: 
i) aggregate beliefs: how accurate were people’s beliefs about take-up at community level? further, show differences in beliefs between base- and endline. 
ii) individual level knowledge: 2 separate tests - person by person, pairs. check for don’t knows. first, greater certainty about other people’s actions. second, greater accuracy in knowledge about other people’s actions. 
—> check if people systematically under- or overestimate other’s actions. 

Karim, should we distinguish between SMS social info and everyone else (i.e. SMS control/reminder, non-phone owners) here? 
Karim, first thing I’d like to show is that people claim to know better i.e. we have a reduction in the DKs. 
