# Analysis steps and outputs (notes from yesterday)

Note: we submit full questionnaires for census, baseline, SMS recruitment, endline survey.  

1. PoT and cluster selection: describe selection of 150 PoT and selected clusters (criteria used, algorithm). Output: (1) map that shows 3 countries, 150 PoTs, their catchment areas; (2) actual distance distribution - slightly different since could not work with schools? PoTs were moved? Karim can you please give me a sentence stating why the distance from PoT is random i.e. we overcome traditional selection problem that people far from e.g. a clinic/school are different from those close to a clinic/school?
Explain: people were monitored at PoT came from set of selected villages that have a buffer so that people are further away from the next PoT than from their own. closest villages is 2.5km, median 4.2km, average 4.2km, max 7.6 distance from next PoT. whereas distance from your own PoT median/mean is 1.1km, 44% are below 1.25km, max to your 2.9km, min to other 2.5km. this allows us to limit spillovers i.e. people switching their assignment (SUTVA).

2. Census: describe random sampling of i) baseline (purely random sample), ii) SMS (phone owners), iii) PoT monitored and endline sample (stratified phone and non-phone owners). collected and will use GPS coordinates for each house. 
Output: use GPS data as controls and for heterogeneity analysis?
(1) model density into take-up specification, create density measure at cluster level i.e. how closely people live. use that as (a) covariate and (b) to explore heterogeneity in take-up (Karim please confirm we will do that?) since density is likely correlated with people’s likelihood of meeting others i.e. the potential of observing the signal and therefore the effectiveness of social signals; 
(2) model density into social knowledge and aggregate beliefs specification, create social proximity measure for each individual link (i.e. respondent A and 10 other individuals given to her in social knowledge table) —> control for social proximity or explicitly estimate heterogeneity (Karim ?) in knowledge and accuracy of beliefs given differential impact of signaling with varying social proximity. B&T model assumes more accurate knowledge about other individuals’ actions, or at least that people believe they have more accurate knowledge. —> look at individuals who are more or less “central” and people that live more or less “far/close” to me.

ANNE: indicate these estimations for endline data. 

3. Baselines survey: use to test for balance in relevant characteristics and report their means and std, prior to sensitization. 
1) Balance for all baseline characteristics: 
— general characteristics: 
- level of education 
- gender 
- floor 
- age
- ethnicity 
- religion 
- phone ownership 
—  specific ones, we added to baseline to ensure balance (likely correlated with outcome)
- knowledge about deworming and externalities of deworming: know_deworm who_worms who_worms_other effect_worms effect_worms_other how_spread how_spread_other effect_worms effect_worms_other stop_worms stop_worms_other when_treat when_treat_other worms_affect w_affect_how neighbours_worms_affect n_w_affect_how
- beliefs about other people’s deworming choice: a) control, b) ink 
- stigma and praise: yes/no, between 0 and 10 - compare deworming stigma/praise to immunization 


SMS
balance table for phone and non-phone owners 
recruitment: reminder, social info check percentage of people that signed up. 
receipt: failed, sent, successful - sentence 
information take-up verified through a) reward pick up: percentage of people responding to message and b) endline data check knowledge about text received and content. 

Sensitization ~ had knowledge about availability of deworming, and incentives:
- day_treat_begin 
- day_treat_end
- 
- also: treatment free, CHV visited etc. 

PoT treatment take-up data: 
- take-up data PoT: is there any effect at all of incentives on deworming take-up? could be none. we are comparing control to each incentive. we are doing that on the sample of EVERYONE (SMS + non SMS treated), composed of the phone owners and non-phone owners (pure control), adjusting for oversampling.
—> conditional on sensitizing people, is there any effect of incentives? = sample of monitored (0 or 1, every monitored and know whether came or not)
output: (1) monitored people, separately for counties, distances and the overall 
(2) same as above for each of 12 days ~dynamics, we might explore further. assumption: static decision-making. non-SMS people. 
—> we stratified. b/c is that there might be no difference comparing control and calendar. 
- distance, counties - total sample. SMS treated, non SMS treated (incl. SMS control - phone owners and non-phone owners). total sample (ignore SMS) split by far and close. 

Graph 1: comparison control, ink, bracelet, calendar
—> we can a) see impact of incentives overall compared to control
—> b) we can see the impact of incentives 

—> Calendar is the control for the bracelet: based on “fact” that calendar is as or more liked than bracelet. Data: endline survey we allowed to choose. Data for the ink and control. Calendar has greater private value than the bracelet —> people if offered both should choose on average calendar more. 

Calendar lower take-up than the bracelet —> same private value, but is less visible, i.e. we are giving people something coming for deworming, normative. Visibility. “did you see bracelet/calendar?”

—> Main confounds: (1) reminded of deworming, (2) social learning, people observe others coming for deworming - receiving information about other people’s choices which could have affect their valuation of deworming or beliefs of its importance etc. - sth. that has nothing to do with a desire to display. 
—> confounds: explain reminder / social info

Descriptive statistics: 
- Sensitization: Knowledge about deworming and externalities 
- Externalities - baseline: balance, endline: after the sensitization. 3 questions.
- Knowledge: who can get storms, deworming pills, experience deworming 

Endline:
- Outcome: 
a) Aggregate belief: compare for every single person - sample, non-text message treatment people (control SMS people - phone, non-phone owners). aggregate beliefs: distribution test.
b) Individual knowledge: different people ~ could be wrong with all of them:
over- and underestimating of beliefs.
Don’t knows ~ across four arms. 
social knowledge surveys: too afraid. 

- Descriptives:
- ink/bracelet/calendar retention
- seen someone with bracelet/ink/calendar
- choice between calendar and bracelet
