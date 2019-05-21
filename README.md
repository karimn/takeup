TakeUp
======

This subfolder contains all the code and data currently used in the TakeUp study.

## Pilot

* [forms_2a_2c.R](forms_2a_2c.R)
* [generate_pilot_indepth_forms.R](generate_pilot_indepth_forms.R )
* [ke_takeup_pilot_schools.R](ke_takeup_pilot_schools.R)
* [observability_table.Rmd](observability_table.Rmd)
* [pilot_summary.R](pilot_summary.R)
* [TakeUp Scoping Summary.Rmd](TakeUp Scoping Summary.Rmd)
* [takeup_pilot_util.R](takeup_pilot_util.R)
* [takeup_scoping.R](takeup_scoping.R)

## RCT

* [analysis_util.R](analysis_util.R), this script has all the functions used in various other scripts/notebooks to assist in analysis.
* [param_dynamic_analysis.Rmd](param_dynamic_analysis.Rmd), presenting analysis for the dynamic parametric Bayesian model.
* [param_static_analysis.Rmd](param_static_analysis.Rmd), presenting analysis for the static parametric Bayesian model.
* [preliminary_findings_report.Rmd](preliminary_findings_report.Rmd)
* [takeup_field_notebook.Rmd](takeup_field_notebook.Rmd), this notebook records all the study monitoring and data quality activities conducted during the study's field work. This is where the first stage of data preparation takes place, using the raw data we got from SurveyCTO.
* [prepare_analysis_data.R](prepare_analysis_data.R), second stage of data preparation for analysis.
* [run_stan.R](run_stan.R), this scripts runs the Baysian multilevel models.
* [secobeliefs.R](secobeliefs.R), this script runs the second-order Bayesian analysis. 
* [social_know.Rmd](social_know.Rmd), presenting analysis from social knowledge tables.
* [takeup-ea-presentation.Rmd](takeup-ea-presentation.Rmd)
* [takeup_analysis.Rmd](takeup_analysis.Rmd), old Frequentist analysis. 
* [takeup_analysis2.Rmd](takeup_analysis2.Rmd), latest Frequentist analysis.
* [takeup_bayesian_analysis.Rmd](takeup_bayesian_analysis.Rmd), old Baysian analysis.
* [takeup_pap.Rmd](takeup_pap.Rmd), pre-analysis plan.
* [takeup_rct_assign_clusters.R](takeup_rct_assign_clusters.R), functions used to randomly select clusters for the study while ensuring sufficient buffering.
* [takeup_rct_design.Rmd](takeup_rct_design.Rmd), so old power analysis code.
* [takeup_rct_prep.R](takeup_rct_prep.R)
* [takeup_rct_target_villages.R](takeup_rct_target_villages.R), selecting villages from clusters.
* [takeup_sms.R](takeup_sms.R), script used to send SMS messages to subjects.
* [takeup_stat_model.Rmd](takeup_stat_model.Rmd)
* [takeup_workingpaper.Rmd](takeup_workingpaper.Rmd), main results document.
* [data/](data/), processed data used in analysis.
* [raw-data/](raw-data/), data from SurveyCTO.
* [images/](images/), image files used in working paper.
* [instruments/](instruments/), survey instruments.
* [slurm_scripts/](slurm_scripts/), scripts used to run analysis on Berkeley servers.
* [stan_analysis_data/](stan_analysis_data/), post-processed Stan analysis data.
* [stanfit/](stanfit/), Stan fits.
* [stan_models/](stan_models/), Stan Bayesian models.
