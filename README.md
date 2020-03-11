# TakeUp: RCT Design, Data Preparation, and Analysis

This repo is used to manage the R and Stan code used to design the TakeUp (see [working paper](takeup_workingpaper.pdf) for details no the study), prepare survey data and conduct analysis.

Here are the main files:

* [rct-design-fieldwork/](rct-design-fieldwork/) contains all scripts used for study design (e.g., randomization) and assisting with fieldwork/study implementation.
* [instruments/](instruments/) contains all the survey instruments and their translations.
* [data/](data/) points to a **private** submodule with survey data. This data has not yet been publicly shared.
* [images/](images/) stores all images used in producing the working paper.
* [stan_models/](stan_models/) contains all the Stan models used for analysis.
* [multilvlr/](multilvlr/) points to a public submodule. This is only used for some utility functions.

* [prepare_analysis_data.R](prepare_analysis_data.R) this script is used to prepare the survey data for analysis.
* [run_stan_dist.R](run_stan_dist.R) is the script that runs the Bayesian analysis for the working paper's main results.
* [postprocess_dist_fit.R](postprocess_dist_fit.R) postprocesses the fit data produced by [run_stan_dist.R](run_stan_dist.R).
* [takeup_workingpaper.Rmd](takeup_workingpaper.Rmd) this is the R markdown file use to produce the working paper PDF file.
* [takeup_workingpaper.pdf](takeup_workingpaper.pdf) this is the latest version of the working paper.
* [takeup_pap.Rmd](takeup_pap.Rmd) this is the R markdown file used to generate the pre-analysis plan PDF.

