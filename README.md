# TakeUp: RCT Design, Data Preparation, and Analysis

This repo is used to manage the R and Stan code used to design the TakeUp, prepare survey data and conduct analysis.

Abstract (see [working paper](takeup_workingpaper.pdf)):

> Can social image concerns motivate adults to internalize health externalities? In collaboration with the Kenyan Government, we implement a new community program that offers free deworming treatment to 200,000 adults and emphasizes the public good aspect of deworming. Importantly, we randomize the introduction of two types of social signals in the form of colorful bracelets and ink applied to the thumb. The bracelets and ink allow adults to signal that they contributed to protecting their community from worms. Further, we exogenously vary the travel distance to treatment locations. To separate reputational utility from private consumption utility and social learning/salience, we combine experimental identification with a structural model's non-experimental identification. We find that (1) bracelets as signals increase deworming take-up by roughly 13 percent, net the private consumption and learning/salience effects; (2) there is no detectable effect for the ink signal, which we attribute to its private disutility which outweighs any reputational utility from signaling; (3) adults are highly sensitive to distance but signaling treatments see little change; (4) the rate-of-change in take-up in response to increased distance is negligibly small compared against the control.  Detailed survey data on first and second-order beliefs shed light on the underlying mechanism: signals reduce information asymmetries, and adults are more likely to think that others have information about their deworming decision.

Here are the main files:

* [rct-design-fieldwork/](rct-design-fieldwork/) contains all scripts used for study design (e.g., randomization) and assisting with fieldwork/study implementation.
* [instruments/](instruments/) contains all the survey instruments and their translations.
* [data/](data/) points to a **private** submodule with survey data. This data has not yet been publicly shared.
* [images/](images/) stores all images used in producing the working paper.
* [stan_models/](stan_models/) contains all the Stan models used for analysis.
    - [takeup_reduced.stan](stan_models/takeup_reduced.stan) Reduced form model.
    - [takeup_struct.stan](stan_models/takeup_struct.stan) Structural model.
* [multilvlr/](multilvlr/) points to a public submodule. This is only used for some utility functions.

* [prepare_analysis_data.R](prepare_analysis_data.R) this script is used to prepare the survey data for analysis.
* [run_stan_dist.R](run_stan_dist.R) is the script that runs the Bayesian analysis for the working paper's main results.
* [postprocess_dist_fit.R](postprocess_dist_fit.R) postprocesses the fit data produced by [run_stan_dist.R](run_stan_dist.R).
* [takeup_workingpaper.Rmd](takeup_workingpaper.Rmd) this is the R markdown file use to produce the working paper PDF file.
* [takeup_workingpaper.pdf](takeup_workingpaper.pdf) this is the latest version of the working paper.
* [takeup_pap.Rmd](takeup_pap.Rmd) this is the R markdown file used to generate the pre-analysis plan PDF.

