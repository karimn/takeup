#!/usr/bin/env bash

VERSION=31
CMDSTAN_ARGS="--cmdstanr --include-paths=~/Code/takeup/stan_models"
MODELS="--models=STRUCTURAL_LINEAR_U_SHOCKS" # REDUCED_FORM_NO_RESTRICT

./run_stan_dist.R prior --outputname=dist_prior${VERSION} ${MODELS} ${CMDSTAN_ARGS} --update
# ./run_stan_dist.R fit --outputname=dist_fit${VERSION} ${MODELS} ${CMDSTAN_ARGS} --update
# ./run_stan_dist.R cv --outputname=dist_kfold${VERSION} --folds=10 ${MODELS} ${CMDSTAN_ARGS} --update
# ./postprocess_dist_fit.R $VERSION 

# Simulation
# ./run_stan_dist_sim.R ${CMDSTAN_ARGS} 