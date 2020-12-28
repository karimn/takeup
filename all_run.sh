#!/usr/bin/env bash

VERSION=30
CMDSTAN_ARGS="--cmdstanr --include-paths=~/Code/takeup/stan_models"
# MODELS="--models=1,3,4"
MODELS="--models=STRUCTURAL_LINEAR"

./run_stan_dist.R fit --outputname=dist_fit${VERSION} ${MODELS} ${CMDSTAN_ARGS} --update
# ./run_stan_dist.R cv --outputname=dist_kfold${VERSION} --folds=10 ${MODELS} ${CMDSTAN_ARGS} --update
# ./postprocess_dist_fit.R $VERSION 

# Simulation
./run_stan_dist_sim.R --cmdstanr --include-paths=~/Code/takeup/stan_models