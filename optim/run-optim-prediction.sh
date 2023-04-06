#!/bin/bash
PRED_DISTANCE="--pred-distance" # --pred-distance
SUPPRESS_REP="" # "--suppress-rep"
OUTPUT_NAME="rep"
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
NUM_POST_DRAWS=200
# RUN_ESTIMATION="--run-estimation"
SINGLE_CHAIN="--single-chain"
VERSION=86


set -e

run_model_components () {
    Rscript ./optim/predict-takeup-for-optim.R \
                                ${VERSION} \
                                $1 \
                                $2 \
                                --output-name=static-cutoff-b-$1-mu-$2-${MODEL} \
                                --to-csv \
                                --num-post-draws=200 \
                                --rep-cutoff=Inf \
                                --dist-cutoff=3500 \
                                --type-lb=-Inf \
                                --type-ub=Inf \
                                --num-cores=12 \
                                --model=${MODEL} \
                                ${PRED_DISTANCE} \
                                --data-input-name=full-experiment.rds \
                                --single-chain \
                                --static-signal-pm \
                                --static-signal-distance=0
}


run_model_components "control" "control"
run_model_components "control" "bracelet"
run_model_components "control" "ink"
run_model_components "control" "calendar"

# run_model_components "bracelet" "bracelet"
# run_model_components "calendar" "calendar"
# run_model_components "ink" "ink"