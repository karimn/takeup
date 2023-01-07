#!/bin/bash
PREDDISTANCE="" # --pred-distance
SUPPRESSREP="" # "--suppress-rep"
OUTPUTNAME="rep"
MODEL="STRUCTURAL_LINEAR_U_SHOCKS"
NUMPOSTDRAWS=200
RUNESTIMATION="--run-estimation"


set -e

Rscript ./optim/predict-takeup-for-optim.R \
                            71 \
                            --output-name=${OUTPUTNAME} \
                            --from-csv \
                            --num-post-draws=${NUMPOSTDRAWS} \
                            --rep-cutoff=Inf \
                            --dist-cutoff=10000 \
                            --num-cores=12 \
                            --type-lb=-Inf \
                            --type-ub=Inf \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_LINEAR_MU_REP \
                            ${PREDDISTANCE} \
                            ${SUPPRESSREP} \
                            ${RUNESTIMATION} \
                            bracelet \
                            control


Rscript ./optim/predict-takeup-for-optim.R \
                            71 \
                            --output-name=${OUTPUTNAME} \
                            --from-csv \
                            --num-post-draws=${NUMPOSTDRAWS} \
                            --rep-cutoff=Inf \
                            --dist-cutoff=10000 \
                            --num-cores=12 \
                            --type-lb=-Inf \
                            --type-ub=Inf \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP \
                            ${PREDDISTANCE} \
                            ${SUPPRESSREP} \
                            ${RUNESTIMATION} \
                            bracelet \
                            control

Rscript ./optim/predict-takeup-for-optim.R \
                            71 \
                            --output-name=${OUTPUTNAME} \
                            --from-csv \
                            --num-post-draws=${NUMPOSTDRAWS} \
                            --rep-cutoff=Inf \
                            --dist-cutoff=10000 \
                            --num-cores=12 \
                            --type-lb=-Inf \
                            --type-ub=Inf \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS \
                            ${PREDDISTANCE} \
                            ${SUPPRESSREP} \
                            ${RUNESTIMATION} \
                            bracelet \
                            control