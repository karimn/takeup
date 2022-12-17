#!/bin/bash

Rscript ./optim/predict-takeup-for-optim.R \
                            71 \
                            --output-name=rep \
                            --to-csv \
                            --num-post-draws=10 \
                            --rep-cutoff=Inf \
                            --dist-cutoff=10000 \
                            --num-cores=12 \
                            --type-lb=-3 \
                            --type-ub=3 \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_LINEAR_MU_REP \
                            --pred-distance \
                            bracelet \
                            control


Rscript ./optim/predict-takeup-for-optim.R \
                            71 \
                            --output-name=rep \
                            --to-csv \
                            --num-post-draws=10 \
                            --rep-cutoff=Inf \
                            --dist-cutoff=10000 \
                            --num-cores=12 \
                            --type-lb=-3 \
                            --type-ub=3 \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP \
                            --pred-distance \
                            bracelet \
                            control

Rscript ./optim/predict-takeup-for-optim.R \
                            71 \
                            --output-name=rep \
                            --to-csv \
                            --num-post-draws=10 \
                            --rep-cutoff=Inf \
                            --dist-cutoff=10000 \
                            --num-cores=12 \
                            --type-lb=-3 \
                            --type-ub=3 \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS \
                            --pred-distance \
                            bracelet \
                            control