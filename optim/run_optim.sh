#!/bin/bash

# Setting arguments
PREDDISTANCE="" # --pred-distance
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP"
NUMPOSTDRAWS=200
SUBSIDY=0.25
TARGETCONSTRAINT=0.32
DRYRUN="" # "--sub-sample", "--fake-data"
VERSION=71
POSTERIORMEDIAN="--posterior-median" # --posterior-median
SKIPPREDICTION=0 # 1


set -e

for VARIABLE in "no-rep" "rep"
do
    if [ $VARIABLE == "rep" ]
    then
        REPVAR=""
    else
        REPVAR="--suppress-rep"
    fi

    if [ ${POSTERIORMEDIAN} == "--posterior-median" ] 
    then 
        POSTVAR="median"
    else
        POSTVAR="post-draws"
    fi
    
    if [ $SKIPPREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    --output-name=no-cutoff-${VARIABLE} \
                                    --from-csv \
                                    --num-post-draws=${NUMPOSTDRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=10000 \
                                    --num-cores=12 \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --model=${MODEL} \
                                    ${PREDDISTANCE} \
                                    ${REPVAR} \
                                    bracelet \
                                    calendar \
                                    control
    fi


    Rscript ./optim/optimal_allocation.R  \
                                ${POSTERIORMEDIAN} \
                                --num-cores=12 \
                                --min-cost  \
                                --target-constraint=${TARGETCONSTRAINT} \
                                --output-path=optim/data \
                                --output-filename=no-cutoff-${VARIABLE} \
                                --input-path=optim/data  \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --time-limit=10000 \
                                --demand-input-filename=pred-demand-dist-fit${VERSION}-no-cutoff-${VARIABLE}.csv


    Rscript ./optim/postprocess_allocation.R  \
                                --min-cost \
                                ${POSTERIORMEDIAN} \
                                --target-constraint=$TARGETCONSTRAINT \
                                --input-path=optim/data \
                                --optim-input-a-filename=no-cutoff-${VARIABLE}-${POSTVAR}-optimal-allocation.rds \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --demand-input-a-filename=pred-demand-dist-fit${VERSION}-no-cutoff-${VARIABLE}.csv \
                                --output-path=optim/plots \
                                --output-basename=no-cutoff-${VARIABLE}-${POSTVAR} \
                                --map-plot \
                                --cutoff-type=no-cutoff

    ## now cutoff
    if [ $SKIPPREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    --output-name=cutoff-${VARIABLE} \
                                    --from-csv \
                                    --num-post-draws=${NUMPOSTDRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=2500 \
                                    --num-cores=12 \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --model="STRUCTURAL_LINEAR_U_SHOCKS" \
                                    ${PREDDISTANCE} \
                                    ${REPVAR} \
                                    bracelet \
                                    calendar \
                                    control
    fi


    Rscript ./optim/optimal_allocation.R  \
                                ${POSTERIORMEDIAN} \
                                --num-cores=12 \
                                --min-cost  \
                                --target-constraint=${TARGETCONSTRAINT} \
                                --output-path=optim/data \
                                --output-filename=cutoff-${VARIABLE} \
                                --input-path=optim/data  \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --time-limit=10000 \
                                --demand-input-filename=pred-demand-dist-fit${VERSION}-cutoff-${VARIABLE}.csv 


    Rscript ./optim/postprocess_allocation.R  \
                                --min-cost \
                                ${POSTERIORMEDIAN} \
                                --target-constraint=$TARGETCONSTRAINT \
                                --input-path=optim/data \
                                --optim-input-a-filename=cutoff-${VARIABLE}-${POSTVAR}-optimal-allocation.rds \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --demand-input-a-filename=pred-demand-dist-fit${VERSION}-cutoff-${VARIABLE}.csv \
                                --output-path=optim/plots \
                                --output-basename=cutoff-${VARIABLE}-${POSTVAR} \
                                --map-plot \
                                --cutoff-type=cutoff
done


