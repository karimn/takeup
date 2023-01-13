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
RUNESTIMATION="--run-estimation"


set -e

if [ ${POSTERIORMEDIAN} == "--posterior-median" ] 
then 
    POSTVAR="median"
else
    POSTVAR="post-draws"
fi

run_no_cutoff_optim () {
    if [ $SKIPPREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    $1 \
                                    $2 \
                                    --output-name=no-cutoff-b-$1-mu-$2 \
                                    --from-csv \
                                    --num-post-draws=${NUMPOSTDRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=10000 \
                                    --num-cores=12 \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --model=${MODEL} \
                                    ${PREDDISTANCE} \
                                    ${RUNESTIMATION} 
    fi

    Rscript ./optim/optimal_allocation.R  \
                                ${POSTERIORMEDIAN} \
                                --num-cores=12 \
                                --min-cost  \
                                --target-constraint=${TARGETCONSTRAINT} \
                                --output-path=optim/data \
                                --output-filename=no-cutoff-b-$1-mu-$2 \
                                --input-path=optim/data  \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --time-limit=10000 \
                                --demand-input-filename=pred-demand-dist-fit${VERSION}-no-cutoff-b-$1-mu-$2.csv

    Rscript ./optim/postprocess_allocation.R  \
                                --min-cost \
                                ${POSTERIORMEDIAN} \
                                --target-constraint=$TARGETCONSTRAINT \
                                --input-path=optim/data \
                                --optim-input-a-filename=no-cutoff-b-$1-mu-$2-${POSTVAR}-optimal-allocation.rds \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --demand-input-a-filename=pred-demand-dist-fit${VERSION}-no-cutoff-b-$1-mu-$2.csv \
                                --output-path=optim/plots \
                                --output-basename=no-cutoff-b-$1-mu-$2-${POSTVAR} \
                                --map-plot \
                                --cutoff-type=no-cutoff

}


run_cutoff_optim () {
    if [ $SKIPPREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    $1 \
                                    $2 \
                                    --output-name=cutoff-b-$1-mu-$2 \
                                    --from-csv \
                                    --num-post-draws=${NUMPOSTDRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=2500 \
                                    --num-cores=12 \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --model=${MODEL} \
                                    ${PREDDISTANCE} \
                                    ${RUNESTIMATION} 
    fi

    Rscript ./optim/optimal_allocation.R  \
                                ${POSTERIORMEDIAN} \
                                --num-cores=12 \
                                --min-cost  \
                                --target-constraint=${TARGETCONSTRAINT} \
                                --output-path=optim/data \
                                --output-filename=cutoff-b-$1-mu-$2 \
                                --input-path=optim/data  \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --time-limit=10000 \
                                --demand-input-filename=pred-demand-dist-fit${VERSION}-cutoff-b-$1-mu-$2.csv

    Rscript ./optim/postprocess_allocation.R  \
                                --min-cost \
                                ${POSTERIORMEDIAN} \
                                --target-constraint=$TARGETCONSTRAINT \
                                --input-path=optim/data \
                                --optim-input-a-filename=cutoff-b-$1-mu-$2-${POSTVAR}-optimal-allocation.rds \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --demand-input-a-filename=pred-demand-dist-fit${VERSION}-cutoff-b-$1-mu-$2.csv \
                                --output-path=optim/plots \
                                --output-basename=cutoff-b-$1-mu-$2-${POSTVAR} \
                                --map-plot \
                                --cutoff-type=cutoff

}



run_cutoff_optim "control" "control"
run_cutoff_optim "control" "bracelet"
run_cutoff_optim "bracelet" "bracelet"

run_no_cutoff_optim "control" "control"
run_no_cutoff_optim "control" "bracelet"
run_no_cutoff_optim "bracelet" "bracelet"
