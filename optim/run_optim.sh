#!/bin/bash

# Setting arguments
PREDDISTANCE="" # --pred-distance
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP"
NUMPOSTDRAWS=200
TARGET_CONSTRAINT="target-cutoff-b-control-mu-control-${MODEL}.csv"
DRYRUN="" # "--sub-sample", "--fake-data"
VERSION=71
POSTERIORMEDIAN="--posterior-median" # --posterior-median
SKIPPREDICTION=0 # 1
SKIP_OA=0 # 1 or 0
SKIP_PP=0 # 1 or 0
RUNESTIMATION="--run-estimation"
WELFARE_FUNCTION="log"
CONSTRAINT_TYPE="agg"
OUTPUT_PATH="optim/data/${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}"
PLOT_OUTPUT_PATH="optim/plots/agg-log"
RUN_TARGET_CREATION=1


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
                                    --output-name=no-cutoff-b-$1-mu-$2-${MODEL} \
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

    if [ ${RUN_TARGET_CREATION} == 1 ]
    then

        Rscript ./optim/create-village-target.R \
            pred-demand-dist-fit${VERSION}-no-cutoff-b-control-mu-control-${MODEL}.csv \
            --input-path=optim/data \
            --output-path=optim/data \
            --output-basename=target-no-cutoff-b-control-mu-control-${MODEL} 
    fi

    if [ $SKIP_OA != 1 ]
    then
        Rscript ./optim/optimal_allocation.R  \
                                    ${POSTERIORMEDIAN} \
                                    --num-cores=12 \
                                    --min-cost  \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --target-constraint=${TARGET_CONSTRAINT} \
                                    --output-path=${OUTPUT_PATH} \
                                    --output-filename=no-cutoff-b-$1-mu-$2-${MODEL} \
                                    --input-path=optim/data  \
                                    --village-input-filename=village-df.csv \
                                    --pot-input-filename=pot-df.csv \
                                    --time-limit=1000 \
                                    --demand-input-filename=pred-demand-dist-fit${VERSION}-no-cutoff-b-$1-mu-$2-${MODEL}.csv \
                                    --welfare-function=${WELFARE_FUNCTION}
    fi

    if [ $SKIP_PP != 1 ]
    then
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIORMEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --demand-input-path="optim/data" \
                                    --optim-input-a-filename=no-cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --village-input-filename=village-df.csv \
                                    --pot-input-filename=pot-df.csv \
                                    --demand-input-a-filename=pred-demand-dist-fit${VERSION}-no-cutoff-b-$1-mu-$2-${MODEL}.csv \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-no-cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=no-cutoff
    fi

}


run_cutoff_optim () {
    if [ $SKIPPREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    $1 \
                                    $2 \
                                    --output-name=cutoff-b-$1-mu-$2-${MODEL} \
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

    if [ ${RUN_TARGET_CREATION} == 1 ]
    then
        Rscript ./optim/create-village-target.R \
            pred-demand-dist-fit${VERSION}-cutoff-b-control-mu-control-${MODEL}.csv \
            --input-path=optim/data \
            --output-path=optim/data \
            --output-basename=target-cutoff-b-control-mu-control-${MODEL} 
    fi
    if [ $SKIP_OA != 1 ]
    then
        Rscript ./optim/optimal_allocation.R  \
                                    ${POSTERIORMEDIAN} \
                                    --num-cores=12 \
                                    --min-cost  \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --target-constraint=${TARGET_CONSTRAINT} \
                                    --output-path=${OUTPUT_PATH} \
                                    --output-filename=cutoff-b-$1-mu-$2-${MODEL} \
                                    --input-path=optim/data  \
                                    --village-input-filename=village-df.csv \
                                    --pot-input-filename=pot-df.csv \
                                    --time-limit=1000 \
                                    --demand-input-filename=pred-demand-dist-fit${VERSION}-cutoff-b-$1-mu-$2-${MODEL}.csv \
                                    --welfare-function=${WELFARE_FUNCTION}
    fi

    if [ $SKIP_PP != 1 ]
    then
    Rscript ./optim/postprocess_allocation.R  \
                                --min-cost \
                                ${POSTERIORMEDIAN} \
                                --constraint-type=${CONSTRAINT_TYPE} \
                                --welfare-function=${WELFARE_FUNCTION} \
                                --optim-input-path=${OUTPUT_PATH} \
                                --demand-input-path="optim/data" \
                                --optim-input-a-filename=cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                --village-input-filename=village-df.csv \
                                --pot-input-filename=pot-df.csv \
                                --demand-input-a-filename=pred-demand-dist-fit${VERSION}-cutoff-b-$1-mu-$2-${MODEL}.csv \
                                --output-path=${PLOT_OUTPUT_PATH} \
                                --output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                --cutoff-type=cutoff
    fi

}


#export -f run_cutoff_optim
#export -f run_no_cutoff_optim
#
#parallel --link run_cutoff_optim \
#    {1} {2} \
#    ::: control control bracelet \
#    ::: control bracelet bracelet

## Cutoff
run_cutoff_optim "control" "control"
run_cutoff_optim "control" "bracelet"

run_cutoff_optim "control" "calendar"

# run_cutoff_optim "bracelet" "bracelet"


# ## No Cutoff
run_no_cutoff_optim "control" "control"
run_no_cutoff_optim "control" "bracelet"

run_no_cutoff_optim "control" "calendar"

# run_no_cutoff_optim "bracelet" "bracelet"
