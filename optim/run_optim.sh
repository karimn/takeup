#!/bin/bash

# Setting arguments
PREDDISTANCE="" # --pred-distance
MODEL="STRUCTURAL_LINEAR_U_SHOCKS"
NUMPOSTDRAWS=50
#TARGET_CONSTRAINT="target-cutoff-b-control-mu-control-${MODEL}.csv"
DRYRUN="" # "--sub-sample", "--fake-data"
VERSION=71
POSTERIORMEDIAN="--posterior-median" # --posterior-median
SKIPPREDICTION=1 # 1
SKIP_OA=0 # 1 or 0
SKIP_PP=0 # 1 or 0
RUN_TARGET_CREATION=0
RUNESTIMATION="--run-estimation"
WELFARE_FUNCTION="log"
CONSTRAINT_TYPE="agg"
OUTPUT_PATH="optim/data/${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-kakamega"
PLOT_OUTPUT_PATH="optim/plots/${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-kakamega"
NUM_CORES=6
DATA_INPUT_NAME="KAKAMEGA-experiment.rds"
CUTOFF="no-" # either no- or empty string
SOLVER="gurobi"

mkdir -p ${OUTPUT_PATH}
mkdir -p ${PLOT_OUTPUT_PATH}

set -e


# Rscript ./optim/create-distance-data.R \
#     --output-name=${DATA_INPUT_NAME} \
#     --num-extra-pots=4 \
#     --county=KAKAMEGA

if [ ${POSTERIORMEDIAN} == "--posterior-median" ] 
then 
    POSTVAR="median"
else
    POSTVAR="post-draws"
fi


run_optim () {

    if [ ${CUTOFF} == "no-" ]
    then
        CUTOFF_DIST=10000
    else
        CUTOFF_DIST=2500
    fi

    # Gurobi doesn't play nice with renv atm
    if [ ${SOLVER} == "gurobi" ]
    then
        VANILLA="--vanilla" 
    else
        VANILLA=""
    fi
    echo $VANILLA
    echo $SOLVER

    if [ $SKIPPREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    $1 \
                                    $2 \
                                    --output-name=${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --from-csv \
                                    --num-post-draws=${NUMPOSTDRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=${CUTOFF_DIST} \
									--num-cores=${NUM_CORES} \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${OUTPUT_PATH} \
                                    --model=${MODEL} \
                                    ${PREDDISTANCE} \
                                    ${RUNESTIMATION} 
    fi

    if [ ${RUN_TARGET_CREATION} == 1 ]
    then

        Rscript ./optim/create-village-target.R \
            pred-demand-dist-fit${VERSION}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv \
            --input-path=${OUTPUT_PATH} \
            --output-path=${OUTPUT_PATH} \
            --output-basename=target-${CUTOFF}cutoff-b-control-mu-control-${MODEL} 
    fi

    if [ $SKIP_OA != 1 ]
    then
        Rscript ${VANILL} ./optim/optimal_allocation.R  \
                                    ${POSTERIORMEDIAN} \
                                    --num-cores=12 \
                                    --min-cost  \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --target-constraint=target-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv \
                                    --output-path=${OUTPUT_PATH} \
                                    --output-filename=${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --input-path=${OUTPUT_PATH}  \
                                    --data-input-path=optim/data \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --time-limit=1000 \
                                    --demand-input-filename=pred-demand-dist-fit${VERSION}-${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}.csv \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --solver=${SOLVER}
    fi

    if [ $SKIP_PP != 1 ]
    then
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIORMEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --demand-input-path=${OUTPUT_PATH} \
                                    --optim-input-a-filename=${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff
    fi

}


# run_cutoff_optim () {
#     if [ $SKIPPREDICTION != 1 ]
#     then

#         Rscript ./optim/predict-takeup-for-optim.R \
#                                     ${VERSION} \
#                                     $1 \
#                                     $2 \
#                                     --output-name=cutoff-b-$1-mu-$2-${MODEL} \
#                                     --from-csv \
#                                     --num-post-draws=${NUMPOSTDRAWS} \
#                                     --rep-cutoff=Inf \
#                                     --dist-cutoff=10000 \
# 									--num-cores=${NUM_CORES} \
#                                     --type-lb=-Inf \
#                                     --type-ub=Inf \
#                                     --data-input-name=$DATA_INPUT_NAME \
#                                     --output-path=${OUTPUT_PATH} \
#                                     --model=${MODEL} \
#                                     ${PREDDISTANCE} \
#                                     ${RUNESTIMATION} 
#         Rscript ./optim/predict-takeup-for-optim.R \
#                                     ${VERSION} \
#                                     $1 \
#                                     $2 \
#                                     --output-name=cutoff-b-$1-mu-$2-${MODEL} \
#                                     --from-csv \
#                                     --num-post-draws=${NUMPOSTDRAWS} \
#                                     --rep-cutoff=Inf \
#                                     --dist-cutoff=2500 \
# 									--num-cores=${NUM_CORES} \
#                                     --type-lb=-Inf \
#                                     --type-ub=Inf \
#                                     --data-input-name=$DATA_INPUT_NAME \
#                                     --model=${MODEL} \
#                                     --output-path=${OUTPUT_PATH} \
#                                     ${PREDDISTANCE} \
#                                     ${RUNESTIMATION} 
#     fi

#     if [ ${RUN_TARGET_CREATION} == 1 ]
#     then
#         Rscript ./optim/create-village-target.R \
#             pred-demand-dist-fit${VERSION}-cutoff-b-control-mu-control-${MODEL}.csv \
#             --input-path=optim/data \
#             --output-path=optim/data \
#             --output-basename=target-cutoff-b-control-mu-control-${MODEL} 
#     fi
#     if [ $SKIP_OA != 1 ]
#     then
#         Rscript ./optim/optimal_allocation.R  \
#                                     ${POSTERIORMEDIAN} \
#                                     --num-cores=12 \
#                                     --min-cost  \
#                                     --constraint-type=${CONSTRAINT_TYPE} \
#                                     --target-constraint=target-cutoff-b-control-mu-control-${MODEL} \
#                                     --output-path=${OUTPUT_PATH} \
#                                     --output-filename=cutoff-b-$1-mu-$2-${MODEL} \
#                                     --input-path=optim/data  \
#                                     --data-input-name=$DATA_INPUT_NAME \
#                                     --time-limit=1000 \
#                                     --demand-input-filename=pred-demand-dist-fit${VERSION}-cutoff-b-$1-mu-$2-${MODEL}.csv \
#                                     --welfare-function=${WELFARE_FUNCTION}
#     fi

#     if [ $SKIP_PP != 1 ]
#     then
#     Rscript ./optim/postprocess_allocation.R  \
#                                 --min-cost \
#                                 ${POSTERIORMEDIAN} \
#                                 --constraint-type=${CONSTRAINT_TYPE} \
#                                 --welfare-function=${WELFARE_FUNCTION} \
#                                 --optim-input-path=${OUTPUT_PATH} \
#                                 --demand-input-path="optim/data" \
#                                 --optim-input-a-filename=cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
#                                 --data-input-name=$DATA_INPUT_NAME \
#                                 --demand-input-a-filename=pred-demand-dist-fit${VERSION}-cutoff-b-$1-mu-$2-${MODEL}.csv \
#                                 --output-path=${PLOT_OUTPUT_PATH} \
#                                 --output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
#                                 --cutoff-type=cutoff
#     fi

# }


#export -f run_cutoff_optim
#export -f run_no_cutoff_optim
#
#parallel --link run_cutoff_optim \
#    {1} {2} \
#    ::: control control bracelet \
#    ::: control bracelet bracelet


CUTOFF=""
## Cutoff
run_optim "control" "control"
run_optim "control" "bracelet"

run_optim "control" "calendar"

# run_optim "bracelet" "bracelet"


CUTOFF="no-"
# ## No Cutoff
# run_optim "control" "control"
# run_optim "control" "bracelet"

# run_optim "control" "calendar"

# run_optim "bracelet" "bracelet"
