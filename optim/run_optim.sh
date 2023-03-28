#!/usr/bin/env bash

LATEST_VERSION=86
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided

NUM_CORES=16

# Setting arguments
PRED_DISTANCE="" # --pred-distance
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
NUM_POST_DRAWS=200
POSTERIOR_MEDIAN="--posterior-median" # --posterior-median
SKIP_PREDICTION=0 # 1
SKIP_OA=0 # 1 or 0
SKIP_PP=0 # 1 or 0
RUN_TARGET_CREATION=1
RUN_ESTIMATION="--run-estimation"
WELFARE_FUNCTION="log"
CONSTRAINT_TYPE="agg"
COUNTY="full"
OUTPUT_PATH="optim/data/${MODEL}/${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${COUNTY}-many-pots" # /many-pots
PLOT_OUTPUT_PATH="optim/plots/${MODEL}/${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${COUNTY}-many-pots" #/many-pots
DATA_INPUT_NAME="${COUNTY}-many-pots-experiment.rds"
CUTOFF="" # either no- or empty string
SOLVER="gurobi"
MANY_POTS="--many-pots" #"--many-pots"
SUPPRESS_REP="" # "suppress-rep-" #suppress-rep-


mkdir -p ${OUTPUT_PATH}
mkdir -p ${PLOT_OUTPUT_PATH}

set -e



Rscript ./optim/create-distance-data.R \
    --output-name=${DATA_INPUT_NAME} \
    --num-extra-pots=100 \
    --county-subset=${COUNTY} \
    --distance-cutoff=3500

if [[ ${POSTERIOR_MEDIAN} == "--posterior-median" ]]
then 
    POSTVAR="median"
else
    POSTVAR="post-draws"
fi


run_optim () {

    if [[ ${CUTOFF} == "no-" ]]
    then
        CUTOFF_DIST=10000
    else
        CUTOFF_DIST=3500
    fi

    if [[ ${SUPPRESS_REP} != "" ]]
    then
        SUP_REP_VAR="--suppress-reputation"
    else
        SUP_REP_VAR=""
    fi

    # Gurobi doesn't play nice with renv atm

    if [ $SKIP_PREDICTION != 1 ]
    then
        Rscript ./optim/predict-takeup-for-optim.R \
                                    ${VERSION} \
                                    $1 \
                                    $2 \
                                    --output-name=${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --to-csv \
                                    --num-post-draws=${NUM_POST_DRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=${CUTOFF_DIST} \
                                    --num-cores=${NUM_CORES} \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${OUTPUT_PATH} \
                                    --model=${MODEL} \
                                    --single-chain \
                                    ${PRED_DISTANCE} \
                                    ${RUN_ESTIMATION} \
                                    ${SUP_REP_VAR}
    fi

    if [ ${RUN_TARGET_CREATION} == 1 ]
    then
    echo "Running target creation"
        Rscript ./optim/create-village-target.R \
            pred-demand-dist-fit${VERSION}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv \
            --input-path=${OUTPUT_PATH} \
            --output-path=${OUTPUT_PATH} \
            --num-cores=${NUM_CORES} \
            --output-basename=target-${CUTOFF}cutoff-b-control-mu-control-${MODEL} 
    fi

    if [ $SKIP_OA != 1 ]
    then
    echo "Running optimization"
        Rscript ./optim/optimal_allocation.R  \
                                    ${POSTERIOR_MEDIAN} \
                                    --num-cores=12 \
                                    --min-cost  \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --target-constraint=target-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv \
                                    --output-path=${OUTPUT_PATH} \
                                    --output-filename=${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --input-path=${OUTPUT_PATH}  \
                                    --data-input-path=optim/data \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --time-limit=10000 \
                                    --demand-input-filename=pred-demand-dist-fit${VERSION}-${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}.csv \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --solver=${SOLVER}
    fi

    if [ $SKIP_PP != 1 ]
    then
        echo "Running postprocessing"
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIOR_MEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --optim-input-a-filename=${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff \
                                    --pdf-output-path=presentations/takeup-${MODEL}-fig
    fi

}

                            # --constraint-type=agg \
                            # --welfare-function=log \
                            # --min-cost \
                            # --optim-input-path=optim/data/agg-log-full-many-pots \
                            # --optim-input-a-filename=cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS-median-optimal-allocation.rds \
                            # --optim-input-b-filename=cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS-median-optimal-allocation.rds \
                            # --comp-output-basename=agg-log-cutoff-cf-b1-control-b2-control-mu1-control-mu2-bracelet-STRUCTURAL_LINEAR_U_SHOCKS-median \
                            # --output-path=optim/plots/agg-log-full-many-pots \
                            # --output-basename=agg-log-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS-median \
                            # --cutoff-type=cutoff
                            # --data-input-name=full-many-pots-experiment.rds \
                            # --posterior-median \
                            # --pdf-output-path=presentations/takeup-fig/optim


compare_option () {
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIOR_MEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --optim-input-a-filename=${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --optim-input-b-filename=${SUPPRESS_REP}${CUTOFF}cutoff-b-$3-mu-$4-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --comp-output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${SUPPRESS_REP}${CUTOFF}cutoff-b1-$1-mu1-$2-b2-$3-mu2-$4-${MODEL}-${POSTVAR} \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff \
                                    --pdf-output-path=presentations/takeup-fig/optim/
}

CUTOFF=""
## Cutoff
run_optim "control" "control"
run_optim "control" "bracelet"
run_optim "control" "calendar"
run_optim "control" "ink"
#
#
run_optim "bracelet" "bracelet"
run_optim "ink" "ink"
run_optim "calendar" "calendar"
# #
# run_optim "bracelet" "control"
# run_optim "ink" "control"
# #

#  run_optim "calendar" "control"
#



# this ed hhhhhhhhhhhhhhhhhhhh
# if [[ ${SUPPRESS_REP} == "" ]]
# then 
#     compare_option "control" "control" "control" "bracelet"
#     compare_option "control" "control" "control" "ink"
#     compare_option "control" "control" "control" "calendar"
# fi

# if [[ ${SUPPRESS_REP} == "suppress-rep-" ]]
# then
#     compare_option "bracelet" "bracelet" "bracelet" "bracelet"
#     compare_option "control" "control" "control" "control"
#     compare_option "ink" "ink" "ink" "ink"
#     compare_option "calendar" "calendar" "calendar" "calendar"
# fi

CUTOFF="no-"

# ## No Cutoff
# run_optim "control" "control"
# # run_optim "control" "bracelet"

# run_optim "control" "calendar"
# run_optim "control" "ink"
# run_optim "bracelet" "bracelet"


# run_optim "bracelet" "control"
# run_optim "ink" "ink"
# run_optim "calendar" "calendar"



# Rscript ./optim/compare-optim.R \ 
#     --input-path=${OUTPUT_PATH} \
#     --output-path=${PLOT_OUTPUT_PATH} \
#     ${MANY_POTS} \
#     ${POSTERIOR_MEDIAN}

