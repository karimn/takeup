
#!/usr/bin/env bash

LATEST_VERSION=86
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided

NUM_CORES=16

# Setting arguments
PRED_DISTANCE="" # --pred-distance
MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
NUM_POST_DRAWS=200
POSTERIOR_MEDIAN="--posterior-median" # --posterior-median
SKIP_PREDICTION=1 # 1
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
CONSTRAINT_TARGET="rep"
STATIC_SIGNAL_PM="" # "--static-signal-pm"
STATIC_SIGNAL_DIST=
DEMAND_NAME="" # "static-"

mkdir -p ${OUTPUT_PATH}
mkdir -p ${PLOT_OUTPUT_PATH}

set -e



Rscript ./optim/create-distance-data.R \
    --output-name=${DATA_INPUT_NAME} \
    --num-extra-pots=100 \
    --county-subset=${COUNTY} \
    --distance-cutoff=3500


# Create experiment target
Rscript ./optim/create-experiment-target.R \
                --constraint-type=agg \
                --welfare-function=log \
                --min-cost \
                --output-path=${OUTPUT_PATH} \
                --output-basename=summ-agg-log \
                --cutoff-type=cutoff \
                --data-input-name=full-many-pots-experiment.rds \
                --posterior-median \
                --demand-input-path=${OUTPUT_PATH} \
                --demand-input-filename=pred-demand-dist-fit${VERSION}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv


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
                                    --output-name=${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --to-csv \
                                    --num-post-draws=${NUM_POST_DRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=${CUTOFF_DIST} \
                                    --num-cores=${NUM_CORES} \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --data-input-name=${DATA_INPUT_NAME} \
                                    --output-path=${OUTPUT_PATH} \
                                    --model=${MODEL} \
                                    --single-chain \
                                    ${STATIC_SIGNAL_PM} \
                                    --static-signal-distance=${STATIC_SIGNAL_DIST} \
                                    ${PRED_DISTANCE} \
                                    ${RUN_ESTIMATION} \
                                    ${SUP_REP_VAR}
    fi

    if [ ${RUN_TARGET_CREATION} == 1 ]
    then

        if [ ${CONSTRAINT_TARGET} == "closest"]
        then
            TMP_CONSTRAINT_VAR="--assign-closest"
        else
            TMP_CONSTRAINT_VAR=""
        fi

        echo "Running target creation"
            Rscript ./optim/create-village-target.R \
                pred-demand-dist-fit${VERSION}-${SUPPRESS_REP}${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv \
                --input-path=${OUTPUT_PATH} \
                --output-path=${OUTPUT_PATH} \
                --num-cores=${NUM_CORES} \
                --output-basename=target-${CONSTRAINT_TARGET}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}  \
                ${TMP_CONSTRAINT_VAR}


    fi

    if [ $SKIP_OA != 1 ]
    then
    echo "Running optimization"
        Rscript ./optim/optimal_allocation.R  \
                                    ${POSTERIOR_MEDIAN} \
                                    --num-cores=12 \
                                    --min-cost  \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --target-constraint=summ-agg-log-experiment-target-constraint.csv \
                                    --output-path=${OUTPUT_PATH} \
                                    --output-filename=target-${CONSTRAINT_TARGET}-${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL} \
                                    --input-path=${OUTPUT_PATH}  \
                                    --data-input-path=optim/data \
                                    --data-input-name=${DATA_INPUT_NAME} \
                                    --time-limit=10000 \
                                    --demand-input-filename=pred-demand-dist-fit${VERSION}-${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}.csv \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --solver=${SOLVER}
                                    # --target-constraint=target-${CONSTRAINT_TARGET}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv \

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
                                    --optim-input-a-filename=target-${CONSTRAINT_TARGET}-${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --data-input-name=${DATA_INPUT_NAME} \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-target-${CONSTRAINT_TARGET}-${WELFARE_FUNCTION}-${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff \
                                    --pdf-output-path=presentations/takeup-${MODEL}-fig
    fi

}



compare_option () {
    if [[ ${CONSTRAINT_TARGET} == "rep-" ]]
    then
        TMP_REP_VAR_A=""
        TMP_REP_VAR_B=""
    else 
        TMP_REP_VAR_A="suppress-rep-"
        TMP_REP_VAR_B=""
    fi
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIOR_MEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --optim-input-a-filename=target-${CONSTRAINT_TARGET}-${DEMAND_NAME}${TMP_REP_VAR_A}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --optim-input-b-filename=target-${CONSTRAINT_TARGET}-${DEMAND_NAME}${TMP_REP_VAR_B}${CUTOFF}cutoff-b-$3-mu-$4-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --comp-output-basename=${CONSTRAINT_TYPE}-target-${CONSTRAINT_TARGET}-${WELFARE_FUNCTION}-${SUPPRESS_REP}${CUTOFF}cutoff-b1-$1-mu1-$2-b2-$3-mu2-$4-${MODEL}-${POSTVAR} \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-target-${CONSTRAINT_TARGET}-${WELFARE_FUNCTION}-${DEMAND_NAME}${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff \
                                    --pdf-output-path=presentations/takeup-${MODEL}-fig/
}



run_optim "control" "control" # run control control
run_optim "control" "bracelet" # counterfactual varying bracelet visibility
run_optim "bracelet" "bracelet" # now bracelet bracelet
# now we suppress reputation completely. 
# this is because I didn't think of creating a treatment variable with 0 visibility
# so we change a global variable woooo
SUPPRESS_REP="suppress-rep-"
run_optim "bracelet" "bracelet"
# Now we swap to static signalling, fixed at d = 0.5
SUPPRESS_REP=""
STATIC_SIGNAL_PM="--static-signal-pm" # "--static-signal-pm"
STATIC_SIGNAL_DIST=500
DEMAND_NAME="static-" 
run_optim "bracelet" "bracelet"


Rscript ./optim/create-presentation-plots.R \
                            --constraint-type=agg \
                            --welfare-function=log \
                            --min-cost \
                            --output-path=${OUTPUT_PATH} \
                            --output-basename=${CONSTRAINT_TYPE}-target-${CONSTRAINT_TARGET}-${WELFARE_FUNCTION}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}-${POSTVAR} \
                            --cutoff-type=cutoff
                            --data-input-name=full-many-pots-experiment.rds \
                            --posterior-median \
                            --pdf-output-path=presentations/takeup-${MODEL}-fig
                            --demand-input-path=optim/data/${MODEL}/agg-log-full-many-pots \
                            --demand-input-filename=pred-demand-dist-fit${VERSION}-cutoff-b-control-mu-control-${MODEL}.csv

Rscript ./optim/create-presentation-plots.R \
                            --constraint-type=agg \
                            --welfare-function=log \
                            --min-cost \
                            --output-path=${OUTPUT_PATH} \
                            --output-basename=${CONSTRAINT_TYPE}-target-${CONSTRAINT_TARGET}-${WELFARE_FUNCTION}-suppress-rep-${CUTOFF}cutoff-b-bracelet-mu-bracelet-${MODEL}-${POSTVAR} \
                            --cutoff-type=cutoff
                            --data-input-name=full-many-pots-experiment.rds \
                            --posterior-median \
                            --pdf-output-path=presentations/takeup-${MODEL}-fig
                            --demand-input-path=optim/data/${MODEL}/agg-log-full-many-pots \
                            --demand-input-filename=pred-demand-dist-fit${VERSION}-suppress-rep-cutoff-b-bracelet-mu-bracelet-${MODEL}.csv