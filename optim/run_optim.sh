#!/usr/bin/env bash

#SBATCH --partition=broadwl
#SBATCH --job-name=optim        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=24       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=6G         # memory per cpu-core (4G is default)
#SBATCH --time=0-06:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-%j.log
#SBATCH --error=temp/log/takeup-%j.log
#SBATCH --export=IN_SLURM=1

LATEST_VERSION=71
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."

  module load midway2 gdal/2.4.1 udunits cmake R/4.2.0 gurobi/952

  NUM_CORES=$SLURM_CPUS_PER_TASK

  echo "Running with ${CORES} cores."
else
  NUM_CORES=16
fi

# Setting arguments
PRED_DISTANCE="" # --pred-distance
MODEL="STRUCTURAL_LINEAR_U_SHOCKS"
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
OUTPUT_PATH="optim/data/${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${COUNTY}/many-pots" # /many-pots
PLOT_OUTPUT_PATH="optim/plots/${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${COUNTY}/many-pots" #/many-pots
DATA_INPUT_NAME="${COUNTY}-experiment.rds"
CUTOFF="no-" # either no- or empty string
SOLVER="gurobi"
MANY_POTS="--many-pots" #"--many-pots"
SUPPRESS_REP="suppress-rep-" #suppress-rep-


mkdir -p ${OUTPUT_PATH}
mkdir -p ${PLOT_OUTPUT_PATH}

set -e



Rscript ./optim/create-distance-data.R \
    --output-name=${DATA_INPUT_NAME} \
    --num-extra-pots=100 \
    --county-subset=${COUNTY} \
    --distance-cutoff=3500

if [ ${POSTERIOR_MEDIAN} == "--posterior-median" ] 
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
        CUTOFF_DIST=3500
    fi

    if [ ${SUPPRESS_REP} != "" ]
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
                                    --from-csv \
                                    --num-post-draws=${NUM_POST_DRAWS} \
                                    --rep-cutoff=Inf \
                                    --dist-cutoff=${CUTOFF_DIST} \
                                    --num-cores=${NUM_CORES} \
                                    --type-lb=-Inf \
                                    --type-ub=Inf \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${OUTPUT_PATH} \
                                    --model=${MODEL} \
                                    ${PRED_DISTANCE} \
                                    ${RUN_ESTIMATION} \
                                    ${SUP_REP_VAR}
    fi

    if [ ${RUN_TARGET_CREATION} == 1 ]
    then

        Rscript ./optim/create-village-target.R \
            pred-demand-dist-fit${VERSION}-${CUTOFF}cutoff-b-control-mu-control-${MODEL}.csv \
            --input-path=${OUTPUT_PATH} \
            --output-path=${OUTPUT_PATH} \
            --num-cores=${NUM_CORES} \
            --output-basename=target-${CUTOFF}cutoff-b-control-mu-control-${MODEL} 
    fi

    if [ $SKIP_OA != 1 ]
    then
        Rscript ./optim/optimal_allocation.R  \
                                    ${POSTERIOR_MEDIAN} \
                                    --num-cores=12 \
                                    --min-cost  \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --target-constraint=target-${CUTOFF}cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS.csv \
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
        Rscript ./optim/postprocess_allocation.R  \
                                    --min-cost \
                                    ${POSTERIOR_MEDIAN} \
                                    --constraint-type=${CONSTRAINT_TYPE} \
                                    --welfare-function=${WELFARE_FUNCTION} \
                                    --optim-input-path=${OUTPUT_PATH} \
                                    --demand-input-path=${OUTPUT_PATH} \
                                    --optim-input-a-filename=${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR}-optimal-allocation.rds \
                                    --data-input-name=$DATA_INPUT_NAME \
                                    --output-path=${PLOT_OUTPUT_PATH} \
                                    --output-basename=${CONSTRAINT_TYPE}-${WELFARE_FUNCTION}-${SUPPRESS_REP}${CUTOFF}cutoff-b-$1-mu-$2-${MODEL}-${POSTVAR} \
                                    --cutoff-type=${CUTOFF}cutoff
    fi

}


CUTOFF=""
## Cutoff
run_optim "control" "control"
# run_optim "control" "bracelet"
# run_optim "control" "calendar"
# run_optim "control" "ink"
#
#
run_optim "bracelet" "bracelet"
run_optim "ink" "ink"
run_optim "calendar" "calendar"
#
# run_optim "bracelet" "control"
# run_optim "ink" "control"
# #

#  run_optim "calendar" "control"
#

CUTOFF="no-"

# ## No Cutoff
# run_optim "control" "control"
# #run_optim "control" "bracelet"
#Jrun_optim "control" "calendar"
#run_optim "control" "ink"
# run_optim "bracelet" "bracelet"
#

#run_optim "bracelet" "control"
# run_optim "ink" "ink"
# run_optim "calendar" "calendar"



Rscript ./optim/compare-optim.R --input-path=$OUTPUT_PATH --output-path=$PLOT_OUTPUT_PATH $MANY_POTS
