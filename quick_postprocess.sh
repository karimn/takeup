#!/usr/bin/env bash

#SBATCH --partition=bigmem2
#SBATCH --job-name=quick-takeup        # create a short name for your job
#SBATCH --nodes=2                # node count
#SBATCH --ntasks=2              # total number of tasks across all nodes
#SBATCH --cpus-per-task=8      # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=10G         # memory per cpu-core (4G is default)
#SBATCH --time=0-10:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-%j.log
#SBATCH --error=temp/log/takeup-%j.log
#SBATCH --export=IN_SLURM=1

LATEST_VERSION=95
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided
CMDSTAN_ARGS="--cmdstanr"
SLURM_INOUT_DIR=~/scratch-midway2

echo "Version: $VERSION"

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."
	

  module load -f midway2 gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0

  POSTPROCESS_INOUT_ARGS="--input-path=${SLURM_INOUT_DIR} --output-path=${SLURM_INOUT_DIR}"

  echo "Running with ${CORES} cores."
  echo "INOUT ARGS: ${POSTPROCESS_INOUT_ARGS}."
else
  POSTPROCESS_INOUT_ARGS="--input-path=data/stan_analysis_data --output-path=temp-data/"
  CORES=8
fi

STAN_THREADS=$((${CORES} / 4))

	models=(
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB" 
		"REDUCED_FORM_NO_RESTRICT"
		)

	for model in "${models[@]}"
	do

postprocess_models () {
    echo "RUNNING: $1"
    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_ate_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            1 2 3 4

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_ate_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            --prior \
            1 2 3 4

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_roc_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            1 2 3 4

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_submodel_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            1 2 3 4
}


srun -n1 postprocess_models "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP" &
srun -n1 postprocess_models "REDUCED_FORM_NO_RESTRICT" &
wait



