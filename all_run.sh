#!/usr/bin/env bash

#SBATCH --partition=broadwl
#SBATCH --job-name=takeup        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=24       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=5G         # memory per cpu-core (4G is default)
#SBATCH --time=0-06:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/takeup-%j.log
#SBATCH --error=temp/log/takeup-%j.log
#SBATCH --export=IN_SLURM=1

LATEST_VERSION=91
VERSION=${1:-$LATEST_VERSION} # Get version from command line if provided
CMDSTAN_ARGS="--cmdstanr"
SLURM_INOUT_DIR=~/scratch-midway2
ITER=400
# MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"

echo "Version: $VERSION"

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."
	

  module load -f midway2 gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0

  OUTPUT_ARGS="--output-path=${SLURM_INOUT_DIR}"
  POSTPROCESS_INOUT_ARGS="--input-path=${SLURM_INOUT_DIR} --output-path=${SLURM_INOUT_DIR}"
  CORES=$SLURM_CPUS_PER_TASK

  echo "Running with ${CORES} cores."
  echo "INOUT ARGS: ${POSTPROCESS_INOUT_ARGS}."
else
  OUTPUT_ARGS="--output-path=data/stan_analysis_data"
  POSTPROCESS_INOUT_ARGS=
  CORES=8
fi

STAN_THREADS=$((${CORES} / 4))

fit_models () {
	models=(
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP" 
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_SD_WTP_VAL"
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_MU_WTP_VAL"
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_WTP_SUBMODEL"
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_BELIEFS_SUBMODEL"
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_SUBMODELS"
		"REDUCED_FORM_NO_RESTRICT_DIST_CTS"
		"REDUCED_FORM_NO_RESTRICT"
		)

	for [model in "${models[@]}"]
	do
		Rscript --no-save \
			--no-restore \
			--verbose \
			run_takeup.R takeup fit \
			--models=${model} \
			${CMDSTAN_ARGS} \
			${OUTPUT_ARGS} \
			--update-output \
			--threads=${STAN_THREADS} \
			--outputname=dist_fit${VERSION} \
			--num-mix-groups=1 \
			--chains=2 \
			--iter=${ITER} \
			--sequential > temp/log/output-${model}-fit${VERSION}.txt 2>&1

		Rscript --no-save \
			--no-restore \
			--verbose \
			run_takeup.R takeup prior \
			--models=${model} \
			${CMDSTAN_ARGS} \
			${OUTPUT_ARGS} \
			--update-output \
			--threads=${STAN_THREADS} \
			--outputname=dist_prior${VERSION} \
			--chains=2 \
			--num-mix-groups=1 \
			--iter=${ITER} \
			--sequential > temp/log/output-${model}-prior${VERSION}.txt 2>&1
	done

		Rscript --no-save \
			--no-restore \
			--verbose \
			postprocess_dist_fit.R \
			${VERSION} \
			--load-from-csv \
			--cores=1  \
			${POSTPROCESS_INOUT_ARGS} > temp/log/postprocess-output${VERSION}.txt 2>&1 

}


fit_models
# WTP 
# Rscript ./run_takeup.R wtp prior --chains=4 --iter 2000 --outputname=wtp_model_prior --output-path=data/stan_analysis_data --multilevel &
# Rscript ./run_takeup.R wtp fit   --chains=4 --iter 2000 --outputname=wtp_model_fit   --output-path=data/stan_analysis_data --multilevel
# wait

# Distance
#Rscript ./run_takeup.R dist prior --chains=4 --iter 800 --outputname=dist_model_prior --output-path=data/stan_analysis_data --num-mix-groups=1 --multilevel &
#Rscript ./run_takeup.R dist fit   --chains=4 --iter 800 --outputname=dist_model_fit   --output-path=data/stan_analysis_data --num-mix-groups=1 --multilevel
#wait

# Beliefs
#Rscript ./run_takeup.R beliefs prior --chains=4 --iter 1000 --outputname=beliefs_prior --output-path=data/stan_analysis_data &
#Rscript ./run_takeup.R beliefs fit   --chains=4 --iter 1000 --outputname=beliefs       --output-path=data/stan_analysis_data --num-mix-groups=1
#wait
