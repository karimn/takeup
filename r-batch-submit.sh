#!/bin/bash
#PBS -N "takeup-all-run-800" 
#PBS -j oe
#PBS -V
#PBS -e /home/edjee/projects/takeup/logs/error.txt
#PBS -o /home/edjee/projects/takeup/logs/output.txt
#PBS -l procs=20,mem=100gb
cd $PBS_O_WORKDIR




fit_models () {

AKARING_OUTPUT_PATH="/share/akaringlab/takeup-output/"

LATEST_VERSION=90
CMDSTAN_ARGS="--cmdstanr --include-paths=stan_models"
VERSION=${1:-$LATEST_VERSION}
#OUTPUT_ARGS="--output-path=$PBS_O_WORKDIR/data/stan_analysis_data"
OUTPUT_ARGS="--output-path=${AKARING_OUTPUT_PATH}/data/stan_analysis_data"
POSTPROCESS_INOUT_ARGS="--input-path=${AKARING_OUTPUT_PATH}/data/stan_analysis_data --output-path=${AKARING_OUTPUT_PATH}/temp-data"
STAN_THREADS=4
ITER=800

#!/bin/bash

# Define an array of strings
fruits=("apple" "banana" "cherry")

# Loop through the array of strings
for fruit in "${fruits[@]}"
do
  echo "I like $fruit"
done

	models=(
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP" 
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_SD_WTP_VAL"
		"STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIGH_MU_WTP_VAL"
		"REDUCED_FORM_NO_RESTRICT_DIST_CTS")

	for model in "${models[@]}"
	do
		Rscript --no-save \
			--no-restore \
			--verbose \
			run_takeup.R takeup fit \
			--models=${model} \
			${CMDSTAN_ARGS} \
			${OUTPUT_ARGS} \
			--update \
			--threads=${STAN_THREADS} \
			--outputname=dist_fit${VERSION} \
			--num-mix-groups=1 \
			--iter=${ITER} \
			--sequential > logs/output-${model}-fit${VERSION}.txt 2>&1

		Rscript --no-save \
			--no-restore \
			--verbose \
			run_takeup.R takeup prior \
			--models=${model} \
			${CMDSTAN_ARGS} \
			${OUTPUT_ARGS} \
			--update \
			--threads=${STAN_THREADS} \
			--outputname=dist_prior${VERSION} \
			--num-mix-groups=1 \
			--iter=${ITER} \
			--sequential > logs/output-${model}-prior${VERSION}.txt 2>&1
	done
#
#		Rscript --no-save \
#			--no-restore \
#			--verbose \
#			run_takeup.R takeup fit \
#			--models=REDUCED_FORM_NO_RESTRICT \
#			${CMDSTAN_ARGS} \
#			${OUTPUT_ARGS} \
#			--update \
#			--threads=${STAN_THREADS} \
#			--outputname=dist_fit${VERSION} \
#			--num-mix-groups=1 \
#			--iter=${ITER} \
#			--sequential > output-rf-fit${VERSION}.txt 2>&1 
#
#		Rscript --no-save \
#			--no-restore \
#			--verbose \
#			run_takeup.R takeup prior \
#			--models=REDUCED_FORM_NO_RESTRICT \
#			${CMDSTAN_ARGS} \
#			${OUTPUT_ARGS} \
#			--update \
#			--threads=${STAN_THREADS} \
#			--outputname=dist_prior${VERSION} \
#			--num-mix-groups=1 \
#			--iter=${ITER} \
#			--sequential > output-rf-prior${VERSION}.txt 2>&1 
	


		Rscript --no-save \
			--no-restore \
			--verbose \
			postprocess_dist_fit.R \
			${VERSION} \
			--load-from-csv \
			--cores=1  \
			${POSTPROCESS_INOUT_ARGS} > logs/postprocess-output${VERSION}.txt 2>&1 

}

export -f fit_models

# call mpi
mpirun -n 1 -machinefile $PBS_NODEFILE  bash -c 'fit_models'
	

#LATEST_VERSION=88
#CMDSTAN_ARGS="--cmdstanr --include-paths=stan_models"
#VERSION=${1:-$LATEST_VERSION}
#OUTPUT_ARGS="--output-path=$PBS_O_WORKDIR/data/stan_analysis_data"
#POSTPROCESS_INOUT_ARGS="--input-path=$PBS_O_WORKDIR/data/stan_analysis_data --output-path=$PBS_O_WORKDIR/temp-data"
#STAN_THREADS=4
#ITER=800
#MODEL="STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
#
#
##mpirun -n 1 -machinefile $PBS_NODEFILE Rscript --no-save \
##			--no-restore \
##			--verbose \
##			postprocess_dist_fit.R \
##			${VERSION} \
##			--load-from-csv \
##			--cores=1  \
##			${POSTPROCESS_INOUT_ARGS}  
#
#
#
#
#mpirun -n 1 -machinefile $PBS_NODEFILE Rscript --no-save \
#			--no-restore \
#			--verbose \
#			postprocess_misc_params.R \
#			${VERSION} \
#			--load-from-csv \
#			--cores=1  \
#			--models=${MODEL} \
#			${POSTPROCESS_INOUT_ARGS} 
