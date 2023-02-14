#PBS -N "takeup-all-run" 
#PBS -j oe
#PBS -V
#PBS -e $PBS_O_WORKDIR/logs/qsub_log.txt
#PBS -l procs=16,mem=40gb
cd $PBS_O_WORKDIR

LATEST_VERSION=73
CMDSTAN_ARGS="--cmdstanr --include-paths=stan_models"
VERSION=${1:-$LATEST_VERSION}
OUTPUT_ARGS="--output-path=$PBS_O_WORKDIR/data/stan_analysis_data"
POSTPROCESS_INOUT_ARGS="--input_path=$PBS_O_WORKDIR/data/stan_analysis_data --output-path=$PBS_O_WORKDIR/temp-data"
STAN_THREADS=3




mpirun -n 1 -machinefile $PBS_NODEFILE 
	# Structural
	Rscript \
		--no-save \
		--no-restore \
		--verbose \
		run_takeup.R takeup fit \
		--models=STRUCTURAL_LINEAR_U_SHOCKS_NO_REP \
		${CMDSTAN_ARGS} \
		${OUTPUT_ARGS} \
		--threads=1 \
		--update \
		--outputname=dist_fit${VERSION} \
		--num-mix-groups=1 \
		--sequential 
		> logs/outputfile-structural-fit.Rout 2>&1 :
	Rscript \
		--no-save \
		--no-restore \
		--verbose \
		run_takeup.R takeup prior \
		--models=STRUCTURAL_LINEAR_U_SHOCKS_NO_REP \
		${CMDSTAN_ARGS} \
		${OUTPUT_ARGS} \
		--threads=1 \
		--update \
		--outputname=dist_prior${VERSION} \
		--num-mix-groups=1 \
		--sequential 
		> logs/outputfile-structural-prior.Rout 2>&1 :
	# Reduced Form
	Rscript \
		--no-save \
		--no-restore \
		--verbose \
		run_takeup.R takeup fit \
		--models=REDUCED_FORM_NO_RESTRICT \
		${CMDSTAN_ARGS} \
		${OUTPUT_ARGS} \
		--threads=1 \
		--update \
		--outputname=dist_fit${VERSION} \
		--num-mix-groups=1 \
		--sequential 
		> logs/outputfile-reduced-fit.Rout 2>&1 :
	Rscript \
		--no-save \
		--no-restore \
		--verbose \
		run_takeup.R takeup prior \
		--models=REDUCED_FORM_NO_RESTRICT \
		${CMDSTAN_ARGS} \
		${OUTPUT_ARGS} \
		--threads=1 \
		--update \
		--outputname=dist_prior${VERSION} \
		--num-mix-groups=1 \
		--sequential 
		> logs/outputfile-reduced-prior.Rout 2>&1 :
	# Postprocessing
	Rscript \
		--no-save \
		--no-restore \
		--verbose \
		postprocess_dist_fit.R \
		${VERSION} \
		${POSTPROCESS_INOUT_ARGS} \ 
		--cores=1 \ 
		--load-from-csv
