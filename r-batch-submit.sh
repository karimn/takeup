#PBS -N "takeup" 
#PBS -j oe
#PBS -V
#PBS -e $PBS_O_WORKDIR/logs/qsub_log.txt
#PBS -l procs=20,mem=10gb
cd $PBS_O_WORKDIR

CMDSTAN_ARGS="--cmdstanr --include-paths=stan_models"
VERSION=${1:-$LATEST_VERSION}
OUTPUT_ARGS="--output-path=$PBS_O_WORKDIR"


mpirun -n 1 -machinefile $PBS_NODEFILE Rscript \
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
	--multilevel > logs/outputfile.Rout 2>&1