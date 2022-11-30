#PBS -N "takeup" 
#PBS -j oe
#PBS -V
#PBS -e $PBS_O_WORKDIR/logs/qsub_log.txt
#PBS -l procs=30,mem=60gb
cd $PBS_O_WORKDIR

LATEST_VERSION=70
CMDSTAN_ARGS="--cmdstanr --include-paths=stan_models"
VERSION=${1:-$LATEST_VERSION}
OUTPUT_ARGS="--output-path=$PBS_O_WORKDIR"


mpirun -n 1 -machinefile $PBS_NODEFILE Rscript \
	--no-save \
	--no-restore \
	--verbose \
	run_takeup.R takeup fit \
	--models=STRUCTURAL_LINEAR_U_SHOCKS \
	${CMDSTAN_ARGS} \
	${OUTPUT_ARGS} \
	--threads=1 \
	--update \
	--outputname=dist_fit${VERSION} \
	--num-mix-groups=1 \
	--sequential \
	--multilevel 
	> logs/outputfile.Rout 2>&1
