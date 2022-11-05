
#PBS -N "takeup-sbc" 
#PBS -j oe
#PBS -V
#PBS -e /home/edjee/projects/takeup/logs/qsub_log.txt
#PBS -l procs=36,mem=100gb
cd $PBS_O_WORKDIR

CMDSTAN_ARGS=" --include-paths=stan_models"
VERSION=${1:-$LATEST_VERSION}
OUTPUT_ARGS="--output-path=$PBS_O_WORKDIR/data/sbc_output_data"


mpirun -n 1 -machinefile $PBS_NODEFILE Rscript \
	--no-save \
	--no-restore \
	--verbose \
	takeup_reduced_sbc.R sbc  \
	--models=REDUCED_FORM_NO_RESTRICT \
	${CMDSTAN_ARGS} \
	${OUTPUT_ARGS} \
	--threads=1 \
	--num-cores=24 \
	--save-mem \
	--num-sbc-draws=500 \
	--outputname=initial_reduced_form_sbc${VERSION} \
	--multilevel > logs/outputfile.Rout 2>&1
