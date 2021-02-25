#!/usr/bin/env bash

#SBATCH --job-name=takeup_postprocess    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=12       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=0-00:30:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=karimn2.0@gmail.com
#SBATCH --output=temp/log/takeup_postprocess-%j.log
#SBATCH --error=temp/log/takeup_postprocess-%j.log
#SBATCH --export=IN_SLURM=1


if [[ -v IN_SLURM ]]; then
  module purge
  module load rh/devtoolset/8 gdal
fi

VERSION=31
CMDSTAN_ARGS="--cmdstanr --include-paths=~/Code/takeup/stan_models"
MODELS="--models=STRUCTURAL_LINEAR_U_SHOCKS" # REDUCED_FORM_NO_RESTRICT
OUTPUT_ARGS="--outputname=dist_prior${VERSION} --output-path=~/Code/takeup/data/stan_analysis_data"

#./run_stan_dist.R prior ${MODELS} ${CMDSTAN_ARGS} ${OUTPUT_ARGS} --update
# ./run_stan_dist.R fit --outputname=dist_fit${VERSION} ${MODELS} ${CMDSTAN_ARGS} ${OUTPUT_ARGS} --update
# ./run_stan_dist.R cv --outputname=dist_kfold${VERSION} --folds=10 ${MODELS} ${CMDSTAN_ARGS} ${OUTPUT_ARGS} --update
Rscript ./postprocess_dist_fit.R ${VERSION} --cores=12

# Simulation
# ./run_stan_dist_sim.R ${CMDSTAN_ARGS} 