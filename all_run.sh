#!/usr/bin/env bash

#SBATCH --job-name=takeup_postprocess    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=12       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=10G        # memory per cpu-core (4G is default)
#SBATCH --time=0-01:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=karimn2.0@gmail.com
#SBATCH --output=temp/log/takeup_postprocess-%j.log
#SBATCH --error=temp/log/takeup_postprocess-%j.log
#SBATCH --export=IN_SLURM=1

VERSION=52
CMDSTAN_ARGS="--cmdstanr --include-paths=~/Code/takeup/stan_models"
MODELS="--models=STRUCTURAL_LINEAR_U_SHOCKS" # REDUCED_FORM_NO_RESTRICT
# MODELS="--models=REDUCED_FORM_NO_RESTRICT"
SLURM_INOUT_DIR=/tigress/kn6838/takeup

if [[ -v IN_SLURM ]]; then
  echo "Running in SLURM..."

  module purge
  module load rh/devtoolset/8 gdal

  OUTPUT_ARGS="--output-path=${SLURM_INOUT_DIR}"
  POSTPROCESS_INOUT_ARGS="--input-path=${SLURM_INOUT_DIR} --output-path=${SLURM_INOUT_DIR}"
  CORES=$SLURM_CPUS_PER_TASK

  echo "Running with ${CORES} cores."
  echo "INOUT ARGS: ${POSTPROCESS_INOUT_ARGS}."
else
  OUTPUT_ARGS="--output-path=~/Code/takeup/data/stan_analysis_data"
  POSTPROCESS_INOUT_ARGS=
  CORES=12
fi

# Distance
# Rscript ./run_takeup.R dist prior --chains=4 --iter 800 --outputname=dist_model_prior --output-path=~/Code/takeup/data/stan_analysis_data --include-paths=~/Code/takeup/stan_models --num-mix-groups=1 --multilevel &
# Rscript ./run_takeup.R dist fit   --chains=4 --iter 800 --outputname=dist_model_fit   --output-path=~/Code/takeup/data/stan_analysis_data --include-paths=~/Code/takeup/stan_models --num-mix-groups=1 --multilevel
# wait

# Beliefs
# Rscript ./run_takeup.R beliefs prior --chains=4 --iter 1000 --outputname=beliefs_prior --output-path=~/Code/takeup/data/stan_analysis_data --include-paths=~/Code/takeup/stan_models &
# Rscript ./run_takeup.R beliefs fit   --chains=4 --iter 1000 --outputname=beliefs       --output-path=~/Code/takeup/data/stan_analysis_data --include-paths=~/Code/takeup/stan_models --num-mix-groups=1
# wait

# Reduced Form 
# Rscript ./run_takeup.R takeup prior --models=REDUCED_FORM_NO_RESTRICT ${CMDSTAN_ARGS} ${OUTPUT_ARGS} --threads=3 --update --outputname=dist_prior${VERSION} --multilevel
# Rscript ./run_takeup.R takeup fit   --models=REDUCED_FORM_NO_RESTRICT ${CMDSTAN_ARGS} ${OUTPUT_ARGS} --threads=3 --update --outputname=dist_fit${VERSION}   --multilevel
# Rscript ./run_takeup.R takeup cv    --models=REDUCED_FORM_NO_RESTRICT ${CMDSTAN_ARGS} ${OUTPUT_ARGS}             --update --outputname=dist_kfold${VERSION} --folds=10

# Structural
Rscript ./run_takeup.R takeup prior --models=STRUCTURAL_LINEAR_U_SHOCKS ${CMDSTAN_ARGS} ${OUTPUT_ARGS} --threads=3 --update --outputname=dist_prior${VERSION} --num-mix-groups=1 --multilevel
# Rscript ./run_takeup.R takeup fit   --models=STRUCTURAL_LINEAR_U_SHOCKS ${CMDSTAN_ARGS} ${OUTPUT_ARGS} --threads=3 --update --outputname=dist_fit${VERSION}   --num-mix-groups=1 
# Rscript ./run_takeup.R takeup cv    --models=STRUCTURAL_LINEAR_U_SHOCKS ${CMDSTAN_ARGS} ${OUTPUT_ARGS}             --update --outputname=dist_kfold${VERSION} --num-mix-groups=1 --folds=10

Rscript ./postprocess_dist_fit.R ${VERSION} ${POSTPROCESS_INOUT_ARGS} --cores=$CORES --load-from-csv

# Simulation
# ./run_stan_dist_sim.R ${CMDSTAN_ARGS} 
