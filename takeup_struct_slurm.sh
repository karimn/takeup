#!/bin/bash

#SBATCH --job-name=takeup_sim    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=32       # cpu-cores per task (>1 if multi-threaded tasks)
#BSATCH --mem=30G                # Total memory                     
#SBATCH --time=1-00:00:00        # maximum time needed (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=karimn2.0@gmail.com
#SBATCH --output=temp/log/takeup_struct-%j.log
#SBATCH --error=temp/log/takeup_struct-%j.log
#SBATCH --export=IN_SLURM=1

if [ -z ${IN_SLURM} ]; then
  module purge
  module load rh/devtoolset/8 gdal
fi

VERSION=42
CMDSTAN_ARGS="--cmdstanr --include-paths=~/Code/takeup/stan_models"
MODELS="--models=STRUCTURAL_LINEAR_U_SHOCKS" # REDUCED_FORM_NO_RESTRICT
SLURM_INOUT_DIR=/tigress/kn6838/takeup
NUM_STAN_THREADS=8

Rscript run_takeup.R takeup prior --sequential --chains=4 --threads=${NUM_STAN_THREADS} --outputname=dist_prior${VERSION} ${MODELS} ${CMDSTAN_ARGS} --output-path=${SLURM_INOUT_DIR} --update 
Rscript run_takeup.R takeup fit   --sequential --chains=4 --threads=${NUM_STAN_THREADS} --outputname=dist_fit${VERSION}   ${MODELS} ${CMDSTAN_ARGS} --output-path=${SLURM_INOUT_DIR} --update 

Rscript postprocess_dist_fit.R ${VERSION} --cores=32 --input-path=${SLURM_INOUT_DIR} --output-path=${SLURM_INOUT_DIR} --load-from-csv 

