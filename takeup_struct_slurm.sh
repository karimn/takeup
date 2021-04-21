#!/bin/bash

#SBATCH --job-name=takeup_struct # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=32       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=50G                # Total memory                     
#SBATCH --time=0-12:00:00        # maximum time needed (HH:MM:SS)
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

VERSION=51
CMDSTAN_ARGS="--cmdstanr --include-paths=~/Code/takeup/stan_models"
# MODELS="--models=STRUCTURAL_LINEAR_U_SHOCKS,REDUCED_FORM_NO_RESTRICT"
SLURM_INOUT_DIR=/tigress/kn6838/takeup
NUM_STAN_THREADS=8
NUM_DIST_MIX=1

Rscript run_takeup.R takeup prior --models=STRUCTURAL_LINEAR_U_SHOCKS --sequential --chains=4 --threads=${NUM_STAN_THREADS} --num-mix-groups=${NUM_DIST_MIX} --outputname=dist_prior${VERSION} ${CMDSTAN_ARGS} --output-path=${SLURM_INOUT_DIR} --update 
Rscript run_takeup.R takeup prior --models=REDUCED_FORM_NO_RESTRICT   --sequential --chains=4 --threads=${NUM_STAN_THREADS} --multilevel                     --outputname=dist_prior${VERSION} ${CMDSTAN_ARGS} --output-path=${SLURM_INOUT_DIR} --update 
Rscript run_takeup.R takeup fit   --models=STRUCTURAL_LINEAR_U_SHOCKS --sequential --chains=4 --threads=${NUM_STAN_THREADS} --num-mix-groups=${NUM_DIST_MIX} --outputname=dist_fit${VERSION} ${CMDSTAN_ARGS} --output-path=${SLURM_INOUT_DIR} --update 
Rscript run_takeup.R takeup fit   --models=REDUCED_FORM_NO_RESTRICT   --sequential --chains=4 --threads=${NUM_STAN_THREADS} --multilevel                     --outputname=dist_fit${VERSION} ${CMDSTAN_ARGS} --output-path=${SLURM_INOUT_DIR} --update 

Rscript postprocess_dist_fit.R ${VERSION} --cores=16 --input-path=${SLURM_INOUT_DIR} --output-path=${SLURM_INOUT_DIR} --load-from-csv 

