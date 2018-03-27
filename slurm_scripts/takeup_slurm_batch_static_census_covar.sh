#!/bin/bash

#SBATCH --job-name=takeup_stan_static_census_covar
#SBATCH --output=/global/scratch/knaguib/slurm-logs/param_static_census_covar-%j.out
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00
#SBATCH --chdir=.. 

module load r 

Rscript ~/Code/dewormtheworld/takeup/run_stan.R static --num-chains=4 --num-iterations=300 --use-census-covar --output-name=param_static_census_covar_$SLURM_JOB_ID --output-dir=/global/scratch/knaguib 
