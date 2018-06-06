#!/bin/bash

#SBATCH --job-name=takeup_stan_static_census_covar
#SBATCH --output=/global/scratch/knaguib/slurm-logs/param_static_census_covar-%j.out
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=1-00:00
#SBATCH --chdir=..
#SBATCH --mail-user=karim.naguib@evidenceaction.org
#SBATCH --mail-type=all 

module load r 

Rscript ~/Code/dewormtheworld/takeup/run_stan.R static --num-chains=8 --num-iterations=600 --adapt-delta=0.9 --max-treedepth=25 --use-census-covar --output-name=param_static_census_covar_$SLURM_JOB_ID --output-dir=/global/scratch/knaguib #  --include-name-matched  
