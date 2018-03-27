#!/bin/bash

#SBATCH --job-name=takeup_stan_dynamic
#SBATCH --output=/global/scratch/knaguib/slurm-logs/param_dynamic-%j.out
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=2-12:00

module load r 

Rscript ~/Code/dewormtheworld/takeup/run_stan.R dynamic --num-chains=8 --num-iterations=300 --include-name-matched --output-name=param_dynamic_$SLURM_JOB_ID --output-dir=/global/scratch/knaguib
