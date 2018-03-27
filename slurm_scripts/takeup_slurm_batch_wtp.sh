#!/bin/bash

#SBATCH --job-name=takeup_stan_wtp
#SBATCH --output=/global/scratch/knaguib/slurm-logs/param_wtp-%j.out
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00

module load r 

Rscript ~/Code/dewormtheworld/takeup/run_stan.R wtp --num-chains=8 --num-iterations=300 --output-name=param_wtp_$SLURM_JOB_ID --output-dir=/global/scratch/knaguib 
