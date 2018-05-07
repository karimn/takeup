#!/bin/bash

#SBATCH --job-name=takeup_stan_know_table
#SBATCH --output=/global/scratch/knaguib/slurm-logs/param_know_table-%j.out
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00
#SBATCH --chdir=..
#SBATCH --mail-type=all
#SBATCH --mail-user=karim.naguib@evidenceaction.org

module load r 

Rscript ~/Code/dewormtheworld/takeup/run_stan.R beliefs --num-chains=8 --num-iterations=600 --adapt-delta=0.99 --output-name=param_beliefs_$SLURM_JOB_ID --output-dir=/global/scratch/knaguib 
