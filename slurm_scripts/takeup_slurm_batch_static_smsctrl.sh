#!/bin/bash

#SBATCH --job-name=takeup_stan_static_smsctrl
#SBATCH --output=/global/scratch/knaguib/slurm-logs/param_static_smsctrl-%j.out
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=2-12:00

module load r 

Rscript ~/Code/dewormtheworld/takeup/run_stan.R static --num-chains=8 --num-iterations=2000 --adapt-delta=0.8 --max-treedepth=20 --include-name-matched --separate-private-value --no-private-value-interact --sms-control-only --output-name=param_static_smsctrl_$SLURM_JOB_ID --output-dir=/global/scratch/knaguib
