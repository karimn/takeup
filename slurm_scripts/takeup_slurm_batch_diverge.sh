#!/bin/bash

#SBATCH --job-name=takeup_stan
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-12:00

module load r 

Rscript ~/Code/dewormtheworld/takeup/run_stan.R --num-chains=1 --num-iterations=300 --adapt-delta=0.999 --output-name=diver_2
