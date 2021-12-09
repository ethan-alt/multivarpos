#!/bin/bash

#SBATCH --job-name=compute_phat
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3g
#SBATCH	--array=1-105
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/phat/Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/phat/Error/%a.err

## add R module
module add r/3.6.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/R/03_compute_phat_v2.R /nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/phat/Rout/phat_$SLURM_ARRAY_TASK_ID.Rout