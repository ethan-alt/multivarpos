#!/bin/bash

#SBATCH --job-name=phat
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH	--array=1-330
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/PoS/Redo/phat/Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/PoS/Redo/phat/Error/%a.err

## add R module
module add r/3.6.0
R CMD BATCH --no-restore ./R/02_compute_phat_v2.R ~/Projects/PoS/Redo/phat/Rout/sample_$SLURM_ARRAY_TASK_ID.Rout