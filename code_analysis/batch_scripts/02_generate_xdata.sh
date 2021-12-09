#!/bin/bash

#SBATCH --job-name=x_vprior
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3g
#SBATCH	--array=1-12
#SBATCH --output=/nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/x_sample/Log/slurmLogFiles%a.out
#SBATCH --error=/nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/x_sample/Error/%a.err

## add R module
module add r/3.6.0

R CMD BATCH --no-restore /nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/R/02_generate_xdata.R /nas/longleaf/home/ethanalt/Projects/PoS/Redo_DataAnalysis/x_sample/Rout/xdata_$SLURM_ARRAY_TASK_ID.Rout