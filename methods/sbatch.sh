#!/bin/bash
#SBATCH --job-name=wolbachia
#SBATCH --array=1-50%50

srun julia runPGG.jl
