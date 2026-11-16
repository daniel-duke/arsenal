#!/bin/bash

#SBATCH --job-name=oxDNA

#SBATCH --partition=scavenger
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out

### modules
module purge
module load slurm

### python
source /hpc/group/aryalab/dhd13/miniconda3/etc/profile.d/conda.sh
conda activate

### run
