#!/bin/bash

#SBATCH --job-name=oxDNA
#SBATCH --account=csd913

#SBATCH --partition=shared
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
source ~/miniconda3/etc/profile.d/conda.sh
conda activate

### run
