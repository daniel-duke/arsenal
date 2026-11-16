#!/bin/bash

#SBATCH --job-name=oxDNA
#SBATCH --account=csd913

#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --gpus=1

#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out

### modules
module purge
module load slurm

### python
source ~/miniconda3/etc/profile.d/conda.sh
conda activate

### run
bash run_stages.sh
