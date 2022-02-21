#!/bin/bash
# JOB HEADERS HERE
#SBATCH --job-name=cellranger612
#SBATCH --nodes=1
#SBATCH --cpus-per-task=164
#SBATCH --mem=1024G
#SBATCH --time=1-00:00:00

echo "singularity"
module load Singularity

echo "Start"

time singularity exec ./cellranger612.sif bash /script.sh

