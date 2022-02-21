#!/bin/bash
# JOB HEADERS HERE
#SBATCH --job-name=demuxlet
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=800G
#SBATCH --time=2-00:00:00
echo "singularity"
module load Singularity
echo "Start"
time singularity exec ./popscle.sif bash /script.sh
