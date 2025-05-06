#!/bin/bash -l
#SBATCH --job-name=kma
#SBATCH --output=kma_idk.txt
#SBATCH -p cpu-preempt 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1:00:00 

module load conda/latest
conda activate mothitor

cd /work/pi_mrobson_smith_edu/mth364/
python kma_nf1.py