#!/bin/bash -l
#SBATCH --job-name=kmb
#SBATCH --output=kmb.txt
#SBATCH -p cpu-preempt 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:30:00 

module load conda/latest
conda activate mothitor

cd /work/pi_mrobson_smith_edu/mth364/
python kmb_nf1.py