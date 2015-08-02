#!/bin/bash

#SBATCH -o /scratch/pr63so/ga25cux2/roughgen/script_output.%j.out
#SBATCH -D /scratch/pr63so/ga25cux2/
#SBATCH -J roughgen
#SBATCH --mail-type=END
#SBATCH --mail-user=johannes.klicpera@tum.de
#SBATCH --export=NONE
#SBATCH --time=13:00:00
#SBATCH --nodes=8
#SBATCH --partition=snb
#SBATCH --begin=now

source /etc/profile.d/modules.sh
module load python

export OMP_NUM_THREADS=16
export mpi_ranks=8

cd /home/hpc/pr63so/ga25cux2/roughgen

python -u ./parallel_gen.py

python -u ./start_sims.py
