#!/bin/bash
#SBATCH -p phi
#SBATCH -o mat.out
#SBATCH --mem-per-cpu=1GB
#SBATCH -t 120:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

#module load lapack/3.7.0/b1

time ./mat_1d.exe