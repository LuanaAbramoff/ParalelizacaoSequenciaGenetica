#!/bin/bash
#SBATCH --job-name=ex2
#SBATCH --output=ex2.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048MB
#SBATCH --time=00:05:00

mpirun -np 4 ./ex2