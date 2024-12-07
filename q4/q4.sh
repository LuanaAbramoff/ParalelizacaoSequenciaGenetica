#!/bin/bash
#SBATCH --job-name=ex4
#SBATCH --output=ex4.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=2048MB
#SBATCH --time=00:05:00

mpirun -np 4 ./ex4 rna_transcrito.txt