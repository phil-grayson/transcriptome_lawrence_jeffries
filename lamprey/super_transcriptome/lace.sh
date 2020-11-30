#!/bin/bash

#SBATCH -J lace 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH --ntasks-per-node=32
#SBATCH -t 7-00:00                # Runtime in D-HH:MM 
#SBATCH --mem=125G            # Memory pool for all cores 
#SBATCH -o Lace.%A.out       # File to which STDOUT will be written 
#SBATCH -e Lace.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-coling_cpu

source ~/ENV3/bin/activate
module load blat
python ~/programs/Lace-1.14.1/Lace/Lace_run.py -t --cores 32 --outputDir $1 Trinity.fasta clusters.txt 
