#!/bin/bash
#SBATCH -J fastqc
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 0-02:59                 # Runtime in D-HH:MM
#SBATCH --mem=400               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o fastqc.%A.out  # File to which STDOUT will be written
#SBATCH -e fastqc.%A.err  # File to which STDERR will be written
#SBATCH --account=def-docker

# This was run on Graham

module load fastqc

for FILE in $(cat ${1}); do fastqc -f fastq $FILE; done

