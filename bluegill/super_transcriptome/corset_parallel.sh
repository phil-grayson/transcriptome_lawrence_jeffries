#!/bin/bash

# Call this script with 2 fastq file2 for reads 1 and 2.
#
# e.g. sbatch [thisscript.sh](http://thisscript.sh) myfile.R1.fastq  myfile.R2.fastq
#
# myfile.R1.fastq is referenced by the variable $1
# myfile.R2.fastq is referenced by the variable $2

#SBATCH -J corset 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -t 0-06:00                # Runtime in D-HH:MM 
#SBATCH --mem=28G            # Memory pool for all cores 
#SBATCH -o Corset.%A.out       # File to which STDOUT will be written 
#SBATCH -e Corset.%A.err       # File to which STDERR will be written 
#SBATCH --account=jwg-245-aa

/home/pgrayson/programs/corset-1.09-linux64/corset -r true-stop $1
