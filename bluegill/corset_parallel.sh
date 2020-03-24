#!/bin/bash

#SBATCH -J corset 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -t 0-06:00                # Runtime in D-HH:MM 
#SBATCH --mem=28G            # Memory pool for all cores 
#SBATCH -o Corset.%A.out       # File to which STDOUT will be written 
#SBATCH -e Corset.%A.err       # File to which STDERR will be written 
#SBATCH --account=jwg-245-aa

# This was run on Grex

/home/pgrayson/programs/corset-1.09-linux64/corset -r true-stop $1
