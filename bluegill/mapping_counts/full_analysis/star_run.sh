#!/bin/bash

#SBATCH -J star_run 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH --ntasks-per-node=8
#SBATCH -t 3-00:00                # Runtime in D-HH:MM 
#SBATCH --mem=20G            # Memory pool for all cores 
#SBATCH -o Star_run.%A.out       # File to which STDOUT will be written 
#SBATCH -e Star_run.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-coling_cpu

#module load star
#module load nixpkgs/16.09  intel/2018.3 star/2.7.3a
module load nixpkgs/16.09  intel/2016.4 star/2.6.1a

STAR --runMode alignReads --twopassMode Basic --readFilesCommand zcat --runThreadN 8 --readFilesIn $1 $2 --outSAMtype BAM SortedByCoordinate --genomeDir $3 --outFileNamePrefix star_map/mapped_${1} --sjdbOverhang 99
