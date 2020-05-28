#!/bin/bash

#SBATCH -J featureCounts
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH --ntasks-per-node=8
#SBATCH -t 0-03:00                # Runtime in D-HH:MM 
#SBATCH --mem=8G            # Memory pool for all cores 
#SBATCH -o FeatureCounts.%A.out       # File to which STDOUT will be written 
#SBATCH -e FeatureCounts.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-kmj477_cpu

module load nixpkgs/16.09 gcc/7.3.0 subread/2.0.0

mkdir Counts_transcripts
featureCounts -T 8 -p -s 2 -B -O --fraction -a ../lamprey_lace_full/SuperDuper.gff -o Counts_transcripts/Lamprey_subsetSuperTrans_counts.txt *.out.bam
