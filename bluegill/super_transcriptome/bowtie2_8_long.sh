#!/bin/bash

# Call this script with 2 fastq file2 for reads 1 and 2.
#
# e.g. sbatch [thisscript.sh](http://thisscript.sh) myfile.R1.fastq  myfile.R2.fastq
#
# myfile.R1.fastq is referenced by the variable $1
# myfile.R2.fastq is referenced by the variable $2

#SBATCH -J Bowtie2 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 8                # Use n cores for one job 
#SBATCH -t 3-00:00                # Runtime in D-HH:MM 
#SBATCH --mem=28000            # Memory pool for all cores 
#SBATCH -o bt2.%A.out       # File to which STDOUT will be written 
#SBATCH -e bt2.%A.err       # File to which STDERR will be written 
#SBATCH --account=wgf-720-aa

# This was run on Grex

module load gcc/5.2
module load intel/14.0.2.144
module load samtools/1.9

/home/pgrayson/programs/bowtie2-2.3.5.1-linux-x86_64/bowtie2 -x $3 -1 $1 -2 $2 -X 2000 -p 8 --all | samtools view -b -S - |samtools sort - -o $1.bam

samtools index $1.bam
