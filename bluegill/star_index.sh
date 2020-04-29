#!/bin/bash

#SBATCH -J star_index 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH --ntasks-per-node=32
#SBATCH -t 0-01:00                # Runtime in D-HH:MM 
#SBATCH --mem=125G            # Memory pool for all cores 
#SBATCH -o Star_index.%A.out       # File to which STDOUT will be written 
#SBATCH -e Star_index.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-docker_cpu

module load star

#First Index
starindex="star_index2"
mkdir $starindex

STAR --runMode genomeGenerate --runThreadN 32 --genomeDir $starindex --genomeFastaFiles SuperDuper.fasta --sjdbGTFfile SuperDuper.gff --sjdbGTFtagExonParentTranscript gene_id --sjdbOverhang 99 --limitGenomeGenerateRAM 120000000000 --limitIObufferSize 30000000000 --genomeChrBinNbits 10
