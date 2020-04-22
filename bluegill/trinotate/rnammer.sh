#!/bin/bash
#SBATCH -N 1
#SBATCH --mem 125G
#SBATCH --ntasks-per-node=32
#SBATCH -t 1-00:00
#SBATCH -J rnammer
#SBATCH -o Rnammer_%j.out
#SBATCH -e Rnammer_%j.err
#SBATCH --account=def-coling_cpu

/home/pgrayson/programs/Trinotate-Trinotate-v3.2.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /home/pgrayson/programs/RNAMMER/rnammer
