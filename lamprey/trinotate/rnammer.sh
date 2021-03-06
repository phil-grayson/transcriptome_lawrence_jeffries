#!/bin/bash
#SBATCH -N 1
#SBATCH --mem 502G
#SBATCH -t 0-12:00
#SBATCH -J rnammer
#SBATCH -o Rnammer_%j.out
#SBATCH -e Rnammer_%j.err
#SBATCH --account=def-kmj477

/home/pgrayson/programs/Trinotate-Trinotate-v3.2.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /home/pgrayson/programs/RNAMMER/rnammer
