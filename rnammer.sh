#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 28G
#SBATCH -t 1-00:00
#SBATCH -J rnammer
#SBATCH -o Rnammer_%j.out
#SBATCH -e Rnammer_%j.err
#SBATCH --account=jwg-245-aa

/home/pgrayson/programs/Trinotate-Trinotate-v3.2.0/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /home/pgrayson/programs/RNAMMER/rnammer
