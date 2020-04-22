#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 1000
#SBATCH -t 0-12:00
#SBATCH -J blastp
#SBATCH -o blastp_%j.out
#SBATCH -e blastp_%j.err
#SBATCH --account=def-coling_cpu

module load nixpkgs/16.09 gcc/7.3.0 blast+/2.10.0

blastp -query longest_orfs.pep.vol.${SLURM_ARRAY_TASK_ID}.fasta -db ../../trinotate_dbs/uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.vol.${SLURM_ARRAY_TASK_ID}.outfmt6
