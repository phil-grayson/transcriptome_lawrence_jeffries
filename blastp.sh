#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 3000
#SBATCH -t 1-00:00
#SBATCH -J trinotate
#SBATCH -o trinotate_%j.out
#SBATCH -e trinotate_%j.err
#SBATCH --account=jwg-245-aa

/home/pgrayson/programs/ncbi-blast-2.10.0+/bin/blastp -query longest_orfs.pep.vol.${SLURM_ARRAY_TASK_ID}.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.vol.${SLURM_ARRAY_TASK_ID}.outfmt6
