#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 5000
#SBATCH -t 3-00:00
#SBATCH -J trinotate
#SBATCH -o trinotate_%j.out
#SBATCH -e trinotate_%j.err
#SBATCH --account=xjc-394-aa

/home/pgrayson/programs/ncbi-blast-2.10.0+/bin/blastx -query Trinity.fasta.vol.${SLURM_ARRAY_TASK_ID}.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.vol.${SLURM_ARRAY_TASK_ID}.outfmt6
