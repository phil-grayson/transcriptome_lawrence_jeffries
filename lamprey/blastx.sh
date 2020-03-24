#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 5000
#SBATCH -t 3-00:00
#SBATCH -J trinotate
#SBATCH -o trinotate_%j.out
#SBATCH -e trinotate_%j.err
#SBATCH --account=def-kmj477 

# This was run on Beluga

module load nixpkgs/16.09 
module load gcc/7.3.0
module load trinotate/3.2.0

blastx -query Trinity.fasta.vol.${SLURM_ARRAY_TASK_ID}.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.vol.${SLURM_ARRAY_TASK_ID}.outfmt6
