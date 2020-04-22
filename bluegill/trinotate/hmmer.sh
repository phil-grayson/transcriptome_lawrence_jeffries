#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem 5000
#SBATCH -t 3-00:00
#SBATCH -J hmmer
#SBATCH -o hmmer_%j.out
#SBATCH -e hmmer_%j.err
#SBATCH --account=def-coling

module load nixpkgs/16.09 
module load gcc/7.3.0
module load trinotate/3.2.0

hmmscan --cpu 16 --domtblout TrinotatePFAM.out ../../trinotate_dbs/Pfam-A.hmm longest_orfs.pep > pfam.log
