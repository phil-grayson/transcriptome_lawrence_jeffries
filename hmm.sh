#!/bin/bash
#SBATCH -n 16
#SBATCH --mem 25G
#SBATCH -t 3-00:00
#SBATCH -J hmm
#SBATCH -o hmm_%j.out
#SBATCH -e hmm_%j.err
#SBATCH --account=wgf-720-aa


hmmscan --cpu 16 --domtblout TrinotatePFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log
