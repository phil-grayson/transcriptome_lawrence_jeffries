#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 5000
#SBATCH -t 3-00:00
#SBATCH -J tmhmm
#SBATCH -o tmhmm_%j.out
#SBATCH -e tmhmm_%j.err
#SBATCH --account=jwg-245-aa

~/programs/tmhmm-2.0c/bin/tmhmm --short < longest_orfs.pep > tmhmm.out
