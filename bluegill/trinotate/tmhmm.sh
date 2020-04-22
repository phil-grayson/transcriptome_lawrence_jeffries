#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 500
#SBATCH -t 0-12:00
#SBATCH -J tmhmm
#SBATCH -o tmhmm_%j.out
#SBATCH -e tmhmm_%j.err
#SBATCH --account=def-coling_cpu

~/programs/tmhmm-2.0c/bin/tmhmm --short < longest_orfs.pep > tmhmm.out
