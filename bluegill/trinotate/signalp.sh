#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 20000
#SBATCH -t 1-00:00
#SBATCH -J signalP
#SBATCH -o sP_%j.out
#SBATCH -e sP_%j.err
#SBATCH --account=def-coling

~/programs/signalp-4.1/signalp -f short -n signalp.out longest_orfs.pep
