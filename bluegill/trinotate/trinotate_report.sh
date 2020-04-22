#!/bin/bash
#SBATCH -N 1
#SBATCH --mem 4G
#SBATCH -t 1-00:00
#SBATCH -J trinityReport
#SBATCH -o TrinityRep_%j.out
#SBATCH -e TrinityRep_%j.err
#SBATCH --account=def-docker_cpu

~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite report > trinotate_annotation_report_bluegill_RF.xls
