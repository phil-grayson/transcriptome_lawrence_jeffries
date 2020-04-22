#!/bin/bash

#SBATCH -J corset 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -t 7-00:00                # Runtime in D-HH:MM 
#SBATCH --mem=50G            # Memory pool for all cores 
#SBATCH -o Corset.%A.out       # File to which STDOUT will be written 
#SBATCH -e Corset.%A.err       # File to which STDERR will be written 
#SBATCH --account=def-docker_cpu

module load corset

corset -D 99999999999 -i corset trim_NS.1271.001.NEBNext_dual_i7_E2---NEBNext_dual_i5_E2.LC12-6G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_131---NEBNext_dual_i5_131.LC12-1L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_143---NEBNext_dual_i5_143.LC24-23L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_B4---NEBNext_dual_i5_B4.LC24-19G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_129---NEBNext_dual_i5_129.LC6-23L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_H1---NEBNext_dual_i5_H1.LC6-21G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_97---NEBNext_dual_i5_97.LN12-3G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_F12---NEBNext_dual_i5_F12.LN12-7L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_104---NEBNext_dual_i5_104.LN24-4G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_172---NEBNext_dual_i5_172.LN24-2L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_F7---NEBNext_dual_i5_F7.LN6-1G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_165---NEBNext_dual_i5_165.LN6-7L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_C6---NEBNext_dual_i5_C6.LT12-7G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_152---NEBNext_dual_i5_152.LT12-6L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_F6---NEBNext_dual_i5_F6.LT24-1G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_161---NEBNext_dual_i5_161.LT24-7L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_C5---NEBNext_dual_i5_C5.LT6-8G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_147---NEBNext_dual_i5_147.LT6-6L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_112---NEBNext_dual_i5_112.LTN12-2G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_187---NEBNext_dual_i5_187.LTN12-22L_R1.fastq.gz.bam.corset-reads,trim_NS.1271.001.NEBNext_dual_i7_110---NEBNext_dual_i5_110.LTN6-9G_R1.fastq.gz.bam.corset-reads,trim_NS.1271.002.NEBNext_dual_i7_179---NEBNext_dual_i5_179.LTN6-8L_R1.fastq.gz.bam.corset-reads -l 5 -x 1000
