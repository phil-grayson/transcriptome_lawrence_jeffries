#!/bin/bash

#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 12000
#SBATCH -t 3-00:00
#SBATCH -J trinotate
#SBATCH -o trinotate_%j.out
#SBATCH -e trinotate_%j.err
#SBATCH --account=def-kmj477

module load nixpkgs/16.09 
module load gcc/7.3.0
module load trinotate/3.2.0

blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query longest_orfs.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log
module load signalp
signalp -f short -n signalp.out longest_orfs.pep
module load tmhmm
tmhmm --short < longest_orfs.pep > tmhmm.out
RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /usr/bin/software/rnammer_v1.2/rnammer
