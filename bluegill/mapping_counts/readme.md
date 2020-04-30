# Scripts and workflow for STAR mapping and featureCounts
- Everything here was carried out on Graham and Cedar in April 2020

- STAR was used to generate an index
	- `sbatch star_index.sh`

- Trimmed reads were mapped back to the Lace superTranscriptome using STAR
	- `ls *R1.f*.gz > files.txt` <br/>
	- `for line in $(cat files.txt); do FILE2=$(echo $line | awk -F '[_]' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_R2.fastq.gz"}'); sbatch star_run.sh $line $FILE2 bluegill_lace_update/star_index/; sleep 0.3; done` <br/>

- STAR bam files were used in two runs of featureCounts
	1. exon based 
	- `sbatch featureCounts.sh`
	2. transcript based	
	- `sbatch featureCounts_transcript_level.sh`
