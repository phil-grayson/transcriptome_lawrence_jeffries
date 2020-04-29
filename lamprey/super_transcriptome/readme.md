# superTranscriptome assembly using Corset and Lace
- Everything was run across Grex and Compute Canada Systems
- Corset requires that reads are mapped back to the all locations within the transcriptome

- Read mapping was carried out in Bowtie2
	- `sbatch bowtie2_makeindex.sh Trinity.fasta lamprey_RF_609610` <br/>
	- `ls trim*R1.fastq.gz > files.txt` <br/>
	- `for line in $(cat files.txt); do FILE2=$(echo $line | awk -F '[_]' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_R2.fastq.gz"}'); sbatch bowtie2_8_long.sh $line $FILE2 lamprey_RF_609610; sleep 0.1; done`

- Following the read mapping, Corset was run independently on each bam to generate a `corset-reads` file
	- `for FILE in $(ls ../*bam); do sbatch corset_parallel.sh $FILE; sleep 0.1; done`

- `corset-reads` files were used in two separate Corset runs.
	- subset used only those 22 libraries that had been chosen for the initial transcriptome assembly 
	- `sbatch corset_bigMem_subset.sh` 
	- full used all 148 libraries and required 1 month of runtime at 500 Gb of memory <br/>
	- `sbatch corset_full_bigMem_lamprey.sh`


- Corset output was then run through Lace to generate the superTranscriptome and annotation file
	- `sbatch lace.sh lamprey_lace_update`

