# Initial QC and read trimming
- Carried out on Compute Canada systems in March of 2020

- Fastqc takes about 7 minutes per sample, so I split the full list of files into 20 batches so that they would fit under the 3 hour time limit
	- `ls *gz > files.txt`
	- `split files.txt -l 20`
	- `ls x* > files.txt`
	- `for list in $(cat files.txt); do sbatch fastqc.sh $list; sleep 0.2; done`
	- a few jobs required a little extra memory, so I generated a file list for those and submitted:
	- `sbatch fastqc_mem.sh makeup`

- I downloaded and installed MultiQC to collect all the important data from the 100's of resulting Fastqc files and then ran:
	- `multiqc .`

- To get trimmomatic to run, I needed a list of generic file names without R1 or R2.  I reduced my terminal width to only allow one file across and ran:
	- `ls NS*R1*gz`
	- I copied this list to BBEdit (Textwrangler) and used find/replace with grep to replace `1\.fastq\.gz` with nothing
	- I copied the resulting list back onto the cluster as `files.txt`
	- `for pair in $(cat files.txt); do sbatch trim.sh ${pair}*gz; sleep 0.1; done`

- Fastqc and Multiqc were then run again (as above) on the trimmed files

- The best sample per condition (as defined by highest read number while accounting for duplication rate) was selected for the trinity assembly
	- This resulted in 26 for bluegill and 22 for lamprey
	- These files are listed in the samples.txt files in the assembly folder and were also used for the corset run with the subset of samples

