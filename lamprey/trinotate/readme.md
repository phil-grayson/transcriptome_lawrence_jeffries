# Trinotate scripts and walkthrough
- Everything was carried out on Graham in April of 2020 with a locally installed Trinotate-Trinotate-v3.2.1 and modules <br/>
`wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.2.1.tar.gz` <br/>
`tar -xf Trinotate-v3.2.1.tar.gz`<br/>

- RF `Trinity.fasta` assemblies and `Trinity.fasta.gene_trans_map` files were scp'd over from biology-01 <br/>
- The trinotate database `Trinotate.sqlite` was initiated with: <br/>
`~/programs/Trinotate-Trinotate-v3.2.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate` <br/>
- This above command also generates `uniprot_sprot.dat` and `Pfam-A.hmm.gz` which are both used for downstream applications <br/>
- These files must be modified as follows to be used as input:<br/>
`module load nixpkgs/16.09 gcc/7.3.0 trinotate/3.2.0`<br/>
`makeblastdb -in uniprot_sprot.pep -dbtype prot`<br/>
`gunzip Pfam-A.hmm.gz`<br/>
`hmmpress Pfam-A.hmm`<br/>

- Transdecoder is used to predict protein coding regions in the `Trinity.fasta` de novo assembly.  
`salloc --time=3:0:0 --ntasks=1 --mem=20G --account=def-coling`
`module load nixpkgs/16.09  gcc/5.4.0 transdecoder/5.5.0` <br/>
`TransDecoder.LongOrfs -t Trinity.fasta &> log.txt &` <br/>

- TransDecoder's output, `longest_orfs.pep`, is used throughout the trinotate pipeline <br/>
- Prior to the blast runs, I used `split_fasta.pl` from [here](https://github.com/gmarnellos/Trinotate_example_supplement/blob/master/split_fasta.pl) to parallelize the blast runs for both `longest_orfs.pep` and `Trinity.fasta`

- Blast arrays for bluegill: <br/>
`sbatch --array=1-50 blastx.sh`<br/>
`sbatch --array=1-32 blastp.sh`<br/>
`cat blastx.vol.*.outfmt6 > blastx.bluegill.sprot.outfmt6`<br/>
`cat blastp.vol.*.outfmt6 > blastp.bluegill.sprot.outfmt6`<br/>

- Blast arrays for lamprey: <br/>
`sbatch --array=1-61 blastx.sh`<br/>
`sbatch --array=1-60 blastp.sh`<br/>
`cat blastp.vol.*.outfmt6 > blastp.lamprey.sprot.outfmt6`<br/>
`cat blastx.vol.*.outfmt6 > blastx.lamprey.sprot.outfmt6`<br/>

- HMMER <br/>
`sbatch hmmer.sh`<br/>

- tmhmm <br/>
`sbatch tmhmm.sh` <br/>

- signalP <br/>
`sbatch signalp.sh` <br/>

- rnammer
`sbatch rnammer.sh`
- Of note, rnammer is incredibly memory intensive (502 Gb of the lamprey) and very difficult to get running.  Briefly: <br/>
- Download rnammer [here](https://services.healthtech.dtu.dk/service.php?RNAmmer-1.2)
- Download hmmer-2.3.2 [here](http://hmmer.org/download.html)
- Edit the rnammer scripts as suggested [here](https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required) and [here](https://blog.karinlag.no/2013/10/rnammer-install/)
- Perl might need modification.  On Compute Canada systems, you follow [this](https://docs.computecanada.ca/wiki/Perl)
`module load perl/5.22.4`<br/>
`module load gcc/5.4.0`<br/>
`cpan`<br/>
`cpan> install XML::Simple`<br/>

