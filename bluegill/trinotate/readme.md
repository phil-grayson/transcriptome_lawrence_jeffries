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

- Blast arrays: <br/>
`sbatch --array=1-50 blastx.sh`<br/>
`sbatch --array=1-32 blastp.sh`<br/>
`cat blastx.vol.*.outfmt6 > blastx.bluegill.sprot.outfmt6`<br/>
`cat blastp.vol.*.outfmt6 > blastp.bluegill.sprot.outfmt6`<br/>

- HMMER <br/>
`sbatch hmmer.sh`<br/>

- tmhmm <br/>
`sbatch tmhmm.sh` <br/>

- signalP <br/>
`sbatch signalp.sh` <br/>

- rnammer <br/>
`sbatch rnammer.sh`
- Of note, rnammer is incredibly memory intensive (502 Gb on the lamprey) and very difficult to get running.  Briefly: <br/>
  - Download rnammer [here](https://services.healthtech.dtu.dk/service.php?RNAmmer-1.2)
  - Download hmmer-2.3.2 [here](http://hmmer.org/download.html)
  - Edit the rnammer scripts as suggested [here](https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required) and [here](https://blog.karinlag.no/2013/10/rnammer-install/)
  - Perl might need modification.  On Compute Canada systems, you follow [this](https://docs.computecanada.ca/wiki/Perl) <br/>
  `module load perl/5.22.4`<br/>
  `module load gcc/5.4.0`<br/>
  `cpan`<br/>
  `cpan> install XML::Simple`<br/>

- With everything run, it's time to collect them into the sqlite database for Trinotate<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep longest_orfs.pep`<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.bluegill.sprot.outfmt6`<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.bluegill.sprot.outfmt6`<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out`<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out`<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_signalp signalp.out`<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite LOAD_rnammer Trinity.fasta.rnammer.gff`<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/Trinotate Trinotate.sqlite report > trinotate_annotation_report_bluegill_RF.xls`<br/>

- GO annotations were predicted as follows:<br/>
`~/programs/Trinotate-Trinotate-v3.2.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate_annotation_report_bluegill_RF.xls -G --include_ancestral_terms > go_annotations_bluegill_RF.txt`
