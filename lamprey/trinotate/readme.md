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
- Prior to the blast runs, I used `split_fasta.pl` from [here](https://github.com/gmarnellos/Trinotate_example_supplement/blob/master/split_fasta.pl) to parallelize the blast runs for both longest_orfs.pep and `Trinity.fasta`


