# BUSCO was run from a docker image on biology-01
- for the original Trinity assembly:<br/>
`sudo docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.1.4_cv1 busco -m transcriptome -c 8 -i Trinity.fasta -o bluegill_euks -l eukaryota_odb10 &> busco_bluegill_eukaryotes.txt &`<br/>
- for the Corset/Lace SuperTranscriptome:<br/>
`sudo docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.1.4_cv1 busco -m transcriptome -c 8 -i SuperDuper.fasta -o bluegill_lace -l eukaryota_odb10 &> busco_bluegill_eukaryotes_superTranscriptome.txt &`
