# BUSCO was run from a docker image on biology-01
`sudo docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.1.4_cv1 busco -m transcriptome -c 8 -i Trinity.fasta -o bluegill_euks -l eukaryota_odb10 &> busco_bluegill_eukaryotes.txt &`
