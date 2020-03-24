# Trinity was run on Biology-01 with the following command
# For lamprey, the run crashed and was restarted, adding 2 to the log file

Trinity --max_memory 250G --seqType fq --samples_file samples.txt --KMER_SIZE 25 --SS_lib_type FR --CPU 24 --bflyCalculateCPU --output full_trinity_assembly &>log.txt &
