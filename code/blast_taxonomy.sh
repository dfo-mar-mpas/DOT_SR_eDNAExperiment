#Blasting fasta seqs from dada2 
#Need to make sure I'm pointing the blast task at the right database, which also has taxdb in that path for getting species names
#blastn and blastx should be in the PATH so don't need to type their paths

export BLASTDB=/mnt/nvme1n1p1/Cody/blastdb


blastn -db nt_euk -query dna-sequences.fasta \
-max_target_seqs 5 \
-out 12Sblast_results.tsv \
-outfmt "6 qseqid sseqid pident length ssciname sblastname staxid" \
-num_threads 20
