#!/usr/bin/bash 

#Processing eDNA metabarcode data using QIIME2 and cutadapt, followed by DADA2 taxonomy assignments using custom fish reference databases
#The qiime tutorials are useful and found at https://docs.qiime2.org/2022.2/tutorials/overview/#useful-points-for-beginners 
#First activate QIIME if it hasn't been, can also reactivate qiime if you close the window 
conda activate qiime2-amplicon-2023.9 && source tab-qiime #activate tab completion

#check currently active conda environment
conda info

# View plugin, used to view any .qzv file in html and export tsv and fasta files
qiime tools view /path_to_file/filename.qzv

#for this project we have 12S fish sequences and COI general biodiversity sequences
#First import our data using a 'manifest' file of all fastq file names

#12S
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-12Smanifest \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path 12S-combined-demux.qza

#COI 
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-DOTCOImanifest \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path COI-combined-demux.qza

  
#12S
qiime demux summarize \
  --i-data 12S-DOT-combined-demux.qza \
  --o-visualization 12S-DOT-demux-subsample.qzv ##save tsv file of per-sample-fastq-counts.tsv for optional step below ##
  
#COI
qiime demux summarize \
  --i-data COI-combined-demux.qza \
  --o-visualization COI-DOT-demux.qzv

  
 ## OPTIONAL: filter out samples with less than 100 reads (can set this to any number) ##
qiime demux filter-samples \
  --i-demux 16S-combined-demux.qza \
  --m-metadata-file /path_to_output_folder/per-sample-fastq-counts.tsv \
  --p-where 'CAST([forward sequence count] AS INT) > 100' \
  --o-filtered-demux /path_to_output_folder/filename_greater100reads.qza
 
  ###########Use fastp to trim
 
  #McInnes primers
'
for i in *_R1.fastq.gz;
do
  SAMPLE=$(echo ${i} | sed "s/_R1.fastq.gz//") 
  echo ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
  cutadapt -g AGCGYAATCACTTGTCTYTTAA \
  -G CRBGGTCGCCCCAACCRAA \
    --discard-untrimmed --nextseq-trim 20 -m 20 -n 2 \
    --match-read-wildcards \
    -o ${SAMPLE}_R1trimmed.fastq \
    -p ${SAMPLE}_R2trimmed.fastq \
    ${SAMPLE}_R1.fastq.gz  ${SAMPLE}_R2.fastq.gz ;
done
'


#can also trim primers/adapters in Qiime with cutadapt

#MiFishU-F - these primers target ~180bp of the 12S rDNA region
#5 GTCGGTAAAACTCGTGCCAGC 3
#MiFishU-R
#3 GTTTGACCCTAATCTATGGGGTGATAC 5 > need to reverse this so its read 5 prime to 3 prime in cutadapt

#12S MiFish primers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences 12S-combined-demux.qza \
--p-cores 40 \
--p-front-f GTCGGTAAAACTCGTGCCAGC \
--p-front-r NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
--p-error-rate 0.15 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 30 \
--o-trimmed-sequences 12s-demux-trimmed-2023-test2.qza \
--output-dir trimmed \
--verbose
#visualize the trimming results
qiime demux summarize --i-data 12s-demux-trimmed-2023-test2.qza \
--o-visualization 12s-trimmed-visual

#mICOIintF GGWACWGGWTGAACWGTWTAYCCYCC
#jgHCO2198 TAIACYTCIGGRTCICCRAARAAYCA

#COI Leray primers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences COI-combined-demux.qza \
--p-cores 40 \
--p-front-f GGWACWGGWTGAACWGTWTAYCCYCC \
--p-front-r TAIACYTCIGGRTCICCRAARAAYCA \
--p-error-rate 0.15 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 50 \
--o-trimmed-sequences coi-demux-dot-trimmed.qza \
--output-dir trimmed \
--verbose

'example output from COI
=== Summary ===

Total read pairs processed:            235,239
  Read 1 with adapter:                 234,765 (99.8%)
  Read 2 with adapter:                 232,812 (99.0%)

== Read fate breakdown ==
Pairs that were too short:              44,105 (18.7%)
Pairs discarded as untrimmed:            2,599 (1.1%)
Pairs written (passing filters):       188,535 (80.1%)
'
#visualize the trimming results
qiime demux summarize --i-data coi-demux-dot-trimmed.qza \
--o-visualization coi-trimmed-visual


#denoise using dada2 which infers ASVs 
#Note: using --p-n-threads = 0 will use all threads available 
### for 16S use --p-trunc-len-f 125 and --p-trunc-len-r 125; 12S use 116 and 108 ###
# can add --p-min-overlap 12 or some other number if need be
#COI
qiime dada2 denoise-paired \
--i-demultiplexed-seqs coi-demux-dot-trimmed.qza  \
--p-trunc-len-f  260 \
--p-trunc-len-r  260 \
--p-n-threads 0 \
--p-n-reads-learn 3000000 \
--p-pooling-method independent \
--output-dir dada2out \
--verbose


#Try COI with some 5' trimming too - this worked better in this case!
qiime dada2 denoise-paired \
--i-demultiplexed-seqs coi-demux-dot-trimmed.qza  \
--p-trunc-len-f  240 \
--p-trunc-len-r  240 \
--p-trim-left-f 5 \
--p-trim-left-r 5 \
--p-n-threads 0 \
--p-n-reads-learn 3000000 \
--p-pooling-method independent \
--output-dir dada2out_test2 \
--verbose


#12S - trying some different r-len truncs
qiime dada2 denoise-paired \
--i-demultiplexed-seqs 12s-DOT-demux-trimmed.qza \
--p-trim-left-r 3 \
--p-trunc-len-f  140 \
--p-trunc-len-r  134 \
--p-n-threads 0 \
--p-min-overlap 8 \
--p-pooling-method independent \
--output-dir DOT12S-denoised-Test2 \
--verbose


#Generate summaries of denoising stats and feature table
#12S
qiime feature-table summarize \
  --i-table dada2out/table.qza \
  --o-visualization dada2out/table.qzv \
  --m-sample-metadata-file 2023-DOT-metadata.tsv &&
qiime feature-table tabulate-seqs \
  --i-data DOT12S-denoised-Test2/representative_sequences.qza \
  --o-visualization DOT12S-denoised-Test2/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file DOT12S-denoised-Test2/denoising_stats.qza \
  --o-visualization DOT12S-denoised-Test2/denoising-stats.qzv
  
#COI
qiime feature-table summarize \
  --i-table dada2out/table.qza \
  --o-visualization dada2out/table.qzv \
  --m-sample-metadata-file 2023-DOT-metadata.tsv &&
qiime feature-table tabulate-seqs \
  --i-data dada2out/representative_sequences.qza \
  --o-visualization dada2out/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file dada2out/denoising_stats.qza \
  --o-visualization dada2out/denoising-stats.qzv
  
  
 qiime tools view /path_to_output_folder/filename_rep_seqs.qzv  ## export the ASV fasta file from the view for input into FuzzyID2 and BLAST



 ### export results to biom formatted file
qiime tools export \
--input-path dada2out/table.qza \
--output-path dada2out/DOTCOI_filtered_table_biom ##specifying a folder output here, this tool will automatically export a file called 'feature-table.biom' to this folder

### convert biom to tsv
biom convert -i dada2out_test2/DOTCOI_filtered_table_biom/feature-table.biom \
-o dada2out_test2/DOTCOI_filtered_table_biom/DOTCOI_feature_table_export.tsv \
--to-tsv

### OPTIONAL filtering after exporting to tsv
## Remove rare ASV's by calculating if an ASV has a read number that is less than 0.1% of the total read number of that ASV across all samples. 
## This is summing across columns in the exported feature table, calculating 0.1% of that sum, and removing all instances where read numbers were less than that number.
 
 #Generate a phylogenetic tree from our data
 cd dada2out/
 qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative_sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
  #now use the rooted tree to generate some biodiversity stats
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1500 \
  --p-n-jobs-or-threads auto \
  --m-metadata-file  ../2023-DOT-metadata.tsv \
  --output-dir COI-core-metrics-results
 
####################
######TAXONOMY######
#####################
## there are lots of ways to do taxonomy, including just blasting, or building a reference database using rescript (below), or using FuzzyID2 with a custom library
  #using rescript to train our classifier
  qiime rescript filter-taxa \
  --i-taxonomy fish-16S-ref-tax.qza \
  --m-ids-to-keep-file fish-16S-ref-seqs-keep.qza \
  --o-filtered-taxonomy fish-16S-ref-taxa-keep.qza
  
  qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-16S-ref-taxa-keep.qza \
 --o-taxonomy-stats fish-16S-ref-tax-keep-eval.qzv
 
 qiime metadata tabulate \
 --m-input-file fish-16S-ref-taxa-keep.qza \
 --o-visualization fish-16S-ref-tax-keep.qzv &&
 qiime rescript evaluate-seqs \
 --i-sequences fish-16S-ref-seqs-keep.qza \
 --p-kmer-lengths 32 16 8 \
 --o-visualization fish-16S-ref-seqs-keep-eval.qzv
 
 #Build and evaluate classifier
 #here, the --o-classifier output is of type TaxonomicClassifier and the -o-observed-taxonomy is FeatureData[Taxonomy] (same as --i-taxonomy)
 qiime rescript evaluate-fit-classifier \
 --i-sequences fish-16S-ref-seqs-keep.qza \
 --i-taxonomy fish-16S-ref-taxa-keep.qza \
 --p-n-jobs -1 \
 --o-classifier ncbi-16S-fish-refseqs-classifier.qza \
 --o-evaluation ncbi-16S-fish-refseqs-classifier-evaluation.qzv \
 --o-observed-taxonomy ncbi-16S-fish-refseqs-predicted-taxonomy.qza \
 --output-dir 16S-Classifier \
 --verbose 
 
 qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-16S-ref-taxa-keep.qza ncbi-16S-fish-refseqs-predicted-taxonomy.qza \
 --p-labels ref-taxonomy predicted-taxonomy \
 --o-taxonomy-stats 16S-ref-taxonomy-evaluation.qzv \
 --verbose
 
 
  qiime rescript evaluate-fit-classifier \
 --i-sequences fish-12S-ref-seqs-keep.qza \
 --i-taxonomy fish-12S-ref-taxa-keep.qza \
 --p-n-jobs -1 \
 --o-classifier ncbi-12S-fish-refseqs-classifier.qza \
 --o-evaluation ncbi-12S-fish-refseqs-classifier-evaluation.qzv \
 --o-observed-taxonomy ncbi-12S-fish-refseqs-predicted-taxonomy.qza \
 --output-dir 12S-Classifier \
 --verbose 
 
 qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-12S-ref-taxa-keep.qza ncbi-12S-fish-refseqs-predicted-taxonomy.qza \
 --p-labels ref-taxonomy predicted-taxonomy \
 --o-taxonomy-stats 12S-ref-taxonomy-evaluation.qzv \
 --verbose
 
#Now back to qiime to do our taxonomy
  qiime feature-classifier classify-sklearn \
  --i-classifier ../../../ReferenceData/Rescript_ReferenceDB/ncbi-16S-fish-refseqs-classifier.qza \
  --i-reads representative_sequences.qza \
  --o-classification 16S-taxonomy.qza

qiime metadata tabulate \
  --m-input-file 16S-taxonomy.qza \
  --o-visualization 16S-taxonomy.qzv
  
  qiime feature-classifier classify-sklearn \
  --i-classifier ../../../ReferenceData/Rescript_ReferenceDB/ncbi-12S-fish-refseqs-classifier.qza \
  --i-reads representative_sequences.qza \
  --o-classification 12S-taxonomy.qza
  
  qiime metadata tabulate \
  --m-input-file 12S-taxonomy.qza \
  --o-visualization 12S-taxonomy.qzv
  

