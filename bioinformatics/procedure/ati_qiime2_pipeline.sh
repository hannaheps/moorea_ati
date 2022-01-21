cd ~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/bioinformatics/procedure
##Work out of the procedures folder & save to input vs. output

source activate qiime2-2020.11

#Raw sequence data are not kept in git repository because the files are too large
#Sequence data coming from CGRB from OSU need the use of a manifest file bc the file names are not formatted normally
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../input/manifest_ati.txt \
  --output-path ../output/ati-demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2

#Try just the forward reads bc the reverse reads were VERY bad quality
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ../input/manifest_ati_single.txt \
  --output-path ../output/ati-demux-single-end.qza \
  --input-format SingleEndFastqManifestPhred33V2


#Visualise the sequence quality summary to detect where to truncate sequences in next step
qiime demux summarize \
  --i-data ../output/ati-demux-paired-end.qza \
  --o-visualization ../output/visualizations/ati-demux-paired-end.qzv
  
#just single end
qiime demux summarize \
  --i-data ../output/ati-demux-single-end.qza \
  --o-visualization ../output/visualizations/ati-demux-single-end.qzv
  
#Dada2 denoising & merging 
#End with representative sequences & a feature table
#Looks like the primers were removed by sequencing facility
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ../output/ati-demux-paired-end.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 100 \
  --o-representative-sequences ../output/ati-rep-seqs.qza \
  --o-denoising-stats ../output/ati-denoising-stats.qza \
  --o-table ../output/ati-table.qza
  
qiime metadata tabulate \
  --m-input-file ../output/ati-denoising-stats.qza \
  --o-visualization ../output/visualizations/ati-denoising-stats.qzv
  
#single end denoising
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ../output/ati-demux-single-end.qza \
  --p-trunc-len 240 \
  --o-table ../output/ati-single-table.qza \
  --o-representative-sequences ../output/ati-single-rep-seqs.qza \
  --o-denoising-stats ../output/ati-single-denoising-stats.qza

qiime metadata tabulate \
  --m-input-file ../output/ati-single-denoising-stats.qza \
  --o-visualization ../output/visualizations/ati-single-denoising-stats.qzv

#Moving forward with the single forward reads only due to low quality reverse reads
#Here we will remove samples that are less than 1000 reads 
qiime feature-table filter-samples \
  --i-table ../output/ati-single-table.qza \
  --p-min-frequency 1000 \
  --o-filtered-table ../output/ati-filtered-table.qza

##Add metadata to these
##Metadata MUST be in tab separated format
qiime feature-table summarize \
  --i-table ../output/ati-filtered-table.qza  \
  --o-visualization ../output/visualizations/ati-filtered-table.qzv  \
  --m-sample-metadata-file ../input/metadata_ati.txt

qiime feature-table tabulate-seqs \
  --i-data ../output/ati-single-rep-seqs.qza  \
  --o-visualization ../output/visualizations/ati-single-rep-seqs.qzv


#Creating/training a Naive Bayes classifier on the SILVA database
#First, import the silva database sequences (.fasta file - make sure sequences are in uppercase)
#and the taxonomy (.txt file) and change both to .qza artifacts
#Check for latest version of SILVA - Nov 2019 it is v132, but v138 should be out soon
#The SILVA files are not kept in the git repository because they are too big

qiime tools import \
 --type 'FeatureData[Sequence]' \
 --input-path ../../../../../../Biocomputing/databases/silva_138_99_16S.fna \
 --output-path ../../../../../../Biocomputing/databases/silva-138-bac.qza

#Extract the reference reads from SILVA for your dataset - this takes a long time!
#This is based on your primer set that you use. 
#The following code is based on 515/806R (Parada - Apprill, updated) primers, but you can substitute for any other primer set

qiime feature-classifier extract-reads \
 --i-sequences ../../../../../../Biocomputing/databases/silva-138-bac.qza \
 --p-f-primer GTGYCAGCMGCCGCGGTAA \
 --p-r-primer GGACTACNVGGGTWTCTAAT \
 --o-reads  ../../../../../../Biocomputing/databases/silva-v4-ref-seqs.qza

#Extract the reference taxonomy
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ../../../../../../Biocomputing/databases/taxonomy_all_levels.txt \
  --output-path ../../../../../../Biocomputing/databases/silva-v4-ref-taxonomy.qza
	
#Next, train a "Naive Bayes classifier" using the reference sequences we just extracted above
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ../../../../../../Biocomputing/databases/silva-v4-ref-seqs.qza \
  --i-reference-taxonomy ../../../../../../Biocomputing/databases/silva-v4-ref-taxonomy.qza \
  --o-classifier ../../../../../../Biocomputing/databases/classifier-silva-138-v4.qza

#Okay now we can pull back into our own git repository to work on our own data
#The above classifier can now be used for other analyses as well
#This next step can take a long time and needs lots of memory!
qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../../../Biocomputing/databases/classifier-silva-138-v4.qza \
  --i-reads ../output/ati-rep-seqs.qza \
  --o-classification ../output/ati-tax.qza

#For a trial with single end use this silva classifier trained on v132 previously
qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../../primer_trial/classifier_silva.qza \
  --i-reads ../output/ati-single-rep-seqs.qza \
  --o-classification ../output/ati-tax.qza


##A bug was found in the output from feature-classifier classify-sklearn where the taxon column has trailing white spaces
#e.g., "D_0__Bacteria;D-1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales. " <=trailing space at end of string
#Use the following code to remove the trailing white spaces in the taxonomy artefact:

qiime tools export \
  --input-path ../output/ati-tax.qza \
  --output-path ../output/taxonomy-with-spaces
  
qiime metadata tabulate \
  --m-input-file ../output/taxonomy-with-spaces/taxonomy.tsv  \
  --o-visualization ../output/taxonomy-with-spaces/taxonomy-as-metadata.qzv

qiime tools export \
 --input-path ../output/taxonomy-with-spaces/taxonomy-as-metadata.qzv \
 --output-path ../output/taxonomy-with-spaces/taxonomy-as-metadata

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path ../output/taxonomy-with-spaces/taxonomy-as-metadata/metadata_nohash.tsv \
  --output-path ../output/ati-tax-without-spaces.qza

#The following code spits out a table of features, or ASVs, present in sample set
#You can view by dragging the .qzv file into view.qiime2.org
#Or use the "qiime tools view" command

#There is a # sign in the qza file?
qiime metadata tabulate \
  --m-input-file ../output/ati-tax-without-spaces.qza \
  --o-visualization ../output/visualizations/ati-tax-without-spaces.qzv
  
### Filtering your table: removing mitochondria and chloroplast reads, PCR/sequencing errors, etc###
##Filter table

qiime taxa filter-table \
  --i-table ../output/ati-filtered-table.qza \
  --i-taxonomy ../output/ati-tax-without-spaces.qza \
  --p-exclude chloroplast,eukaryote \
  --o-filtered-table ../output/ati-filtered-noeuk-table.qza
   
qiime feature-table summarize \
  --i-table ../output/ati-filtered-noeuk-table.qza \
  --o-visualization ../output/visualizations/ati-filtered-noeuk-table.qzv
  
##Create a phylogenetic tree (With large datasets this takes ALOT of RAM - give it 50G)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ../output/ati-single-rep-seqs.qza \
  --o-alignment ../output/tree-building/ati-aligned-rep-seqs.qza \
  --o-masked-alignment ../output/tree-building/ati-masked-aligned-rep-seqs.qza \
  --o-tree ../output/tree-building/ati-unrooted-tree.qza \
  --o-rooted-tree ../output/tree-building/ati-rooted-tree.qza
 
  