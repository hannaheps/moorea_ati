cd ~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/bioinformatics/procedure
##Work out of the procedures folder & save to input vs. output

source activate qiime2-2020.11

#Raw sequence data are not kept in git repository because the files are too large
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../../../rawdata/ati_water_16S \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ../output/ati-demux-paired-end.qza

#Visualise the sequence quality summary to detect where to truncate sequences in next step
qiime demux summarize \
  --i-data ../output/ATI-demux-paired-end.qza \
  --o-visualization ../output/visualization/ati-demux-paired.qzv
  
#Dada2 denoising & merging 
#End with representative sequences & a feature table
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ../output/ati-demux-paired-end.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 17 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --o-representative-sequences ../output/ati-rep-seqs.qza \
  --o-table ../output/ati-table.qza

qiime feature-table filter-samples \
  --i-table ../output/ati-table.qza \
  --p-min-frequency 5000 \ #check what this can be - remove samples that have below this # of reads
  --o-filtered-table ../output/ati-filtered-table.qza

##Add metadata to these
##Metadata MUST be in tab separated format
qiime feature-table summarize \
  --i-table ../output/ati-filtered-table.qza  \
  --o-visualization ../output/visualization/ati-filtered-table.qzv  \
  --m-sample-metadata-file ../input/ati-metadata.txt

qiime feature-table tabulate-seqs \
  --i-data ../output/ati-rep-seqs.qza  \
  --o-visualization ../output/visualizations/ati-rep-seqs.qzv


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
  --input-path ../output/taxonomy-with-spaces/taxonomy-as-metadata/metadata.tsv \
  --output-path ../output/ati-tax-without-spaces.qza

#The following code spits out a table of features, or ASVs, present in sample set
#You can view by dragging the .qzv file into view.qiime2.org
#Or use the "qiime tools view" command

qiime metadata tabulate \
  --m-input-file ../output/ati-tax-without-spaces.qza \
  --o-visualization ../output/visualizations/ati-tax-without-spaces.qzv
  
### Filtering your table: removing mitochondria and chloroplast reads, PCR/sequencing errors, etc###
##Filter table

qiime taxa filter-table \
  --i-table ../output/ati-filtered-table.qza \
  --i-taxonomy ../output/ati-tax-without-spaces.qza \
  --p-exclude chloroplast,eukaryote,mitochondria \ 
  --o-filtered-table ../output/ati-filtered-noeuk-table.qza
   
qiime feature-table summarize \
  --i-table ../output/ati-filtered-noeuk-table.qza \
  --o-visualization ../output/visualizations/ati-filtered-noeuk-table.qzv
  
##Create a phylogenetic tree (With large datasets this takes ALOT of RAM - give it 50G)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ..output/ati-rep-seqs.qza \
  --o-alignment ..output/tree-building/ati-aligned-rep-seqs.qza \
  --o-masked-alignment ..output/tree-building/ati-masked-aligned-rep-seqs.qza \
  --o-tree ..output/tree-building/ati-unrooted-tree.qza \
  --o-rooted-tree ..output/tree-building/ati-rooted-tree.qza
 
  