### 16S QIIME2 Pipeline for Around-the-Island (ATI) microbial community analyses###
### 2021 ONLY ####

#Make sure you are in the github repository "moorea_ati" and cd to bioinformatics/new_analysis/procedure

cd /moorea_ati/bioinformatics/new_analysis/procedure
##Starat in the 'procedure' folder
##Work from input and save to output

source activate qiime2-2020.11

##Important Note:
#Raw sequence data are not kept in git repository because the files are TOO large
#For raw sequence data, please contact either H. Epstein or D. Silva at OSU (this will be replaced with public repository on publication)
#Sequence data coming from the CQLS from OSU need the use of a manifest file bc the file names are not formatted normally

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../input/manifest_microbe_ati_2021.txt \
  --output-path ../output/ati-2021-demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2


#Visualise the sequence quality summary to detect where to truncate sequences in next step
qiime demux summarize \
  --i-data ../output/ati-2021-demux-paired-end.qza \
  --o-visualization ../output/visualizations/ati-2021-demux-paired-end.qzv

qiime tools view ../output/visualizations/ati-2021-demux-paired-end.qzv


##Spot checked the fasta files and it looks like the primers have already been removed by sequencing facility
##Might as well double check for all samples using the computer...
#Use plug-in cutadapt to make sure the primers are removed prior to proceeding to dada2
#This in general is better than simply trimming to the length of the primers as it searches specifically for
#primer sets that are provided. Here we used the following 515/806 primers (Parada & Apprill)
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ../output/ati-2021-demux-paired-end.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --o-trimmed-sequences ../output/ati-2021-demux-paired-end-trim.qza

#Always good practice to visualise and make sure all looks good from this step
qiime demux summarize \
  --i-data ../output/ati-2021-demux-paired-end-trim.qza \
  --o-visualization ../output/visualizations/ati-2021-demux-paired-end-trim.qzv

qiime tools view ../output/visualizations/ati-2021-demux-paired-end-trim.qzv

##Looks good, let's continue!

#Dada2 denoising & merging 
#You should end with representative sequences & a feature table
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ../output/ati-2021-demux-paired-end-trim.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 220 \
  --o-representative-sequences ../output/ati-2021-rep-seqs.qza \
  --o-denoising-stats ../output/ati-2021-denoising-stats.qza \
  --o-table ../output/ati-2021-table.qza
  
qiime metadata tabulate \
  --m-input-file ../output/ati-2021-denoising-stats.qza \
  --o-visualization ../output/visualizations/ati-2021-denoising-stats.qzv
  
#Check the status of the denoising
qiime tools view ../output/visualizations/ati-2021-denoising-stats.qzv
  
## filter out samples less than 1000 reads 
## this step can also be done in R, especially if negative controls were low read count
#All samples are above 1000 reads, except for one of the negative controls
#This means I will not use the following code block & hashed it out:

#qiime feature-table filter-samples \
#  --i-table ../output/ati-2021-table.qza \
#  --p-min-frequency 1000 \
#  --o-filtered-table ../output/ati-2021-filtered-table.qza

##Add metadata to these to visualize 
##Metadata MUST be in tab separated format
qiime feature-table summarize \
  --i-table ../output/ati-2021-table.qza  \
  --o-visualization ../output/visualizations/ati-2021-table.qzv  \
  --m-sample-metadata-file ../../../metadata/new_metadata/input/metadata_microbe_ati_2021.txt

qiime feature-table tabulate-seqs \
  --i-data ../output/ati-2021-rep-seqs.qza  \
  --o-visualization ../output/visualizations/ati-2021-rep-seqs.qzv

##Here is a good place to check one or two of the seqs via blast - are they bacteria??
qiime tools view ../output/visualizations/ati-2021-table.qzv
qiime tools view ../output/visualizations/ati-2021-rep-seqs.qzv
#yes! #1 top most abundant taxon in our samples is Prochlorococcus!


#Creating/training a Naive Bayes classifier on the SILVA database
#First, import the silva database sequences (.fasta file - make sure sequences are in uppercase)
#and the taxonomy (.txt file) and change both to .qza artifacts
#Check for latest version of SILVA - currently v138 
#The SILVA files are not kept in the git repository because they are TOO big

##SILVA v138 is in a new file format, so instead of creating a classifier locally
#we are using the one created by QIIME2 - this is not ideal, but currently best option
#(find it here: https://docs.qiime2.org/2022.11/data-resources/#taxonomy-classifiers-for-use-with-q2-feature-classifier)

#located locally here:

/Users/hannah/Documents/OSUDocs/Biocomputing/databases/SILVA_138/silva-138-99-seqs-515-806.qza
 
 ##Skip the following codes that are commented out below

#qiime tools import \
# --type 'FeatureData[Sequence]' \
# --input-path ../../../../../../Biocomputing/databases/silva_138_99_16S.fna \
# --output-path ../../../../../../Biocomputing/databases/silva-138-bac.qza

#Extract the reference reads from SILVA for your dataset - this takes a long time!
#This is based on your primer set that you use. 
#The following code is based on 515/806R (Parada - Apprill, updated) primers, but you can substitute for any other primer set

#qiime feature-classifier extract-reads \
# --i-sequences ../../../../../../Biocomputing/databases/silva-138-bac.qza \
# --p-f-primer GTGYCAGCMGCCGCGGTAA \
# --p-r-primer GGACTACNVGGGTWTCTAAT \
# --o-reads  ../../../../../../Biocomputing/databases/silva-v4-ref-seqs.qza

#Extract the reference taxonomy
#qiime tools import \
#  --type 'FeatureData[Taxonomy]' \
#  --input-format HeaderlessTSVTaxonomyFormat \
#  --input-path ../../../../../../Biocomputing/databases/taxonomy_all_levels.txt \
#  --output-path ../../../../../../Biocomputing/databases/silva-v4-ref-taxonomy.qza
	
#Next, train a "Naive Bayes classifier" using the reference sequences we just extracted above
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads ../../../../../../Biocomputing/databases/silva-v4-ref-seqs.qza \
#  --i-reference-taxonomy ../../../../../../Biocomputing/databases/silva-v4-ref-taxonomy.qza \
#  --o-classifier ../../../../../../Biocomputing/databases/classifier-silva-138-v4.qza

#Okay now we can pull back into our own git repository to work on our own data
#The above classifier can now be used for other analyses as well
#Note: This next step can take a long time and needs lots of memory!
qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../../../../Biocomputing/databases/SILVA_138/classifier-silva-138.qza \
  --i-reads ../output/ati-2021-rep-seqs.qza \
  --o-classification ../output/ati-2021-tax.qza

##Trial skipping the next steps first, the bug might be gone! 
##Looks like the bug is gone so no need for the following steps (hashed out)

##A bug was found in the output from feature-classifier classify-sklearn where the taxon column has trailing white spaces
#e.g., "D_0__Bacteria;D-1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales. " <=trailing space at end of string
#Use the following code to remove the trailing white spaces in the taxonomy artefact:

#qiime tools export \
#  --input-path ../output/ati-tax.qza \
#  --output-path ../output/taxonomy-with-spaces
  
#qiime metadata tabulate \
#  --m-input-file ../output/taxonomy-with-spaces/taxonomy.tsv  \
#  --o-visualization ../output/taxonomy-with-spaces/taxonomy-as-metadata.qzv

#qiime tools export \
# --input-path ../output/taxonomy-with-spaces/taxonomy-as-metadata.qzv \
# --output-path ../output/taxonomy-with-spaces/taxonomy-as-metadata

#qiime tools import \
#  --type 'FeatureData[Taxonomy]' \
#  --input-path ../output/taxonomy-with-spaces/taxonomy-as-metadata/metadata_nohash.tsv \
#  --output-path ../output/ati-2021-tax-without-spaces.qza

#The following code spits out a table of features, or ASVs, present in sample set
#You can view by dragging the .qzv file into view.qiime2.org
#Or use the "qiime tools view" command

#There is a # sign in the qza file?? Check this file pls
qiime metadata tabulate \
  --m-input-file ../output/ati-2021-tax.qza \
  --o-visualization ../output/visualizations/ati-2021-tax.qzv

qiime tools view ../output/visualizations/ati-2021-tax.qzv

  
### Filtering your table: removing mitochondria and chloroplast reads, PCR/sequencing errors, etc###
##Filter table

qiime taxa filter-table \
  --i-table ../output/ati-2021-table.qza \
  --i-taxonomy ../output/ati-2021-tax.qza \
  --p-exclude Eukaryota \
  --o-filtered-table ../output/ati-2021-noeuk-table.qza
   
qiime feature-table summarize \
  --i-table ../output/ati-2021-noeuk-table.qza \
  --o-visualization ../output/visualizations/ati-2021-noeuk-table.qzv

qiime tools view ../output/visualizations/ati-2021-noeuk-table.qzv

  
##Create a phylogenetic tree (With large datasets this takes ALOT of RAM - give it 50G)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ../output/ati-2021-rep-seqs.qza \
  --o-alignment ../output/tree-building/ati-2021-aligned-rep-seqs.qza \
  --o-masked-alignment ../output/tree-building/ati-2021-masked-aligned-rep-seqs.qza \
  --o-tree ../output/tree-building/ati-2021-unrooted-tree.qza \
  --o-rooted-tree ../output/tree-building/ati-2021-rooted-tree.qza
  
##Check for alpha rarefaction
##NOTE:make sure your metadata file does NOT have a variable column named "depth"
##Otherwise it will return an error: "Plugin error from diversity: cannot insert depth, already exists"
##as the code will add a depth variable for rarefaction
##Thus, I changed the metadata column to "water_depth" so that it does not interfere

qiime diversity alpha-rarefaction \
  --i-table ../output/ati-2021-noeuk-table.qza \
  --i-phylogeny ../output/tree-building/ati-2021-rooted-tree.qza \
  --p-max-depth 12000 \
  --m-metadata-file ../../../metadata/new_metadata/input/metadata_microbe_ati_2021.txt \
  --o-visualization ../output/ati-2021-alpha-rarefaction.qzv

qiime tools view ../output/ati-2021-alpha-rarefaction.qzv

qiime taxa barplot \
  --i-table ../output/ati-2021-noeuk-table.qza \
  --i-taxonomy ../output/ati-2021-tax.qza \
  --m-metadata-file ../../../metadata/new_metadata/input/metadata_microbe_ati_2021.txt \
  --o-visualization ../output/ati-2021-barplot.qzv

qiime tools view ../output/ati-2021-barplot.qzv


##For downstream analyses you need the following four files:
#1. ati-2021-filtered-noeuk-table.qza
#2. ati-2021-rooted-tree.qza
#3. ati-2021-tax.qza
#4. metadata_microbe_ati_2021.txt 

 
  