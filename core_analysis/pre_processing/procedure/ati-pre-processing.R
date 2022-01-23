### R Script for Pre-Processing Around-the-Island Water Microbe Samples (16S) ###

#load libraries
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(decontam)

#Set working directory
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/pre_processing/procedure/")

#use qiime2R to upload data into a phyloseq object

#For final use the following code
physeq <- qza_to_phyloseq("../../../bioinformatics/output/ati-filtered-noeuk-table.qza", #feature table
                          "../../../bioinformatics/output/tree-building/ati-rooted-tree.qza", #tree
                          "../../../bioinformatics/output/ati-tax-without-spaces.qza", #taxonomy reference
                          "../../../bioinformatics/input/metadata_ati.txt") #mapping file


#Check the taxonomic classification is correct
rank_names(physeq)
tax_table(physeq)


##Check the Mitochondrial reads
mito <- subset_taxa(physeq, Family == "Mitochondria")
mito_taxa <- rownames(otu_table(mito))
#Read in your fasta file of representative sequences
#Filter based on your mito_taxa list 
#And export the new fasta file so that you can use Blastn to identify them
library(seqinr)
#Have to find these fasta files in the large output file (Not located on GitHub due to file size)
fasta.repseqs <- read.fasta(file = "../../../large_output/sequences.fasta", seqtyp = "DNA", as.string = TRUE, whole.header = TRUE)
mito.fasta <- fasta.repseqs[c(which(names(fasta.repseqs) %in% mito_taxa))]
write.fasta(mito.fasta, names(mito.fasta), "../../../large_output/mito_sequences.fasta")
#After this we need to make a list of the bad sequences by identifier

#Here check ASVs on BLASTn for mitochondrial vs. bacterial reads
#Keep mitochondrial reads in CSV or txt and re-upload as "bad taxa" to remove
bad.taxa <- read_table("../output/mitochondrial_removal/mitos_to_remove.txt")
#remove mitos them from the dataset
all.taxa <- taxa_names(physeq)
all.taxa <- all.taxa[!(all.taxa %in% bad.taxa)]
physeq.nm <- prune_taxa(all.taxa, physeq)

#make a sample data frame
sample.data <- as(sample_data(physeq.nm), "data.frame") #89 observations, 14 variables

#Decontamination cannot be done because control sample had 2 reads only! Clean samples! 

#Remove singletons
physeq.nm <- prune_taxa(taxa_sums(physeq.nm) > 1, physeq.nm)

#Check final numbers
tax_table(physeq.nm) 
sample_data(physeq.nm) #89 samples total, 14 variables
#Count Sequences
sum(sample_sums(physeq.nm)) #1,112,461

##Make two finalized phyloseq objects, one non-rarefied and one rarefied
physeq_nonrare <- physeq.nm
data.nonrare <- as(sample_data(physeq_nonrare), "data.frame")
saveRDS(physeq_nonrare, file = "../output/ati-physeq.RDS", compress = TRUE)

#Make a rarefied phyloseq object for alpha diversity analyses
physeq_rare <- rarefy_even_depth(physeq_nonrare, sample.size = 1000, rngseed = 711) #Set seed to be reproducible
sample_sums(physeq_rare) #Double check that the sub-sampling worked, this should report 1000 for each sample
data.rare <- as(sample_data(physeq_rare), "data.frame")
saveRDS(physeq_rare, file = "../output/ati-physeq-1000.RDS", compress = TRUE)



