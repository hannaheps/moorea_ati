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
                          "../../../bioinformatics/output/ati-tax.qza", #taxonomy reference
                          "../../../bioinformatics/input/metadata_ati.txt") #mapping file


#Check the taxonomic classification is correct
rank_names(physeq)
tax_table(physeq)


##Check the Mitochondrial reads
mito <- subset_taxa(physeq, Family == "Mitochondria")
mito_taxa <- rownames(otu_table(mito))
write.csv(mito_taxa, "../output/mitochondrial_removal/mito_asvs.csv")
#Here check ASVs on BLAST for mitochondrial vs. bacterial reads
#Keep mitochondrial reads in CSV and re-upload as "bad taxa" to remove
bad.taxa <- read.csv("mitochondrial_removal/output/sequences_to_remove.csv", header = FALSE)
bad.taxa <- levels(bad.taxa$V1)
#remove them from the dataset
all.taxa <- taxa_names(physeq.filt)
all.taxa <- all.taxa[!(all.taxa %in% bad.taxa)]
physeq.nm <- prune_taxa(all.taxa, physeq.filt)

#make a sample data frame
sample.data <- as(sample_data(physeq.nm), "data.frame") #423 observations, 170 variables

#Decontamination cannot be done because control sample had 2 reads only 


#Check final numbers
tax_table(physeq.nm) 
sample_data(physeq.nm) #89 samples total, 14 variables
#Count Sequences
sum(sample_sums(physeq.nm)) #1,112,461

#Make another phyloseq object that is rarefied
#Rarefy to 1000 reads
physeq_rare <- rarefy_even_depth(physeq, sample.size = 1000, rngseed = 711) 

saveRDS(physeq.nm, "../output/ati-physeq.RDS")
saveRDS(physeq_rare, "../output/ati-physeq-1000.RDS")


