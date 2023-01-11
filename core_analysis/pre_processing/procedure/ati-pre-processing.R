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

#Make sure you set working directory into the "pre-processing/procedure folder
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/pre_processing/procedure/")

#use qiime2R to upload data into a phyloseq object
#metadata comes from the output of the combine_meta.R script where sequencing metadata is combined with the full ati dataset


#For final use the following code
physeq <- qza_to_phyloseq("../../../bioinformatics/new_analysis/output/ati-2021-noeuk-table.qza", #feature table
                          "../../../bioinformatics/new_analysis/output/tree-building/ati-2021-rooted-tree.qza", #tree
                          "../../../bioinformatics/new_analysis/output/ati-2021-tax.qza", #taxonomy reference
                          "../../../metadata/new_metadata/output/metadata_ati_full_2021.txt") #mapping file

#Check the taxonomic classification is correct
rank_names(physeq) #7 ranks
tax_table(physeq) #37542 taxa

#Remove unassigned reads (unassigned at the kingdom level) & for good measure the euks
physeq <- subset_taxa(physeq,  Kingdom != "Unassigned")
physeq <- subset_taxa(physeq, Kingdom != "d__Eukaryota")

##All Family == Mitochondria are Rickettsiales, so no need to do extra Mitochondrial removals

#make a sample data frame
sample.data <- as(sample_data(physeq), "data.frame") #191 observations, 18 variables
View(sample.data)

##Decontamination
#Start by looking at the library size
sample.data$LibrarySize <- sample_sums(physeq)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))
#visualize library size
ggplot(data = sample.data, aes(x=Index, y=LibrarySize, color = Sample_type)) +
  geom_point()
##Library size of negative controls are way lower than samples (yay!)

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_type == "control"
contamdf.prev <- isContaminant(physeq, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
#Only 10 potential contaminants!

#Remove contaminants :) 
physeq.noncont <- prune_taxa(!contamdf.prev$contaminant, physeq)

#Remove control samples now because we already dealt with potential contaminants from them & are no longer needed
physeq.noncont <- subset_samples(physeq.noncont, is.neg != "TRUE")

#Remove singletons
physeq.final <- prune_taxa(taxa_sums(physeq.noncont) > 1, physeq.noncont)

#Check final numbers
tax_table(physeq.final) #37454 taxa
sample_data(physeq.final) #195 samples total, 50 variables


#Count Sequences
sum(sample_sums(physeq.final)) #6,715,723
#check sample sums per sample, sorted in ascending order
sort(sample_sums(physeq.final)) #lowest is sample 75 with 10,534, highest is 65,409

##Make two finalized phyloseq objects, one non-rarefied and one rarefied
#physeq.final is non rarefied
saveRDS(physeq.final, file = "../output/ati-2021-physeq-nonrare.RDS", compress = TRUE)

#Make a rarefied phyloseq object for alpha diversity analyses
physeq_rare <- rarefy_even_depth(physeq.final, sample.size = 10534, rngseed = 711) #Set seed to be reproducible
#Note: 8472 ASVs were removed due to subsampling as they were no longer present
sample_sums(physeq_rare) #Double check that the sub-sampling worked, this should report 10534 for each sample
saveRDS(physeq_rare, file = "../output/ati-2021-physeq-rare.RDS", compress = TRUE)

##Export the raw otu and taxonomy file for use in downstream analyses
otu <- otu_table(physeq.final)
tax <- tax_table(physeq.final)
otu.rare <- otu_table(physeq_rare)
tax.rare <- tax_table(physeq_rare)
write.csv(otu, "../output/ati-2021-asv-table.csv")
write.csv(tax, "../output/ati-2021-tax-table.csv")
write.csv(otu, "../output/ati-2021-rarefied-asv-table.csv")
write.csv(tax, "../output/ati-2021-rarefied-tax-table.csv")
