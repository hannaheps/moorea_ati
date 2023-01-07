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
physeq <- qza_to_phyloseq("../../../bioinformatics/new_analysis/output/ati-filtered-noeuk-table.qza", #feature table
                          "../../../bioinformatics/new_analysis/output/ati-rooted-tree.qza", #tree
                          "../../../bioinformatics/new_analysis/output/ati-tax.qza", #taxonomy reference
                          "../../../metadata/new_metadata/output/metadata_with_turb.txt") #mapping file

sample.data <- as(sample_data(physeq), "data.frame")
View(sample.data)
#Check the taxonomic classification is correct
rank_names(physeq)
tax_table(physeq)

physeq <- subset_taxa(physeq,  Kingdom != "Unassigned")

##ON 6 January 2023 these data use green genes so the mitochondria issue is no longer valid
##Check the Mitochondrial reads
#mito <- subset_taxa(physeq, Family == "Mitochondria")
#mito_taxa <- rownames(otu_table(mito))
#Read in your fasta file of representative sequences
#Filter based on your mito_taxa list 
#And export the new fasta file so that you can use Blastn to identify them
#library(seqinr)
#Have to find these fasta files in the large output file (Not located on GitHub due to file size)
#fasta.repseqs <- read.fasta(file = "../../../large_output/sequences.fasta", seqtyp = "DNA", as.string = TRUE, whole.header = TRUE)
#mito.fasta <- fasta.repseqs[c(which(names(fasta.repseqs) %in% mito_taxa))]
#write.fasta(mito.fasta, names(mito.fasta), "../../../large_output/mito_sequences.fasta")
#After this we need to make a list of the bad sequences by identifier

#Here check ASVs on BLASTn for mitochondrial vs. bacterial reads
#Keep mitochondrial reads in CSV or txt and re-upload as "bad taxa" to remove
#bad.taxa <- read_table("../output/mitochondrial_removal/mitos_to_remove.txt")
#remove mitos them from the dataset
#all.taxa <- taxa_names(physeq)
#all.taxa <- all.taxa[!(all.taxa %in% bad.taxa)]
#physeq.nm <- prune_taxa(all.taxa, physeq)

#To maintain consistency in the code, rename physeq to physeq.nm
physeq.nm <- physeq
#make a sample data frame
sample.data <- as(sample_data(physeq.nm), "data.frame") #191 observations, 18 variables
View(sample.data)

##Decontamination
##Currently as of Jan 6 2023 we cannot run decontamination because there are no controls in the dataset
#potential causes of this: they were accidentally filtered out during the qiime2 pipeline - ask Denise to share pipeline to check this & fix

sample.data$LibrarySize <- sample_sums(physeq.nm)
sample.data <- sample.data[order(sample.data$LibrarySize),]
sample.data$Index <- seq(nrow(sample.data))
#visualize
ggplot(data = sample.data, aes(x=Index, y=LibrarySize, color = sample_type)) +
  geom_point()

#Next check for contaminants using prevalence and threshold 0.5 (more conservative)
sample_data(physeq.nm)$is.neg <- sample_data(physeq.nm)$sample_type == "blank"
contamdf.prev <- isContaminant(physeq.nm, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
# X number of contaminants (can't get a # until we include contaminants in phyloseq object)

#Remove contaminants
physeq.noncont <- prune_taxa(!contamdf.prev$contaminant, physeq.nm)

#Remove control samples now because we already dealt with potential contaminants from them & are no longer needed
physeq.noncont <- subset_samples(physeq.noncont, is.neg != "TRUE")

#Remove singletons
physeq.nm <- prune_taxa(taxa_sums(physeq.nm) > 1, physeq.nm)

#Check final numbers
tax_table(physeq.nm) 
sample_data(physeq.nm) #191 samples total, 18 variables


#Count Sequences
sum(sample_sums(physeq.nm)) #6,700,972
sample_sums(physeq.nm) #check the sample sums for very sample, lowest = 11146 (this will change with decontamination)

##Make two finalized phyloseq objects, one non-rarefied and one rarefied
#physeq_nonrare <- physeq.nm
#data.nonrare <- as(sample_data(physeq_nonrare), "data.frame")
saveRDS(physeq.nm, file = "../output/ati-physeq.RDS", compress = TRUE)

#Make a rarefied phyloseq object for alpha diversity analyses
physeq_rare <- rarefy_even_depth(physeq.nm, sample.size = 11146, rngseed = 711) #Set seed to be reproducible
sample_sums(physeq_rare) #Double check that the sub-sampling worked, this should report 1000 for each sample
data.rare <- as(sample_data(physeq_rare), "data.frame")
saveRDS(physeq_rare, file = "../output/ati-physeq-11146.RDS", compress = TRUE)

##Export the raw otu and taxonomy file
otu <- otu_table(physeq.nm)
tax <- tax_table(physeq.nm)
write.csv(otu, "../output/ati-moorea-asv-table.csv")
write.csv(tax, "../output/ati-moorea-tax-table.csv")

