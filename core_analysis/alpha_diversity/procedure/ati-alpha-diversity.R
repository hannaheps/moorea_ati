###R Script for Alpha Diversity: Around-the-island Water Microbe Samples (16S)###

##This is a little awkward because our sample data already has microbiome diversity
#due to the iterative process of the collaborative data collation
#however, this should produce the same numbers as what are already in the sample data frame.

#Load all the libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(cowplot)

## Set working directory
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/alpha_diversity/procedure/")

#Read in the phyloseq object (output from pre-processing R script)
#Using the rarefied data for alpha diversity
physeq_rare <- readRDS("../../pre_processing/output/ati-2021-physeq-rare.RDS") 
data.rare <- as(sample_data(physeq_rare), "data.frame")

###ALPHA DIVERSITY###
#Pre-processing
source("estimate_richness_wPD.R") #reads in the function for calculating alpha diversity, including Faith's Phylogenetic Div. 
erich <- estimate_richness_wPD(physeq_rare, measures = c("Observed", "Shannon", "FaithPD"))

##Add evenness as a function of shannon
erich$evenness <- erich$Shannon/log(erich$Observed)

#Fill out the rest of the data from the phyloseq object
erich <- cbind(erich, data.rare)
rownames(erich) <- rownames(data.rare)
View(erich)

##check the two microbiome richness columns to check that they are in fact the same!
View(df <- cbind(erich$Observed, erich$Microbial_Species_Richness))
#Thank god


##I will not change the names of these because we already have them in the data frame.
##But if it is for any reason necessary to save the metadata with all the alpha diversity metrics
#use the following code:

#Output the RDS file and the metadata as a CSV with alpha diversity metrics 
saveRDS(erich, file = "../output/ati-erich.RDS", compress = TRUE)
#erich <- readRDS("../output/ati-erich.RDS")
write.csv(erich, "../output/ati-2021-metadata-with-alphadiv.csv" )
