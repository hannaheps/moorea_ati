
## Set working directory
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/parrotfish_microbes/procedure/")

#Load all the libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(RColorBrewer)

##read in the rarefied phyloseq object
physeq_rare <- readRDS("../../pre_processing/output/ati-physeq-1000.RDS") 
data.rare <- as(sample_data(physeq_rare), "data.frame")

#Subset to taxa that overlap with Leila's parrotfish data
physeq.overlap <- subset_taxa(physeq_rare, Family == c("Alcanivoracaceae", "Bacteriovoracaceae", "Bacteroidaceae",
                                                       "Clostridiaceae", ""))
