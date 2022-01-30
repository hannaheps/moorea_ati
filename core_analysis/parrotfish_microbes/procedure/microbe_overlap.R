
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
                                                       "Clostridiaceae", "Colwelliaceae", "Corynebacteriaceae", "Cryomorphaceae",
                                                       "Dermacoccaceae", "Desulfobulbaceae", "Dietziaceae", "Flammeovirgaceae", 
                                                       "Flavobacteriaceae", "Hahellaceae", "Hyphomicrobiaceae", "Lachnospiraceae",
                                                       "Moraxellaceae", "Paenibacillaceae", "Pirellulaceae", "Piscirickettsiaceae",
                                                       "Prevotellaceae", "Pseudoalteromonadaceae", "Pseudomonadaceae", "Puniceicoccaceae",
                                                       "Rhodobacteraceae", "Rhodospirillaceae", "Simkaniaceae", "Sphingomonadaceae",
                                                       "Spirochaetaceae", "Streptococcaceae", "Vibrionaceae", "Xanthomonadaceae", "Xenococcaceae"))



physeq.overlap <- subset_taxa(physeq_rare, Family == "Alcanivoracaceae" | Family == "Bacteriovoracaceae" | Family == "Bacteroidaceae" | 
                                Family == "Clostridiaceae" | Family == "Colwelliaceae" | Family =="Corynebacteriaceae" | Family == "Cryomorphaceae" |
                                Family == "Dermacoccaceae" | Family == "Desulfobulbaceae" | Family == "Dietziaceae" | Family == "Flammeovirgaceae" |
                                Family ==  "Flavobacteriaceae" | Family == "Hahellaceae" | Family == "Hyphomicrobiaceae" | Family == "Lachnospiraceae" |
                                Family == "Moraxellaceae" | Family == "Paenibacillaceae" | Family == "Pirellulaceae" | Family == "Piscirickettsiaceae" | 
                                Family == "Prevotellaceae" | Family == "Pseudoalteromonadaceae" | Family == "Pseudomonadaceae" | Family == "Puniceicoccaceae" |
                                Family == "Rhodobacteraceae" | Family == "Rhodospirillaceae" | Family == "Simkaniaceae" | Family == "Sphingomonadaceae" |
                                Family == "Spirochaetaceae" | Family == "Streptococcaceae" | Family == "Vibrionaceae" | Family == "Xanthomonadaceae" |
                                Family == "Xenococcaceae")


physeq.overlap.ra <- transform_sample_counts(physeq.overlap, function(x) x/ sum(x))
physeq.overlap.ra.melt <- psmelt(physeq.overlap.ra)
View(physeq.overlap.ra.melt)

keeps <- c("Site_number","Abundance", "Family")
family.abund <- physeq.overlap.ra.melt[ , (names(physeq.overlap.ra.melt) %in% keeps)]
  
rownames(family.abund) <- NULL
spread(family.abund, Family, Abundance)

# Julianna
trial <- family.abund %>% 
  group_by(Family, Site_number) %>% 
  summarize(abundance = sum(Abundance)) %>% view()

write.csv(family.abund, "../output/family_abund.csv")
library(sjmisc)
trial <- family.abund %>% rotate_df(family.abund)
View(trial)


