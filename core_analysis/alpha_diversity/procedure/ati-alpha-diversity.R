###R Script for Alpha Diversity: Around-the-island Water Microbe Samples (16S)###

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
physeq_rare <- readRDS("../../pre_processing/output/ati-physeq-1000.RDS") 
data.rare <- as(sample_data(physeq_rare), "data.frame")

###ALPHA DIVERSITY###
#Pre-processing
source("estimate_richness_wPD.R") #reads in the function for calculating alpha diversity, including Faith's Phylogenetic Div. 
erich <- estimate_richness_wPD(physeq_rare, measures = c("Observed", "Shannon", "FaithPD"))

rownames(erich) <- rownames(data.rare)
erich$sample.id <- rownames(erich)
erich$location <- data.rare$location
erich$site.number <- data.rare$site.number
erich$lat <- data.rare$lat
erich$long <- data.rare$long
erich$sample_type <- data.rare$sample_type
erich$collection_time <- data.rare$collection_time
erich$collection_date <- data.rare$collection_date
erich$Phosphate <- data.rare$Phosphate
erich$Silicate <- data.rare$Silicaate
erich$Nitrite_plus_Nitrate <- data.rare$Nitrite_plus_Nitrate
erich$Ammonia <- data.rare$Ammonia

#Output the RDS file and the metadata as a CSV with alpha diversity metrics 
saveRDS(erich, file = "../output/ati-erich.RDS", compress = TRUE)
write.csv(erich, "../output/ati-metadata-with-alphadiv.csv" )
#erich <- readRDS("../output/ati-erich.RDS")


###Can we do a super quick visual correlation of alpha diversity and nutrients??
p1 <- ggplot(erich, aes(x = Nitrite_plus_Nitrate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Nitrite + Nitrate") +
  ylab("Observed Microbial Species Richness") +
  theme_bw()

p2 <- ggplot(erich, aes(x = Phosphate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Phosphate") +
  ylab("Observed Microbial Species Richness") +
  theme_bw()

p3 <- ggplot(erich, aes(x = Ammonia, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Ammonia") +
  ylab("Observed Microbial Species Richness") +
  theme_bw()

##How about adding an N:P ratio so we can plot this too?
erich$N.P <- (erich$Ammonia + erich$Nitrite_plus_Nitrate)/erich$Phosphate

p4 <- ggplot(erich, aes(x = N.P, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("N:P Ratio") +
  ylab("Observed Microbial Species Richness") +
  theme_bw()

plot_grid(p1, p2, p3, p4) #from library "cowplot"
ggsave("../output/plots/alpha_vs_nutrients.pdf", plot = last_plot())



    
    
