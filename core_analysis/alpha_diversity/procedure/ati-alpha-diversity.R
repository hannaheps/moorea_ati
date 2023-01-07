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
physeq_rare <- readRDS("../../pre_processing/output/ati-physeq-11146.RDS") 
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

#Output the RDS file and the metadata as a CSV with alpha diversity metrics 
saveRDS(erich, file = "../output/ati-erich.RDS", compress = TRUE)
#erich <- readRDS("../output/ati-erich.RDS")
write.csv(erich, "../output/ati-metadata-with-alphadiv.csv" )
#erich <- readRDS("../output/ati-erich.RDS")


###Can we do a super quick visual correlation of alpha diversity and nutrients??
p1 <- ggplot(erich, aes(x = water_silicate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Silicate") +
  ylab("Microbial Species Richness") +
  theme_bw()

p2 <- ggplot(erich, aes(x = water_silicate, y = Shannon)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Silicate") +
  ylab("Microbial Shannon Diversity") +
  theme_bw()

p3 <- ggplot(erich, aes(x = water_silicate, y = FaithPD)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Silicate") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()

p4 <- ggplot(erich, aes(x = water_silicate, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Silicate") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p1, p2, p3, p4) #from library "cowplot"
ggsave("../output/plots/microbialdiv_vs_silicate.pdf", plot = last_plot())


#Trial the same as above but with the %N

p5 <- ggplot(erich, aes(x = turb_percent_N, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("%N") +
  ylab("Microbial Species Richness") +
  theme_bw()

p6 <- ggplot(erich, aes(x = turb_percent_N, y = Shannon)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("%N") +
  ylab("Microbial Shannon Diversity") +
  theme_bw()

p7 <- ggplot(erich, aes(x = turb_percent_N, y = FaithPD)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("%N") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()

p8 <- ggplot(erich, aes(x = turb_percent_N, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("%N") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p5, p6, p7, p8) #from library "cowplot"

ggsave("../output/plots/microbialdiv_percN.pdf", plot = last_plot())

#Observed vs. water Nutrients
p9 <- ggplot(erich, aes(x = water_silicate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Silicate") +
  ylab("Microbial Species Richness") +
  theme_bw()

p10 <- ggplot(erich, aes(x = water_phosphate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Phosphate") +
  ylab("Microbial Species Richness") +
  theme_bw()

p11 <- ggplot(erich, aes(x = water_ammonia, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Ammonia") +
  ylab("Microbial Species Richness") +
  theme_bw()

p12 <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Nitrite plus Nitrate") +
  ylab("Microbial Species Richness") +
  theme_bw()

plot_grid(p9, p10, p11, p12) #from library "cowplot"
ggsave("../output/plots/observed_water_nutrients.pdf", plot = last_plot())

##Relationship with nitrite vs nitrate so check all microbial metrics
p13 <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Nitrite plus Nitrate") +
  ylab("Microbial Species Richness") +
  theme_bw()

p14 <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = Shannon)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Nitrite plus Nitrate") +
  ylab("Microbial Shannon Diversity") +
  theme_bw()

p15 <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = FaithPD)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Nitrite plus Nitrate") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()

p16 <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Nitrite plus Nitrate") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p13, p14, p15, p16) #from library "cowplot"

ggsave("../output/plots/microbialdiv_nitrite_nitrate.pdf", plot = last_plot())
