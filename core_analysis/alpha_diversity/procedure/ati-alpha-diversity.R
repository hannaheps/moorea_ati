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

##Add evenness as a function of shannon
erich$evenness <- erich$Shannon/log(erich$Observed)

#Fill out the rest of the data from the phyloseq object
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
erich <- readRDS("../output/ati-erich.RDS")
write.csv(erich, "../output/ati-metadata-with-alphadiv.csv" )
#erich <- readRDS("../output/ati-erich.RDS")


###Can we do a super quick visual correlation of alpha diversity and nutrients??
p1 <- ggplot(erich, aes(x = Nitrite_plus_Nitrate, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Nitrite + Nitrate") +
  ylab("Microbial Evenness") +
  theme_bw()

p2 <- ggplot(erich, aes(x = Phosphate, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Phosphate") +
  ylab("Microbial Evenness") +
  theme_bw()

p3 <- ggplot(erich, aes(x = Ammonia, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Ammonia") +
  ylab("Microbial Evenness") +
  theme_bw()

##How about adding an N:P ratio so we can plot this too?
erich$N.P <- (erich$Ammonia + erich$Nitrite_plus_Nitrate)/erich$Phosphate

p4 <- ggplot(erich, aes(x = N.P, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("N:P Ratio") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p1, p2, p3, p4) #from library "cowplot"
ggsave("../output/plots/evenness_vs_nutrients.pdf", plot = last_plot())

##Can we add the Turbinaria Nutrient data? 
turb.nut <- read.csv("../../../metadata/Turbinaria_CHN_May_2021_compiled.csv")
View(turb.nut)
keeps <- c("Site_number","Weight_ug", "Percent_C", "Percent_H", "Percent_N", "C_to_N_ratio")
turb.nut <- turb.nut[ , (names(turb.nut) %in% keeps)]
View(turb.nut)    

#Merge the turb.nut with the erich data! 
erich$Site_number <- erich$site.number
rich.turb <- merge(erich, turb.nut, by = "Site_number")    
View(rich.turb)

#Trial the same as above but with the Turb

p5 <- ggplot(rich.turb, aes(x = Percent_N, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("/%N") +
  ylab("Microbial Species Richness") +
  theme_bw()

p6 <- ggplot(rich.turb, aes(x = Percent_N, y = Shannon)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("/%N") +
  ylab("Microbial Shannon Diversity") +
  theme_bw()

p7 <- ggplot(rich.turb, aes(x = Percent_N, y = FaithPD)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("/%N") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()

p8 <- ggplot(rich.turb, aes(x = Percent_N, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("/%N") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p5, p6, p7, p8) #from library "cowplot"




