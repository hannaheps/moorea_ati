###Extraction of PCoA Axes for use as a consensus microbiome measurement##
##Combo of alpha and beta diversity metrics to use as a covariate in larger Structural Equation Model##

#Load Libraries
#Load all the libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(RColorBrewer)

#Set working directory to beta_diversity/procedures
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/beta_diversity/procedure/")

##pull in the either the "erich" output from the alpha diversity script
#erich <- readRDS("../../alpha_diversity/output/ati-erich.RDS")
##OR
##pull in just the metadata from Nyssa Silbiger's github because it already has all the 
##alpha div metrics:
ati <- read.csv("https://raw.githubusercontent.com/njsilbiger/ATI_NutrientRegimes/main/Data/AllNutrientData_clusters.csv")
erich <- subset(ati, as.factor(Year) == "2021") #I call this "erich" to maintain consistency in my code


#Read in the phyloseq object (output from pre-processing R script)
#Using the rarefied data for alpha diversity
physeq_rare <- readRDS("../../pre_processing/output/ati-2021-physeq-rare.RDS") 

#Make into relative abundance to run ordinations
physeq.r.ra <- transform_sample_counts(physeq_rare, function(x) x/ sum(x))
data.r.ra <- as(sample_data(physeq.r.ra), "data.frame")

##Run a bray curtis pcoa and plot
bc.ord <- phyloseq::ordinate(physeq.r.ra, "PCoA", "bray")

p.ord <- plot_ordination(physeq.r.ra, bc.ord, type = "Site", color = "Habitat", title = "Bray-Curtis PCoA")
p.ord + geom_point(size = 4, aes(shape = Island_shore)) +
  theme_bw()
ggsave("../output/plots/bray_curtis_pcoa.pdf", plot = last_plot())

##Let's pull out axis 1 and 2 (which explian 41.7% of the variance and 9.8% of the variance, respectively)
bc_pcoa_df <- data.frame(pcoa1 = bc.ord$vectors[,1], 
                          pcoa2 = bc.ord$vectors[,2])

erich.pcoa <- cbind(erich, bc_pcoa_df)
View(erich.pcoa)

write.csv(erich.pcoa, "../output/ati_2021_metadata_with_micro.csv")
#erich.pcoa <- read.csv("../output/ati_2021_metadata_with_micro.csv")

##Can we run a UniFrac ordination? 

wuf.ord <- phyloseq::ordinate(physeq.r.ra, "PCoA", "wunifrac")

p.ord.wuf <- plot_ordination(physeq.r.ra, wuf.ord, type = "Site", color = "Habitat", title = "Weighted UniFrac PCoA")
p.ord.wuf + geom_point(size = 4, aes(shape = Island_shore)) +
  theme_bw()
ggsave("../output/plots/weighted_unifrac_pcoa.pdf", plot = last_plot())

##What if we tax_glom to genus and then run the pcoa?
physeq_genus <- physeq_rare %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) 

bc.ord.gen <- phyloseq::ordinate(physeq_genus, "PCoA", "bray")

p.ord.gen <- plot_ordination(physeq_genus, bc.ord.gen, type = "Site", color = "Habitat", title = "Bray-Curtis PCoA")
p.ord.gen + geom_point(size = 4, aes(shape = Island_shore)) +
  theme_bw()
ggsave("../output/plots/bray_curtis_bygenus_pcoa.pdf", plot = last_plot())
#This looks almost the same, but the axis 1 and 2 explain slightly more of the variance 


