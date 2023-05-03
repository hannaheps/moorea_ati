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

set2 <- brewer.pal(n = 5, name = "Set2")

p.ord <- plot_ordination(physeq.r.ra, bc.ord, type = "Site", color = "Habitat", title = "Bray-Curtis PCoA")
p.ord + geom_point(size = 4, aes(shape = Island_shore)) +
  scale_color_manual(values=set2) +
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



##How about looking at dispersion
##I think we want to look at bray curtis ~ nutrient regime (Turbinaria_groups_Wet)
physeq.r.ra.trim <- subset_samples(physeq.r.ra, Turbinaria_groups != "NA")
data.trim <- as(sample_data(physeq.r.ra.trim), "data.frame")
bc <- phyloseq::distance(physeq.r.ra.trim, method = "bray")
adonis2(bc ~ Turbinaria_groups, data = data.trim, method = "bray")
#adonis2(formula = bc ~ Turbinaria_groups, data = data.trim, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#Turbinaria_groups   8   2.6563 0.10503 2.5818  0.001 ***
#Residual          176  22.6343 0.89497                  
#Total             184  25.2905 1.00000  

disp.bc <- betadisper(bc, data.trim$Turbinaria_groups, type = "centroid")
pdf("../output/plots/dispersion_turbinaria_group_bray_curtis.pdf")
boxplot(disp.bc, ylab = "Distance to Centroid", xlab = "Turbinaria Quantile Group")
dev.off()

TukeyHSD(disp.bc, which = "group", ordered = FALSE,
         conf.level = 0.95)

#How about Means separately?
adonis2(bc ~ Turbinaria_quantile_meanName, data = data.trim, method = "bray")
#Turbinaria_quantile_meanName   2   1.7756 0.07021 6.8712  0.001***
#Residual                     182  23.5150 0.92979              
#Total                        184  25.2905 1.00000                       184  26.3659 1.0000      
disp.bc <- betadisper(bc, data.trim$Turbinaria_quantile_meanName, type = "centroid")
pdf("../output/plots/dispersion_turbinaria_mean_bray_curtis.pdf")
boxplot(disp.bc, ylab = "Distance to Centroid", xlab = "Turbinaria Mean")
dev.off()
TukeyHSD(disp.bc, which = "group", ordered = FALSE,
         conf.level = 0.95)
#diff         lwr         upr     p adj
#Low-High -0.10064629 -0.15365560 -0.04763698 0.0000380*
#Med-High -0.01906509 -0.07207440  0.03394422 0.6725000
#Med-Low   0.08158120  0.02814611  0.13501629 0.0011563*

##How about variance separately?
adonis2(bc ~ Turbinaria_quantile_variName, data = data.trim, method = "bray")
#                              Df SumOfSqs     R2      F Pr(>F)  
#Turbinaria_quantile_variName   2   0.4823 0.01907 1.7691  0.068
#Residual                     182  24.8082 0.98093              
#Total                        184  25.2905 1.00000                                  
disp.bc <- betadisper(bc, data.trim$Turbinaria_quantile_variName, type = "centroid")
pdf("../output/plots/dispersion_turbinaria_variance_bray_curtis.pdf")
boxplot(disp.bc, ylab = "Distance to Centroid", xlab = "Turbinaria Variance")
dev.off()


##Can we visualise how these look in ordination space?
turb.ord.plot <- plot_ordination(physeq.r.ra, bc.ord, type = "UniqueID", color = "Turbinaria_quantile_meanName", title = "Bray Curtis PCoA")
turb.ord.plot + 
  geom_point(size = 3) +
  stat_ellipse(aes(colour = Turbinaria_groups), linetype = 2) +
  theme_bw()



#we can a distance based linear model to  look at continuous variables:
##Examples -
bc.full <- phyloseq::distance(physeq.r.ra, method = "bray")
adonis2(bc.full ~ Nitrite_plus_Nitrate, data = data.r.ra, method = "bray")
adonis2(bc.full ~ HIX, data = data.r.ra, method = "bray")


