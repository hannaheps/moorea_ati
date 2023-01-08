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



###Regressions###

#First let's bring in the PCoA data (from beta diversity directory) and include it into erich
erich <- readRDS("../output/ati-erich.RDS")
pcoa_df <- read.csv("../../beta_diversity/output/bray_curtis_PCoA_Axes_data.csv")
erich$pcoa1 <- pcoa_df$pcoa1  

##1. Turbinaria Nutrients

# C to N ratio
p.cton.obs <- ggplot(erich, aes(x = turb_C_to_N_ratio, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Turbinaria C:N") +
  ylab("Microbial Species Richness") +
  theme_bw()
lm.cton <- lm(Observed ~ turb_C_to_N_ratio, erich)
summary(lm.cton)
#Multiple R-squared:  0.01361,	Adjusted R-squared:  0.008338 
#F-statistic: 2.581 on 1 and 187 DF,  p-value: 0.1099
plot(lm.cton)

p.cton.even <- ggplot(erich, aes(x = turb_C_to_N_ratio, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Turbinaria C:N") +
  ylab("Microbial Evenness") +
  theme_bw()
lm.cton.even <- lm(evenness ~ turb_C_to_N_ratio, erich)
summary(lm.cton.even)
#Multiple R-squared:  0.01588,	Adjusted R-squared:  0.01061 
#F-statistic: 3.017 on 1 and 187 DF,  p-value: 0.08405
plot(lm.cton.even)

p.cton.pcoa <- ggplot(erich, aes(x = turb_C_to_N_ratio, y = pcoa1)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Turbinaria C:N") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
lm.cton.pcoa <- lm(pcoa1 ~ turb_C_to_N_ratio, erich)
summary(lm.cton.pcoa)
#Multiple R-squared:  0.07527,	Adjusted R-squared:  0.07033 
#F-statistic: 15.22 on 1 and 187 DF,  p-value: 0.0001333
plot(lm.cton.pcoa)
  

p.cton.fpd <- ggplot(erich, aes(x = turb_C_to_N_ratio, y = FaithPD)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Turbinaria C:N") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()
lm.cton.fpd <- lm(FaithPD ~ turb_C_to_N_ratio, erich)
summary(lm.cton.fpd)
#Multiple R-squared:  0.0192,	Adjusted R-squared:  0.01395 
#F-statistic:  3.66 on 1 and 187 DF,  p-value: 0.05725
plot(lm.cton.fpd)

plot_grid(p.cton.obs, p.cton.even, p.cton.fpd, p.cton.pcoa)
ggsave("../output/plots/C_to_N_microbial_corr.pdf", plot = last_plot())

##2. Water Nutrients

#Nitrite + Nitrate
p.nn.obs <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Water Nitrite + Nitrate") +
  ylab("Microbial Species Richness") +
  theme_bw()
lm.nn.obs <- lm(Observed ~ water_nitrite_plus_nitrate, erich)
summary(lm.nn.obs)
#Multiple R-squared:  0.1674,	Adjusted R-squared:  0.163 
#F-statistic:    38 on 1 and 189 DF,  p-value: 4.185e-09


p.nn.even <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = evenness)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Water Nitrite + Nitrate") +
  ylab("Microbial Evenness") +
  theme_bw()
lm.nn.even <- lm(evenness ~ water_nitrite_plus_nitrate, erich)
summary(lm.nn.even)
#Multiple R-squared:  0.01681,	Adjusted R-squared:  0.01161 
#F-statistic: 3.231 on 1 and 189 DF,  p-value: 0.07384


p.nn.pcoa <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = pcoa1)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Water Nitrite + Nitrate") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
lm.nn.pcoa <- lm(pcoa1 ~ water_nitrite_plus_nitrate, erich)
summary(lm.nn.pcoa)
#Multiple R-squared:  0.2525,	Adjusted R-squared:  0.2486 
#F-statistic: 63.85 on 1 and 189 DF,  p-value: 1.292e-13
plot(lm.nn.pcoa)

p.nn.fpd <- ggplot(erich, aes(x = water_nitrite_plus_nitrate, y = FaithPD)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Water Nitrite + Nitrate") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()
lm.nn.fpd <- lm(FaithPD ~ water_nitrite_plus_nitrate, erich)
summary(lm.nn.fpd)
#Multiple R-squared:  0.2105,	Adjusted R-squared:  0.2063 
#F-statistic: 50.38 on 1 and 189 DF,  p-value: 2.48e-11
plot(lm.nn.pcoa)

plot_grid(p.nn.obs, p.nn.even, p.nn.fpd, p.nn.pcoa)
ggsave("../output/plots/nitrite_nitrate_microbial_corr.pdf", plot = last_plot())



##Silicate
p.sil.pcoa <- ggplot(erich, aes(x = log(water_silicate), y = pcoa1)) +
  geom_point(size=2) +
  #geom_text(label = erich$sample_id, nudge_x = 0.005, nudge_y = 0.005) +
  geom_smooth(method=lm) +
  xlab("log(Water Silicate)") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
lm.sil.pcoa <- lm(pcoa1 ~ log(water_silicate), erich)
summary(lm.sil.pcoa)
#Multiple R-squared:  0.3244,	Adjusted R-squared:  0.3208 
#F-statistic: 90.74 on 1 and 189 DF,  p-value: < 2.2e-16
plot(lm.sil.pcoa)



###Can we add fDOM marine Humic-like things to our regressions?
#erich <- readRDS("../output/ati-erich.RDS")
#pcoa_df <- read.csv("../../beta_diversity/output/bray_curtis_PCoA_Axes_data.csv")
#erich$pcoa1 <- pcoa_df$pcoa1  
fdom_df <- read.csv("../../../metadata/new_metadata/input/fDOM_may2021.proc.csv")
fdom_df$sample_id <- fdom_df$SampleName
erich.fdom <- merge(erich, fdom_df, by = "sample_id", all = TRUE)

p.humic.pcoa <- ggplot(erich.fdom, aes(x = MarineHumic.like, y = pcoa1)) +
  geom_point(size=2) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Marine Humic Like") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
lm.humic.pcoa <- lm(pcoa1 ~ MarineHumic.like, erich.fdom)
summary(lm.humic.pcoa)
#Multiple R-squared:  0.1518,	Adjusted R-squared:  0.1473 
#F-statistic: 33.66 on 1 and 188 DF,  p-value: 2.743e-08
plot(lm.humic.pcoa)

plot_grid(p.humic.pcoa, p.sil.pcoa)
ggsave("../output/plots/silicante_humic_microbial_corr_label.pdf", plot = last_plot())
ggsave("../output/plots/silicate_humic_microbial_corr.pdf", plot = last_plot())


##other fDOM measures Linda said were interesting = HIX (humification index) & BIX (biological index), M.C (Peak C to Peak M ratio)
p.hix.pcoa <- ggplot(erich.fdom, aes(x = HIX, y = pcoa1)) +
  geom_point(size=2) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Humification Index (HIX)") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
lm.hix.pcoa <- lm(pcoa1 ~ HIX, erich.fdom)
summary(lm.hix.pcoa)
#Multiple R-squared:  0.3335,	Adjusted R-squared:  0.3299 
#F-statistic: 94.06 on 1 and 188 DF,  p-value: < 2.2e-16
plot(lm.hix.pcoa)

p.bix.pcoa <- ggplot(erich.fdom, aes(x = BIX, y = pcoa1)) +
  geom_point(size=2) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Biological Index (BIX)") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
lm.bix.pcoa <- lm(pcoa1 ~ BIX, erich.fdom)
summary(lm.bix.pcoa)
#Multiple R-squared:  0.008266,	Adjusted R-squared:  0.002991 
#F-statistic: 1.567 on 1 and 188 DF,  p-value: 0.2122
plot(lm.bix.pcoa)

p.mc.pcoa <- ggplot(erich.fdom, aes(x = M.C, y = pcoa1)) +
  geom_point(size=2) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Ratio of M to C (MC)") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
lm.mc.pcoa <- lm(pcoa1 ~ M.C, erich.fdom)
summary(lm.mc.pcoa)
#Multiple R-squared:  0.2643,	Adjusted R-squared:  0.2604 
#F-statistic: 67.53 on 1 and 188 DF,  p-value: 3.294e-14
plot(lm.mc.pcoa)

plot_grid(p.humic.pcoa, p.hix.pcoa, p.bix.pcoa, p.mc.pcoa)
ggsave("../output/plots/fdom_microbial_corr.pdf", plot = last_plot())



##Add site data?
erich.fdom
site.md <- read.csv("../../../metadata/new_metadata/input/ati_site_metadata.csv")
erich.full <- merge(erich.fdom, site.md, by = "sample_id", all = TRUE)

##Trial colouring dots by reef habitat type

ggplot(erich.full, aes(x = water_nitrite_plus_nitrate, y = pcoa1)) +
  geom_point(size=3, aes(colour = erich.full$Habitat)) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Water Nitrite + Nitrate") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()
ggsave("../output/plots/nn_microbial_corr_byhabitat.pdf", plot = last_plot())

write.csv(erich.full, "../output/ati-metadata-with-alphadiv-fDOM-site.csv" )
erich.full <- read.csv("../output/ati-metadata-with-alphadiv-fDOM-site.csv")

