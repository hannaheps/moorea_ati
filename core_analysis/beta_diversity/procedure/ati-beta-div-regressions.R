###R Script for Beta Diversity Regressions: Around-the-Island Water Microbe Samples (16S)###

#Load all the libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(cowplot)

#set your working directory (should be in /core_analysis/beta_diversity/procedure)
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/beta_diversity/procedure/")


###Regressions###

#Let's bring that in from Nyssa's finalized dataset

ati <- read.csv("https://raw.githubusercontent.com/njsilbiger/ATI_NutrientRegimes/main/Data/AllNutrientData_clusters.csv")
erich <- subset(ati, as.factor(Year) == "2021") #I call this "erich" to maintain consistency in my code

##1. Turbinaria Nutrients

# C to N ratio
ggplot(erich, aes(x = HIX, y = Microbial_PCoA1)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("HIX") +
  ylab("Microbial PCoA1 Axis") +
  theme_bw()
summary(lm(Microbial_PCoA1 ~ HIX, erich))
#Multiple R-squared:  0.3227,	Adjusted R-squared:  0.3192
#F-statistic: 91.49 on 1 and 192 DF,  p-value: < 2.2e-16
plot(lm(Microbial_PCoA1 ~ HIX, erich))

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
p.hix.pcoa <- ggplot(erich.full, aes(x = HIX, y = pcoa1)) +
  geom_point(size=3, aes(colour  = Island_shore)) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Humification Index (HIX)") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw() +
  facet_wrap(erich.full$Habitat)

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


ggplot(erich.full, aes(x = turb_percent_N, y = MarineHumic.like)) +
  geom_point(size=3) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Percent N") +
  ylab("Marine Humic Like") +
  theme_bw()

summary(lm(MarineHumic.like ~ turb_percent_N, erich.full))

##Add Nury's regimes to this data file
nury.regime <- read.csv("../../../metadata/new_metadata/input/regime_meta.csv")
nury.regime$sample_id <- nury.regime$Site
erich.full.regime <- merge(erich.full, nury.regime, by = "sample_id", all = TRUE)
erich.full.regime$Regime <- as.factor(erich.full.regime$Regime)


ggplot(erich.full.regime, aes(x = water_nitrite_plus_nitrate, y = pcoa1)) +
  geom_point(size=3, aes(colour = erich.full.regime$Regime)) +
  #geom_text(label = erich.fdom$sample_id, nudge_x = 0.001, nudge_y = 0.001) +
  geom_smooth(method=lm) +
  xlab("Water Nitrite + Nitrate") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()


ggplot(erich.full.regime, aes(x = Regime, y = pcoa1))+
  geom_boxplot() +
  xlab("Regime Shift") +
  ylab("Microbial Bray-Curtis PCoA1 Axis") +
  theme_bw()






####Mapping

library(sf)
library(sp)
library(OpenStreetMap)
library(viridisLite)

#6 fdom parameters, mess around with their correlations with 3 microbial parameters
#relationships between fdom parameters - check for orthogonal (what has meaning in fdom bc of covariation)
#pca of fdom (normalize e.g., transform best or z score the fdom parameters before pca)



# melted data frame
# fdom indices are collapsed into two columns: index_name and index_value



# filter NAs
map_long_dat = l_dat[!is.na(map_long_dat$Latitude)]

# make spacial points data frame
map_sp = SpatialPointsDataFrame(data = map_long_dat,
                                coords = list(map_long_dat$Latitude,
                                              map_long_dat$Longitude))


# Check geographic range of sampling points
limits = c(
  min(map_long_dat$Longitude),
  min(map_long_dat$Latitude), 
  max(map_long_dat$Longitude),
  max(map_long_dat$Latitude) 
)

# define a bounding box with a small cushion around the minimum and maximum
bbox = list(
  xmin = limits[1] - 0.03,
  ymin = limits[2] - 0.04,
  xmax = limits[3] + 0.03,
  ymax = limits[4] + 0.04
)

# get basemap
sa_map <- openmap(c(bbox$ymax, bbox$xmin),
                  c(bbox$ymin, bbox$xmax),
                  type = "stamen-terrain",
                  mergeTiles = TRUE)

sa_map2 <- openproj(sa_map)


mo_map = function(an_index, outlier_n){
  
  # remove highest and lowest n samples
  ind_map_dat = map_long_dat[map_long_dat$index_name == an_index &
                               !(is.infinite(map_long_dat$log_index_value))]
  ind_map_dat[order(ind_map_dat$ log_index_value)]
  ind_map_dat = head(ind_map_dat, n = -(outlier_n))
  ind_map_dat = tail(ind_map_dat, n = -(outlier_n))
  
  # actual map code
  sa_map2_plt = OpenStreetMap::autoplot.OpenStreetMap(sa_map2)+
    geom_point(data = ind_map_dat,
               aes(x = Longitude,
                   y = Latitude,
                   color = log_index_value))+
    labs(title = paste("Log",an_index),
         x = "Lon",
         y = "Lat")+
    scale_color_gradient(low = "green", high = "red")+
    guides(fill=guide_legend(title=NULL))
  
  print(sa_map2_plt)
  
  ggsave(filename = paste0("plots/exploration/ATI/map_",an_index,".png"), plot = sa_map2_plt))

