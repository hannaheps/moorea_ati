###PARROTFISH MICROBES###

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/parrotfish_microbes/procedure/")

##Read in the fish data & the micro alpha diversity data so that we can combine the two

fish <- read.csv("../input/AroundIsland_FishData_01302022.csv")
rich <- readRDS("../../alpha_diversity/output/ati-erich.RDS")

#Make the site number data column the same name 
rich$Site_number <- rich$site.number

#Merge the two datasets by Site_numbeer
fish.rich <- merge(fish, rich, by = "Site_number")
#Get rid of unuseful columns and/or repeated columns 
drops <- c("site.number","sample_type", "lat", "long", "collection_time", "collection_date", "sample.id")
fish.rich <- fish.rich[ , !(names(fish.rich) %in% drops)]
View(fish.rich)

#Add microbial evenness
fish.rich$evenness <- fish.rich$Shannon/log(fish.rich$Observed)

##Let's add the Turbinaria nutrient data that overlaps 
turb.nut <- read.csv("../../../metadata/Turbinaria_CHN_May_2021_compiled.csv")
View(turb.nut)
keeps <- c("Site_number","Weight_ug", "Percent_C", "Percent_H", "Percent_N", "C_to_N_ratio")
turb.nut <- turb.nut[ , (names(turb.nut) %in% keeps)]
View(turb.nut)    

#Merge the turb.nut with the erich data! 
fish.rich <- merge(fish.rich, turb.nut, by = "Site_number")    
View(fish.rich)

#Output the csv
write.csv(fish.rich, "../output/fish-rich.csv")


##Correlations!
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(cowplot)

##Total Fish & Microbe Alpha
p1 <- ggplot(fish.rich, aes(x = FISH_TOT, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Total Fish Abundance/Minute") +
  ylab("Microbial Species Richness") +
  theme_bw()

p2 <- ggplot(fish.rich, aes(x = FISH_TOT, y = Shannon)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Fish Abundance/Minute") +
  ylab("Microbial Shannon Diversity") +
  theme_bw()

p3 <- ggplot(fish.rich, aes(x = FISH_TOT, y = FaithPD)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Fish Abundance/Minute") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()

p4 <- ggplot(fish.rich, aes(x = FISH_TOT, y = evenness)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Fish Abundance/Minute") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p1, p2, p3, p4)
ggsave("../output/alpha_vs_fishabund.pdf", plot = last_plot())


##Total Fish & Microbe Alpha
p5 <- ggplot(fish.rich, aes(x = FISH_TOT, y = Percent_N)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Total Fish Abundance/Minute") +
  ylab("%N") +
  theme_bw()

p6 <- ggplot(fish.rich, aes(x = FISH_TOT, y = C_to_N_ratio)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Fish Abundance/Minute") +
  ylab("C/N ratio") +
  theme_bw()


plot_grid(p5, p6)
ggsave("../output/TurbN_vs_fishabund.pdf", plot = last_plot())



##Can we use the full fish dataset?
fish.turb<- merge(fish, turb.nut, by = "Site_number")    
View(fish.turb)
write.csv(fish.turb, "../output/fish_turbinaria.csv")




