###PARROTFISH MICROBES###

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/parrotfish_microbes/procedure/")

##Read in the fish data & the micro alpha diversity data so that we can combine the two

fish <- read.csv("../input/fish_summaries.csv")
rich <- readRDS("../../alpha_diversity/output/ati-erich.RDS")

#Make the site number data column the same name 
rich$Site <- rich$site.number

#Merge the two datasets by Site_numbeer
fish.rich <- merge(fish, rich, by = "Site")
#Get rid of unuseful columns and/or repeated columns 
drops <- c("site.number","Site_name", "sample_type", "lat", "long", "collection_time", "collection_date", "sample.id")
fish.rich <- fish.rich[ , !(names(fish.rich) %in% drops)]
View(fish.rich)

#Output the csv
write.csv(fish.rich, "../output/fish-rich.csv")


##Correlations!
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(cowplot)

##Total Herbivore/Min & Microbe Alpha
p1 <- ggplot(fish.rich, aes(x = Herbivore_total, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Total Herbivorous Fish Counts/Minute") +
  ylab("Microbial Species Richness") +
  theme_bw()

p2 <- ggplot(fish.rich, aes(x = Herbivore_total, y = Shannon)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Herbivorous Fish Counts/Minute") +
  ylab("Microbial Shannon Diversity") +
  theme_bw()

p3 <- ggplot(fish.rich, aes(x = Herbivore_total, y = FaithPD)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Herbivorous Fish Counts/Minute") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()

p4 <- ggplot(fish.rich, aes(x = Herbivore_total, y = evenness)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Herbivorous Fish Counts/Minute") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p1, p2, p3, p4)
ggsave("../output/alpha_vs_herbabund.pdf", plot = last_plot())



##Total Corallivores/Min & Microbe Alpha
#Remove the weird outlier @site 173
fish.rich.outrm <- subset(fish.rich, Site !="173")
p5 <- ggplot(fish.rich.outrm, aes(x = Corallivore_total, y = Observed)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Total Corallivorous Fish Counts/Minute") +
  ylab("Microbial Species Richness") +
  theme_bw()

p6 <- ggplot(fish.rich.outrm, aes(x = Corallivore_total, y = Shannon)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Corallivorous Fish Counts/Minute") +
  ylab("Microbial Shannon Diversity") +
  theme_bw()

p7 <- ggplot(fish.rich.outrm, aes(x = Corallivore_total, y = FaithPD)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Corallivorous Fish Counts/Minute") +
  ylab("Microbial Phylogenetic Diversity") +
  theme_bw()

p8 <- ggplot(fish.rich.outrm, aes(x = Corallivore_total, y = evenness)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Corallivorous Fish Counts/Minute") +
  ylab("Microbial Evenness") +
  theme_bw()

plot_grid(p5, p6, p7, p8)
ggsave("../output/alpha_vs_corallabund_outrm.pdf", plot = last_plot())


##Total Fish & Turbinaria nutrient correlations
p9 <- ggplot(fish.rich, aes(x = Herbivore_total, y = Percent_N)) +
  geom_point(size=2) +
  geom_smooth(method=lm) +
  xlab("Total Herbivore Abundance/Minute") +
  ylab("%N") +
  theme_bw()

p10 <- ggplot(fish.rich, aes(x = Corallivore_total, y = Percent_N)) +
  geom_point(size = 2) +
  geom_smooth(method = lm) +
  xlab("Total Coralllivorous Fish Abundance/Minute") +
  ylab("%N") +
  theme_bw()


plot_grid(p9, p10)
ggsave("../output/TurbN_vs_fishabund.pdf", plot = last_plot())

##Can we use the full fish dataset?
fish.turb<- merge(fish, turb.nut, by = "Site_number")    
View(fish.turb)
write.csv(fish.turb, "../output/fish_turbinaria.csv")




