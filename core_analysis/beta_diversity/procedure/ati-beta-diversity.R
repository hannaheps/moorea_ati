###R Script for Beta Diversity: Around-the-island Water Microbe Samples (16S)###

#Load all the libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)
library(gridExtra)
library(vegan)
library(RColorBrewer)

## Set working directory
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/beta_diversity/procedure/")

#Read in the phyloseq object (output from pre-processing R script)
physeq <- readRDS("../../pre_processing/output/ati-physeq.RDS")
data <- as(sample_data(physeq), "data.frame")
##For now I am going to add the Turbinaria nutrient data here, but eventually this needs to be aded
#to the metadata and put into the pre-processing script
turb.nut <- read.csv("../../../metadata/Turbinaria_CHN_May_2021_compiled.csv")
View(turb.nut)
keeps <- c("Site_number","Weight_ug", "Percent_C", "Percent_H", "Percent_N", "C_to_N_ratio")
turb.nut <- turb.nut[ , (names(turb.nut) %in% keeps)]
View(turb.nut)  
#Combine data + turb
data.turb <- merge(data, turb.nut, by = "Site_number")    
View(data.turb)

#Using the rarefied data for alpha diversity
physeq_rare <- readRDS("../../pre_processing/output/ati-physeq-1000.RDS") 
data.rare <- as(sample_data(physeq_rare), "data.frame")
#Add turb
data.rare.turb <- merge(data.rare, turb.nut, by = "Site_number")    
View(data.rare.turb)


#Make relative abundance object from unrarefied
physeq_ra <- transform_sample_counts(physeq, function(x) x/ sum(x))
data.ra <- as(sample_data(physeq), "data.frame")
data.ra.turb <- merge(data.ra, turb.nut, by = "Site_number")
#rarefied relative abundances
physeq_r_ra <- transform_sample_counts(physeq_rare, function(x) x/ sum(x))
data.r.ra <- as(sample_data(physeq_r_ra), "data.frame")
data.rra.turb <- merge(data.r.ra, turb.nut, by = "Site_number")


#How about just look at the bar plots by sample, since there are only 89. 
percent.trial <- physeq_rare %>% 
  tax_glom(taxrank = "Phylum", NArm=TRUE) %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) 
perc.melt <- psmelt(percent.trial)
sum <- ddply(perc.melt, c("Phylum", "site.number"),summarise,
             N = length(Abundance), 
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)


sum.as.factor <- sum
sum.as.factor$Phylum <- as.factor(sum.as.factor$Phylum)
levels(sum.as.factor$Phylum)

nb.cols <- 34
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
ggplot(sum, aes(x = site.number, y = mean, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mycolors) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("../output/plots/phylum_barplot_by_sample.pdf", plot = last_plot())


##What about by top 10 families? Probably more useful
#Something wrong with this code below, need to sort it out
physeq.rra.melt <- psmelt(physeq_r_ra)

x = tapply(physeq.rra.melt$Abundance, physeq.rra.melt$Family, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x) #add top 14 to code below to set NAs
#Clade I and Clade II = SAR11
physeq.rra.melt$Family <-  factor(as.character(physeq.rra.melt$Family, levels=names(x)))
physeq.rra.melt$col_family <- physeq.rra.melt$Family
View(physeq.rra.melt$col_family)
#Set everything that is super low to NA so that we can call them "other"
#sum.fam$col_family <- as.character(sum.fam$col_family)
#sum.fam$col_family <- ifelse(is.na(sum.fam$col_family), 
                            #'Unassigned', sum.fam$col_family)
physeq.rra.melt$col_family <- as.factor(physeq.rra.melt$col_family)

physeq.rra.melt$col_family[sum.fam$col_family != "Cryomorphaceae" &
                    sum.fam$col_family != "Cyanobiaceae" & 
                    sum.fam$col_family != "Arcobacteraceae" &
                    sum.fam$col_family != "Rhodobacteraceae" &
                    sum.fam$col_family != "Clade II" &
                    sum.fam$col_family != "Clade I" &
                    sum.fam$col_family != "Flavobacteriaceae" &
                    sum.fam$col_family != "Nitrincolaceae" &
                    sum.fam$col_family != "NS9 marine group" &
                    sum.fam$col_family != "Halieaceae" &
                    sum.fam$col_family != "SAR116 clade" &
                    sum.fam$col_family != "Halomonadaceae" &
                    sum.fam$col_family != "Moraxellaceae" &
                    sum.fam$col_family != "Actinomarinaceae" &
                    sum.fam$col_family != "AEGEAN-169 marine group"] <- NA


levels(physeq.rra.melt$col_family)
# add new factor
physeq.rra.melt$col_family <- factor(physeq.rra.melt$col_family, levels = c(levels(physeq.rra.melt$col_family), "Other"))
# convert NAs to other
physeq.rra.melt$col_family [is.na(physeq.rra.melt$col_family)] = "Other"
physeq.rra.melt$col_family <- factor(x = physeq.rra.melt$col_family, levels = c("Cryomorphaceae", "Cyanobiaceae", "Arcobacteraceae",
                                                                                "Rhodobacteraceae", "Clade II", "Clade I",
                                                              "Flavobacteriaceae", "Nitrincolaceae", "NS9 marine group",
                                                              "Halieaceae", "SAR116 clade", "Halomonadaceae", "Moraxellaceae", 
                                                              "Actinomarinaceae", "AEGEAN-169 marine group", "Other" ))
#Make a colour scheme
ale.colors <-  c( "#80B1D3","#FFFFB3","#ABA3E6","#FB8072","#8DD3C7", "#FDB462", "#B3DE69",
                  "#FCCDE5", "#BC80BD", "#FFED6F", "#219EBC","#E06F1F","#90BE6D", "#0A9396","#9B2226", "#A3A5A8" )

#plot!
ggplot(physeq.rra.melt, aes(x = site.number, y = Abundance, fill = col_family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=ale.colors) +
  ylab("Relative Abundance") +
  xlab("Sample by Site #") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("../output/plots/family_barplot_by_sample.pdf", plot = last_plot())




write.csv(x, "../output/genus_by_relabund.csv")


physeq.rra.melt <- psmelt(physeq_r_ra)
x = tapply(physeq.rra.melt$Abundance, physeq.rra.melt$Genus, function(x) max(x))
x = sort(x, TRUE) #sort phyla 
View(x)
write.csv(x, "../output/genus_by_relabund.csv")


##Can we pull out psychrobacter relative abundance and then map that? 
psychro <- subset_taxa(physeq_r_ra, Genus == "Psychrobacter")
plot_bar(psychro)

#Make it prettier
psychro.melt <- psmelt(psychro)
ggplot(psychro.melt, aes(x = site.number, y = Abundance))+
  geom_bar(stat = "identity", fill = "turquoise") +
  xlab("\nSite Number") +
  ylab("Mean Relative Abundance\n") +
  ggtitle("Psychrobacter: Mo'orea Around-the-island") +
  theme_bw()

#Export the data for Tom to plot around the island
write.csv(psychro.melt, "../output/psychrobacter_relabund.csv")
as.factor(psychro.melt$OTU)


##Can we do an adonis with a continous variable? 

bc <- phyloseq::distance(physeq_rare, method = "bray")
adonis(bc ~ Percent_N, data = data.rare.turb)
