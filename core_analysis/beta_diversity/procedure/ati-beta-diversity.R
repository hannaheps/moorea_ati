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
physeq <- readRDS("../../pre_processing/output/ati-2021-physeq-nonrare.RDS")
data <- as(sample_data(physeq), "data.frame")
data$sample.id <- rownames(data)

#Using the rarefied data for alpha diversity
physeq_rare <- readRDS("../../pre_processing/output/ati-2021-physeq-rare.RDS") 
data.rare <- as(sample_data(physeq_rare), "data.frame")
data.rare$sample.id <- rownames(data.rare)

#Make relative abundance object from unrarefied
physeq_ra <- transform_sample_counts(physeq, function(x) x/ sum(x))
data.ra <- as(sample_data(physeq), "data.frame")

#rarefied relative abundances
physeq_r_ra <- transform_sample_counts(physeq_rare, function(x) x/ sum(x))
data.r.ra <- as(sample_data(physeq_r_ra), "data.frame")


#How about just look at the bar plots by sample
percent.trial <- physeq_rare %>% 
  tax_glom(taxrank = "Phylum", NArm=TRUE) %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) 
perc.melt <- psmelt(percent.trial)
sum <- ddply(perc.melt, c("Phylum", "Site"),summarise,
             N = length(Abundance), 
             mean = mean(Abundance),
             sd = sd(Abundance), 
             se = sd/sqrt(N)
)


sum.as.factor <- sum
sum.as.factor$Phylum <- as.factor(sum.as.factor$Phylum)
levels(sum.as.factor$Phylum)

nb.cols <- 66
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
ggplot(sum, aes(x = Site, y = mean, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mycolors) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("../output/plots/phylum_barplot_by_sample.pdf", plot = last_plot())



##What about by top 10 families? Probably more useful
#Something wrong with this code below, need to sort it out
physeq.rra.melt <- psmelt(physeq_r_ra)

x <-  tapply(physeq.rra.melt$Abundance, physeq.rra.melt$Family, function(x) max(x))
x <-  sort(x, TRUE) #sort phyla 
View(x) #add top 14 to code below to set NAs
#Clade I and Clade II = SAR11
physeq.rra.melt$Family <-  factor(as.character(physeq.rra.melt$Family, levels=names(x)))
physeq.rra.melt$col_family <- physeq.rra.melt$Family
View(physeq.rra.melt$col_family)
sum.fam <- as.data.frame(physeq.rra.melt)
#Set everything that is super low to NA so that we can call them "other"
sum.fam$col_family <- as.character(sum.fam$col_family)
sum.fam$col_family <- ifelse(is.na(sum.fam$col_family), 
                             'Unassigned', sum.fam$col_family)
physeq.rra.melt$col_family <- as.factor(physeq.rra.melt$col_family)

physeq.rra.melt$col_family[sum.fam$col_family != "Cryomorphaceae" &
                             sum.fam$col_family != "Cyanobiaceae" & 
                             sum.fam$col_family != "Alteromonadaceae" &
                             sum.fam$col_family != "Rhodobacteraceae" &
                             sum.fam$col_family != "Moraxellaceae" &
                             sum.fam$col_family != "Nitrincolaceae" &
                             sum.fam$col_family != "Arcobacteraceae" &
                             sum.fam$col_family != "Halomonadaceae" &
                             sum.fam$col_family != "NS9_marine_group" &
                             sum.fam$col_family != "Flavobacteriaceae" &
                             sum.fam$col_family != "Bacillaceae" &
                             sum.fam$col_family != "Litoricolaceae" &
                             sum.fam$col_family != "Clade_I" &
                             sum.fam$col_family != "Clade_II"] <- NA


levels(physeq.rra.melt$col_family)
# add new factor
physeq.rra.melt$col_family <- factor(physeq.rra.melt$col_family, levels = c(levels(physeq.rra.melt$col_family), "Other"))
# convert NAs to other
physeq.rra.melt$col_family [is.na(physeq.rra.melt$col_family)] = "Other"
physeq.rra.melt$col_family <- factor(x = physeq.rra.melt$col_family, levels = c("Cryomorphaceae", "Cyanobiaceae", "Arcobacteraceae",
                                                                                "Rhodobacteraceae", "Clade_II", "Clade_I",
                                                                                "Flavobacteriaceae", "Nitrincolaceae", "NS9_marine_group",
                                                                                "Halieaceae", "SAR116 clade", "Halomonadaceae", "Moraxellaceae", 
                                                                                "Litoricolaceae", "Bacillaceae", "Other" ))
#Make a colour scheme
ale.colors <-  c( "#80B1D3","#FFFFB3","#ABA3E6","#FB8072","#8DD3C7", "#FDB462", "#B3DE69",
                  "#FCCDE5", "#BC80BD", "#FFED6F", "#219EBC","#E06F1F","#90BE6D", "#0A9396","#9B2226", "#A3A5A8" )

#plot!
ggplot(physeq.rra.melt, aes(x = Site, y = Abundance, fill = col_family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=ale.colors) +
  ylab("Relative Abundance") +
  xlab("Sample by Site #") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("../output/plots/family_barplot_by_sample.pdf", plot = last_plot())



##Top most abundant genera
physeq.rra.melt.comb <- physeq.rra.melt
physeq.rra.melt.comb$Taxonomy <- as.factor(paste(physeq.rra.melt.comb$Family, physeq.rra.melt.comb$Genus, sep = "_")) 

x = tapply(physeq.rra.melt.comb$Abundance, physeq.rra.melt.comb$Taxonomy, function(x) max(x))
x = sort(x, TRUE) #sort by Family/Genus
View(as.data.frame(x))
View(as.data.frame(physeq.rra.melt.comb))

#summarize and pull out abundances of all genera by site
sum.gen.abund <- ddply(physeq.rra.melt.comb, c("Taxonomy", "Site"),summarise,
                       N = sum(Abundance)
)
sum.gen.abund <- sum.gen.abund %>% spread(Taxonomy, N)
print(colnames(sum.gen.abund))
write.csv(sum.gen.abund, "../output/site_taxononmy_abund.csv", row.names = FALSE)

##Can we subset this file to just the top abundant genera (greater than 5% of the community)
sum.gen.abund.top <- subset(sum.gen.abund, select = c("Cryomorphaceae_uncultured", "Cyanobiaceae_Synechococcus_CC9902",
                                                      "Alteromonadaceae_Alteromonas", "Rhodobacteraceae_HIMB11", "Cyanobiaceae_Prochlorococcus_MIT9313",
                                                      "Moraxellaceae_Acinetobacter", "Nitrincolaceae_Marinobacterium", "Arcobacteraceae_uncultured",
                                                      "Halomonadaceae_Halomonas", "NS9_marine_group_NS9_marine_group", "Flavobacteriaceae_NS5_marine_group",
                                                      "Bacillaceae_Bacillus", "Clade_I_Clade_Ia", "Litoricolaceae_Litoricola", "Moraxellaceae_Psychrobacter",
                                                      "Arcobacteraceae_NA", "Clade_II_Clade_II", "SAR86_clade_SAR86_clade", "Halieaceae_OM60(NOR5)_clade"))

sum.gen.abund.top$Site <- rownames(sum.gen.abund.top)
#Cool, let's add it to our erich dataset and make a heatmap

erich <- read.csv("../../alpha_diversity/output/ati-2021-metadata-with-alphadiv.csv")

erich.topabund <- merge(erich, sum.gen.abund.top, by = "Site")
#Remove the duplicated column
erich.topabund <- erich.topabund[, -2]
#export the metadata

write.csv(erich.topabund, "../output/metadata_plus_topabund_genera.csv", row.names = FALSE)



##################

#Below code is old, potentially deprecated.
#If revived, will move above the hashes



##Can we pull out psychrobacter relative abundance and then map that? 
psychro <- subset_taxa(physeq_r_ra, Genus == "Psychrobacter")
plot_bar(psychro)

#Make it prettier
psychro.melt <- psmelt(psychro)
ggplot(psychro.melt, aes(x = site.number, y = Abundance))+
  geom_bar(stat = "identity", fill = "turquoise") +
  xlab("\nSite Number") +
  ylab("Mean Relative Abundance\n") +
  ggtitle("Psychrobacter: Moorea Around-the-island") +
  theme_bw()

#Export the data for Tom to plot around the island
write.csv(psychro.melt, "../output/psychrobacter_relabund.csv")
as.factor(psychro.melt$OTU)


##Can we do an adonis with a continous variable? YES!
#It's called a distance-based linear model (distlm) when it's with continuous
#Remove all the samples with NAs in Percent_N
physeq.percN <- subset_samples(physeq_rare, Site_number !=  "83" & Site_number != "Cook's 1")
percN.data <- as(sample_data(physeq.percN), "data.frame")

bc <- phyloseq::distance(physeq.percN, method = "bray")
adonis(bc ~ Turb_Percent_N, data = percN.data, method = "bray")

#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Turb_Percent_N  1    0.4025 0.40252  2.3777 0.02721  0.023 *
#Residuals      85   14.3895 0.16929         0.97279         
#Total          86   14.7920                 1.00000         

bray_curtis_pcoa <- ecodist::pco(bc)
# All components could be found here: 
# bray_curtis_pcoa$vectors
# But we only need the first two to demonstrate what we can do:
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])


bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PCo1",
       y = "PCo2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller


bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df,
                                   percN = percN.data$Turb_Percent_N)


# Creates a plot
bray_curtis_percN_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = percN)) + 
  geom_point(size = 3) +
  labs(x = "PC1",
       y = "PC2",
       title = "PCoA of Microbial Communities vs  %N") +
  theme_bw()

bray_curtis_percN_plot
ggsave("../output/PCoA_percN.pdf", plot = last_plot())

