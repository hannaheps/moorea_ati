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

#Using the rarefied data for alpha diversity
physeq_rare <- readRDS("../../pre_processing/output/ati-physeq-11146.RDS") 
data.rare <- as(sample_data(physeq_rare), "data.frame")
data.rare$sample.id <- rownames(data.rare)

#Make relative abundance object from unrarefied
physeq_ra <- transform_sample_counts(physeq, function(x) x/ sum(x))
data.ra <- as(sample_data(physeq), "data.frame")

#rarefied relative abundances
physeq_r_ra <- transform_sample_counts(physeq_rare, function(x) x/ sum(x))
data.r.ra <- as(sample_data(physeq_r_ra), "data.frame")


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


##JANUARY 2023
##Running the PCoA

##First up is Bray Curtis
bc <- phyloseq::distance(physeq_rare, method = "bray")
#If we want to look at relationships
#adonis(bc ~ Turb_Percent_N, data = percN.data, method = "bray")

bray_curtis_pcoa <- ecodist::pco(bc)

print(bray_curtis_pcoa)
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


ggsave("../output/plots/PCoA_bray.pdf", plot = last_plot())


##Can we split by reef location
bray_curtis_pcoa_data_df <- cbind(bray_curtis_pcoa_df,
                            data.rare)

write.csv(bray_curtis_pcoa_data_df, "../output/bray_curtis_PCoA_Axes_data.csv")
# Creates a plot
bray_curtis_N_plot <- ggplot(data = bray_curtis_pcoa_data_df, aes(x=pcoa1, y=pcoa2, color = turb_C_to_N_ratio)) + 
  geom_point(size = 3) +
  labs(x = "PC1",
       y = "PC2",
       title = "PCoA of Microbial Communities vs  C:N") +
  theme_bw()

bray_curtis_percN_plot
ggsave("../output/PCoA_percN.pdf", plot = last_plot())


###Weighted UniFrac

wuf <- phyloseq::distance(physeq_rare, method = "wunifrac")
#If we want to look at relationships
#adonis(bc ~ Turb_Percent_N, data = percN.data, method = "bray")

wunifrac_pcoa <- ecodist::pco(wuf)
print(wunifrac_pcoa)
# All components could be found here: 
# bray_curtis_pcoa$vectors
# But we only need the first two to demonstrate what we can do:
wunifrac_pcoa_df <- data.frame(pcoa1 = wunifrac_pcoa$vectors[,1], 
                               pcoa2 = wunifrac_pcoa$vectors[,2])


wunifrac_plot <- ggplot(data = wunifrac_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PCo1",
       y = "PCo2", 
       title = "Weighted UniFrac PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller


ggsave("../output/plots/PCoA_WeightedUniFrac.pdf", plot = last_plot())


wunifrac_pcoa_data_df <- cbind(wunifrac_pcoa_df,
                                  data.rare)

# Creates a plot
wunifrac_N_plot <- ggplot(data = wunifrac_pcoa_data_df, aes(x=pcoa1, y=pcoa2, color = turb_C_to_N_ratio)) + 
  geom_point(size = 3) +
  labs(x = "PCo1",
       y = "PCo2",
       title = "W UniFrac PCoA of Microbial Communities vs  C:N") +
  theme_bw()

bray_curtis_percN_plot
ggsave("../output/plots/wunifrac_PCoA_percN.pdf", plot = last_plot())


##How about nmds
ord.bc <- ordinate(physeq.r.coral, "NMDS", "bray", trymax = 500) #Run 20 stress 0.07221613
stressplot(ord.bc)
scores.bc <- as.data.frame(scores(ord))
scores.bc <- cbind(scores.bc, data.rare)


plot.br <- ggplot(scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = motu, shape = island.side), size = 4) +
  #geom_text(label = scores.outrm$sample.name) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = motu), linetype = 2) +
  theme_bw()

plot.br+ facet_grid(cols = vars(motu))












#distance-based linear model (distlm) when it's with continuous
bc <- phyloseq::distance(physeq_rare, method = "bray")
View(data.rare)
adonis(bc ~ Water_Phosphate, data = data.rare, method = "bray")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Water_Phosphate  1    0.2442 0.24417  1.4123 0.01597  0.164
#Residuals       87   15.0414 0.17289         0.98403       
#Total           88   15.2855                 1.00000         


##distlm for fishy things
fish <- read.csv("../../parrotfish_microbes/input/fish_summaries.csv")
fish$Site_number <- fish$Site
data.full <- left_join(data.rare, fish, by = "Site_number")
sample_data(physeq_rare)$Herbivore_total <- data.full$Herbivore_total
sample_data(physeq_rare)$Corallivore_total <- data.full$Corallivore_total
#Remove samples that don't have any fishy data
physeq.fish <- subset_samples(physeq_rare, Herbivore_total !=  "NA")
fish.data <- as(sample_data(physeq.fish), "data.frame")

#Let's do the adonis
bc.fish <- phyloseq::distance(physeq.fish, method = "bray")
adonis(bc.herb ~ Herbivore_total, data = fish.data, method = "bray")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Herbivore_total  1    0.2621 0.26210  1.6706 0.04819  0.096 .
#Residuals       33    5.1774 0.15689         0.95181         
#Total           34    5.4395                 1.00000  

adonis(bc.fish ~ Corallivore_total, data = fish.data, method = "bray")
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Corallivore_total  1    0.1716 0.17155  1.0747 0.03154   0.35
#Residuals         33    5.2680 0.15963         0.96846       
#Total             34    5.4395                 1.00000   

##Can we extract PCoA vector 1 & 2 to plot around the island to help Nyssa's group?
bc <- phyloseq::distance(physeq_rare, method = "bray")
bray_curtis_pcoa <- ecodist::pco(bc)
# bray_curtis_pcoa$vectors
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df,
                             Site_number = data.rare$Site_number)

write.csv(bray_curtis_pcoa_df, "../output/pcoa_vectors.csv")
