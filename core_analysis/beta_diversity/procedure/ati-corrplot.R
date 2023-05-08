#R script for correlation matrix. Around-the-island 16S water microbe samples - temp 
library(reshape2)
library(dplyr)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(pals)
library(colorRamps)

all.df <- read_csv("https://raw.githubusercontent.com/njsilbiger/ATI_NutrientRegimes/main/Data/AllNutrientData_clusters.csv")

###this is all v all numerical, but can clean up/arrange when we decide what features we want (or if we want...?)

comm.micro.na <- subset(all.df, select = c("Microbial_PCoA1", "Microbial_PCoA2", "Microbial_Species_Richness", "Microbial_Shannon_Diversity", "Microbial_Phylogenetic_Diversity", "Microbial_Evenness", 
"C_to_N_ratio", "Percent_C", "Percent_N", "Phosphate", "Silicate", "Nitrite_plus_Nitrate", "Ammonia", "Turbinaria_N", "Turbinaria_V", "Turbinaria_N_Wet", "Turbinaria_V_Wet", "Nutrient_PC1", "Nutrient_PC2", 
"Ultra.Violet.Humic.like", "Tyrosine.like", "Visible.Humic.like", "Marine.Humic.like", "Tryptophan.like", "Lignin.like", "Distance_to_crest", "Distance_to_shore", "Distance_to_pass", 
"Distance_to_deep_lagoon_water", "Distance_to_population_center"))

comm.micro = na.omit(comm.micro.na)



micro.cor <- cor(comm.micro)
micro.rcorr <- rcorr(as.matrix(comm.micro)) #didn't specify type = c("pearson","spearman")) but lmk
micro.coeff = as.matrix(micro.rcorr$r)
	#micro.coeff
micro.P = as.matrix(micro.rcorr$P)
	#micro.P


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

#cleaned table of correlation/significance
micro_hyptest <-flattenCorrMatrix(micro.coeff, micro.P)
print(micro_hyptest$micro.P)


#####annoying amount of temp figure options

###base - all vs. all
#hierarchical clustering
	pdf("ati_corrplot_hierarchical_f1.pdf") 
	corrplot(micro.coeff, p.mat = micro_hyptest$micro.P, sig.level = 0.05, insig = "blank", order="hclust", outline = TRUE, type="lower", tl.col = "black", tl.cex = 0.8, diag = FALSE, col1 = c)
	dev.off()
#group-by-type
	pdf("ati_corrplot_groupbytype_f2.pdf") 
	corrplot(micro.coeff, p.mat = micro_hyptest$micro.P, sig.level = 0.05, insig = "blank", order='original', outline = TRUE, type="lower", tl.col = "black", tl.cex = 0.8, diag = FALSE, col1 = c)
	dev.off()

###remove all variables except micro v. all
#group-by-type cropped
	pdf("ati_corrplot_groupbytype_f3.pdf") 
	corrplot(micro.coeff[,1:6], order='original', p.mat = micro_hyptest$micro.P[,1:6], sig.level = 0.05, insig = "blank", outline = TRUE, type="lower", tl.col = "black", tl.cex = 0.8, cl.pos="r", cl.ratio = 0.4, diag = FALSE)
	dev.off()

### recolor using pub scheme
p1 <- cubicl(255)
cubebasis <- pal.compress(cubicl)
#hierarchical x recolored 
	pdf("ati_corrplot_hierarchical_f4.pdf") 
	corrplot(micro.coeff, p.mat = micro_hyptest$micro.P, sig.level = 0.05, insig = "blank", order="hclust", outline = TRUE, type="lower", tl.col = "black", tl.cex = 0.8, diag = FALSE, col=colorRampPalette(cubebasis)(255))
	dev.off()
#group-by-type x recolored
	pdf("ati_corrplot_groupbytype_f5.pdf") 
	corrplot(micro.coeff, p.mat = micro_hyptest$micro.P, sig.level = 0.05, insig = "blank", order='original', outline = TRUE, type="lower", tl.col = "black", tl.cex = 0.8, diag = FALSE, col=colorRampPalette(cubebasis)(255))
	dev.off()
#group-by-type x cropped - recolored
	pdf("ati_corrplot_groupbytype_f6.pdf") 
	corrplot(micro.coeff[,1:6], order='original', p.mat = micro_hyptest$micro.P[,1:6], sig.level = 0.05, insig = "blank", outline = TRUE, type="lower", tl.col = "black", tl.cex = 0.8, cl.pos="r", cl.ratio = 0.4, diag = FALSE, col=colorRampPalette(cubebasis)(255))
	dev.off()
