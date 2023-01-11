###Regressions using both alpha and beta div###

#use the combined erich data + pcoa axes to build regressions

#Set working directory to regression directory
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/core_analysis/regressions/procedure/")

##pull in the micro data + metadata file
erich.pcoa <- read.csv("../../beta_diversity/output/ati_2021_metadata_with_micro.csv")

#now we have all the data: micro alpha div plus all the metadata!

##Regressions!!##

hist(log(erich.pcoa$HIX))
ggplot(erich.pcoa, aes(x = log(HIX), y = pcoa1)) + geom_point(size=2) + geom_smooth(method=lm) +
  xlab("log(HIX)") + ylab("Bray-Curtis PCoA Axis 1") +
  theme_bw()
lm.pcoa.HIX <- lm(pcoa1 ~ log(HIX), erich.pcoa)
summary(lm.pcoa.HIX)
(lm.pcoa.HIX)$effects
#Multiple R-squared:  0.3119,	Adjusted R-squared:  0.3083
#F-statistic: 87.04 on 1 and 192 DF,  p-value: < 2.2e-16
plot(lm.pcoa.HIX)

hist(sqrt(erich.pcoa$FishNH4_l_min_m2))
ggplot(erich.pcoa, aes(x = sqrt(FishNH4_l_min_m2), y = pcoa1)) + 
  geom_point(size=2) + 
  geom_smooth(method=lm) +
  xlab("Sqrt(Fish NH4 Excretion)") + 
  ylab("Bray-Curtis PCoA Axis 1") +
  theme_bw()
lm.pcoa.fishexc <- lm(pcoa1 ~ FishNH4_l_min_m2, erich.pcoa)
summary(lm.pcoa.fishbio)


ggplot(erich.pcoa, aes(x = as.factor(Regime), y = pcoa1)) +
  geom_violin() +
  geom_point(aes(colour = Habitat), alpha = 0.8, position = position_jitter(0.05)) +
  xlab("") +
  ylab("Bray-Curtis PCoA Axis 1") +
  theme_bw()
ggsave("../output/plots/violin_regime.pdf", plot = last_plot())
summary(aov(pcoa1 ~ as.factor(Regime), data = erich.pcoa))
TukeyHSD(aov(pcoa1 ~ as.factor(Regime), data = erich.pcoa))


hist(log(erich.pcoa$Distance_to_crest))
ggplot(erich.pcoa, aes(x = log(Silicate), y = pcoa1)) + 
  geom_point(size=2, aes(color = Lagoon)) + 
  geom_smooth(method=loess) +
  xlab("Distance to Population Center") + 
  ylab("Bray-Curtis PCoA Axis 1") +
  theme_bw()
?geom_smooth
lm.pcoa.fishexc <- lm(pcoa1 ~ FishNH4_l_min_m2, erich.pcoa)
summary(lm.pcoa.fishbio)


