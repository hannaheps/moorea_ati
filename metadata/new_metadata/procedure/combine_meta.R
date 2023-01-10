###Quick script to combine the sequencing metadata with the full Around-the-Island data
#The output of this script will be used in downstream scripts and is exported as both csv & txt

#Make sure you are in the "metadata/new_metadata/procedure" directory
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/metadata/new_metadata/procedure/")

#read in the sequencing metadata and the ati full dataset 
#As of Jan 2023, the full ati dataset is called "All_dat_withFish.csv"
#This dataframe will likely change in the future due to the addition of extra data
#It currently includes: all site metadata, water nutrient data, turbinaria nutrient data,
#fDOM, [microbe data, to be input], and fish biomass/excretion data
#To check for updated version, go here: https://drive.google.com/drive/folders/1rHozhz0wUJtu-h9-9YORd8IhtvFD6rSd

meta <- read.csv("../input/metadata_microbe_ati_2021.csv")
ati <- read.csv("../input/All_dat_withFish.csv")
#merge (= left join) the turbinaria with the metadata, specify all = TRUE to introduce NAs where there are no overlapping data
meta.ati <- merge(meta, ati, by = "Unique_id", all = TRUE)

#keep in mind that sample_id is not useful as a unique identifier. 
#For all analyses, use Unique_id for this purpose

#subset to 2021 only
meta.ati$Year <- as.factor(meta.ati$Year)
meta.ati.2021 <- subset(meta.ati, Year == "2021")

#Duplicate Site columns means the merge will put a Site.x and Site.y
#Rename one & remove the other:
names(meta.ati.2021)[names(meta.ati.2021) == 'Site.x'] <- 'Site'
meta.ati.2021 <- meta.ati.2021[,-7]

#We can also combine with categorical site regimes on a one year basis
#These were developed from 1) long-term Turbinaria data collected from the MCR LTER & ATI teams (called "classified_turbData.csv")
#And 2) Nury's Cluster Analysis (called "Regime_meta.csv")
#Updates can be found in the above google drive folder

turb.regime <-  read.csv("../input/classified_turbData.csv")
turb.regime$Year <- as.factor(turb.regime$Year)
turb.regime.2021 <- subset(turb.regime, Year == "2021")
nut.regime <- read.csv("../input/Regime_meta.csv")


meta.ati.2021.full <- merge(meta.ati.2021, turb.regime.2021, by = "Site", all = TRUE)
meta.ati.2021.full <- merge(meta.ati.2021.full, nut.regime, by = "Site", all = TRUE)

View(meta.ati.2021.full)
names(meta.ati.2021.full)[names(meta.ati.2021.full) == 'Year.x'] <- 'Year'
meta.ati.2021.full <- meta.ati.2021.full[,-45]

#Write to CSV and txt
write.csv(meta.ati.2021.full, "../output/metadata_ati_full_2021.csv")
write.table(meta.ati.2021.full, "../output/metadata_ati_full_2021.txt")


