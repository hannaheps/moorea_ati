###Quick script to streamline sequencing data with the full Around-the-Island data 
#They currently do not meld due to the use of differing unique ids during sequencing
#The output of this script will be used in downstream scripts and is exported as both csv & txt

#First bring in the metadata from Nyssa Silbiger's github


#Make sure you are in the "metadata/new_metadata/procedure" directory
setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/metadata/new_metadata/procedure/")

#read in the sequencing metadata and the ati full dataset 
#This dataframe will likely change in the future due to the addition of extra data
#It currently includes: all site metadata, water nutrient data, turbinaria nutrient data,
#fDOM, [microbe data, to be input], and fish biomass/excretion data
#To check for updated version, go here: https://drive.google.com/drive/folders/1rHozhz0wUJtu-h9-9YORd8IhtvFD6rSd

meta <- read.csv("../input/metadata_microbe_ati_2021.csv")
ati <- read.csv("https://raw.githubusercontent.com/njsilbiger/ATI_NutrientRegimes/main/Data/AllNutrientData_clusters.csv")

#subset ati to 2021

ati.2021 <- subset(ati, as.factor(Year) == "2021")
#merge (= left join) the ati with the metadata, specify all = TRUE to introduce NAs where there are no overlapping data
meta.ati.2021 <- merge(meta, ati.2021, by = "Site", all = TRUE)
View(meta.ati.2021)

#Remove duplicate columns & fix names
meta.ati.2021 <- meta.ati.2021[,-3]
rownames <- meta.ati.2021$UniqueID

#Write to CSV and txt for use in downstream processes
write.csv(meta.ati.2021, "../output/metadata_ati_full_2021.csv", row.names = FALSE)
write.table(meta.ati.2021, "../output/metadata_ati_full_2021.txt", row.names = FALSE, sep="\t")
##Check these outputs for index columns & row name corrections

