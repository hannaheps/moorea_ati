###Quick script to combine the Turbinaria Nutrient Data with the rest of the metadata for microbes

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/metadata/procedure/")

#read in the metadata and the turbinaria data
meta <- read.csv("../metadata_ati.csv")
turb <- read.csv("../Turbinaria_CHN_May_2021_compiled.csv")
#Separate out the columns we want to add
keeps <- c("Site_number","Weight_ug", "Percent_C", "Percent_H", "Percent_N", "C_to_N_ratio")
turb.nut <- turb[ , (names(turb) %in% keeps)]
#merge (= left join) the turbinaria with the metadata
meta.with.turb <- merge(meta, turb.nut, by = "Site_number")

write.csv(meta.with.turb, "../metadata_with_turb.csv")
