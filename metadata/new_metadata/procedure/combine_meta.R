###Quick script to combine the Turbinaria Nutrient Data with the rest of the metadata for microbes

setwd("~/Documents/OSUDocs/Projects/French_Polynesia/Around_the_island/moorea_ati/metadata/new_metadata/procedure/")

#read in the metadata and the turbinaria data
meta <- read.csv("../input/metadata_ati.csv")
turb <- read.csv("../input/Turbinaria_CHN_May_2021_compiled.csv")
#Separate out the columns we want to add
keeps <- c("Site_number","Weight_ug", "Percent_C", "Percent_H", "Percent_N", "C_to_N_ratio")
turb.nut <- turb[ , (names(turb) %in% keeps)]
#merge (= left join) the turbinaria with the metadata
meta.with.turb <- merge(meta, turb.nut, by = "Site_number", all = TRUE)
#For some reason the above code duplicates sample from site 80 (line 171)
#so I am removing the row as follows & will check back on code later

meta.with.turb <- meta.with.turb[-171,]

write.csv(meta.with.turb, "../output/metadata_with_turb.csv")
