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

###rename the columns to ensure we know which nutrients came from water and which from turbs
names(meta.with.turb)[names(meta.with.turb) == 'phosphate'] <- 'water_phosphate'
names(meta.with.turb)[names(meta.with.turb) == 'silicate'] <- 'water_silicate'
names(meta.with.turb)[names(meta.with.turb) == 'nitrite_plus_nitrate'] <- 'water_nitrite_plus_nitrate'
names(meta.with.turb)[names(meta.with.turb) == 'ammonia'] <- 'water_ammonia'
names(meta.with.turb)[names(meta.with.turb) == 'Weight_ug'] <- 'turb_weight_ug'
names(meta.with.turb)[names(meta.with.turb) == 'Percent_C'] <- 'turb_percent_C'
names(meta.with.turb)[names(meta.with.turb) == 'Percent_H'] <- 'turb_percent_H'
names(meta.with.turb)[names(meta.with.turb) == 'Percent_N'] <- 'turb_percent_N'
names(meta.with.turb)[names(meta.with.turb) == 'C_to_N_ratio'] <- 'turb_C_to_N_ratio'


write.csv(meta.with.turb, "../output/metadata_with_turb.csv")
