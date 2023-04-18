# R script for Beta Diversity: Around-the-island 16S water microbe samples #
# map pcoa, richness - code and shapefiles derived from T. Adams, M. Donovan #
library(ggplot2)
library(maptools)
library(reshape2)
library(PBSmapping)
library(sp)
library(mapproj)
library(plyr)
library(data.table)
library(rgeos)
library(ggmap)
library(rgdal)
library(spatstat)
library(vegan)
library(GGally)
library(raster)
library(rgdal)
library(dplyr)
library(sf)
library(kriging)
library(ggnewscale)
library(ggspatial)
library(colorRamps)
library(viridis)
library(tidyverse)
library(pals)
library(svglite)

setwd("~/Around_the_island/moorea_ati/core_analysis/beta_diversity/procedure/")
all.df <- read_csv("https://raw.githubusercontent.com/njsilbiger/ATI_NutrientRegimes/main/Data/AllNutrientData_clusters.csv") 
comm.micro <- subset(all.df, select = c("UniqueID", "Microbial_Species_Richness", "Microbial_PCoA1", "Microbial_PCoA2", "Latitude", "Longitude")) #set dataframe (df)
comm.micro <- na.omit(comm.micro) #clean up df 
  #comm.micro #check df
  #dim(comm.micro) #get dimensions 

########## layer raster 
		#external shapefile or other: dtm_merged_5m.tif / dem_land.tif / Isle_outline_gcs84.shp / nut_boundary2.csv 

dem <- raster('./TM/dtm_merged_5m.tif')
dem_land <- raster('./TM/dem_land.tif') 
	#dem_land <- raster("dem_land.tif")
mo.shp <- readOGR('./TM/Isle_outline_gcs84.shp')
mo.shp.dem <- spTransform(mo.shp,projection(dem_land)) # transform to same projection as dem
mo.sf <- st_as_sf(mo.shp.dem)

## add land
  ## add hill shade
slp <- terrain(dem_land,opt='slope') 
asp <- terrain(dem_land,opt='aspect')
hill <- hillShade(slp,asp)# compute hillshade 
plot(hill)

## transform rasters for ggplot
dem.p  <-  rasterToPoints(dem_land) #conform raster to points 
dem.df <-  data.frame(dem.p) 
colnames(dem.df) = c("x", "y", "alt") #place into the df

hill.p  <-  rasterToPoints(hill) 
hill.df <-  data.frame(hill.p)
colnames(hill.df) = c("x", "y", "alt") #same - place hill points into df

########## add coordinants to UTM
spdf<-SpatialPointsDataFrame(coords = comm.micro[c('Longitude','Latitude')], data = comm.micro, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
nutsUTM<-spTransform(spdf, projection(mo.shp.dem))
nutsUTM.df <- as.data.frame(nutsUTM) #add nutrient data into lat/long 
nutsUTM.df #add nutrient data into lat/long 

ashapem <- read.csv("./TM/nut_boundary2.csv", header=T)  #nut_boundary2.csv file is housed by Tom (?) = sites where nutrient data is sampled (?) [nut/lat/long]
ashapem <- as.matrix(ashapem[,2:3]) #borderpolygon object is a list that will be part of the kriging function
borderpolygon <- list(data.frame(ashapem[,1], ashapem[,2]))
ashapem #vis ashapem
borderpolygon #vis 

########## krig

krig1 <- kriging(nutsUTM.df$Longitude.1, nutsUTM.df$Latitude.1, nutsUTM.df$Microbial_PCoA1, pixels=1000,polygons=borderpolygon) #pcoa1
  #str(krig1)
krig2 <- krig1$map
  #head(krig2)

krig3 <- kriging(nutsUTM.df$Longitude.1, nutsUTM.df$Latitude.1, nutsUTM.df$Microbial_PCoA2, pixels=1000,polygons=borderpolygon) #pcoa2
krig4 <- krig3$map

krig5 <- kriging(nutsUTM.df$Longitude.1, nutsUTM.df$Latitude.1, nutsUTM.df$Microbial_Species_Richness, pixels=1000,polygons=borderpolygon) #richness
krig6 <- krig5$map


########## PCoA1

pcoa1_1 <- ggplot() + theme_void() +	geom_point(data=krig2, aes(x=x, y=y, colour=pred), size=4) + 
  scale_color_viridis(begin = 0, end = 1, direction = 1, discrete = FALSE, alpha = 1, option = "D") +
  geom_point(data= nutsUTM.df, aes(Longitude.1, Latitude.1), size = 1.2, shape = 1, fill = "grey") +
  geom_raster(data=hill.df, aes(x=x,y=y,fill = alt),show.legend=FALSE) + scale_fill_gradientn(colors = grey(20:100/100)) +
  new_scale("fill") +
  geom_sf(data = mo.sf, fill=NA, lwd = 0.2) +
  coord_sf() +
	theme(legend.title=element_text(size=14),axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), 
axis.title.y=element_blank()) +  
  labs(color = "PCoA1") +
  theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5), plot.background=element_rect(fill='white'))

ggsave("../output/plots/pcoa1moorea0323.viridis.svg", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)
ggsave("../output/plots/pcoa1_moorea0323.viridis.png", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)

pcoa1_2 <- ggplot() + theme_void() +	geom_point(data=krig2, aes(x=x, y=y, colour=pred), size=4) + 
  scale_colour_gradientn(name="",colours = cubicl(10)) +
  geom_point(data= nutsUTM.df, aes(Longitude.1, Latitude.1), size = 1.2, shape = 1, fill = "grey") +
  geom_raster(data=hill.df, aes(x=x,y=y,fill = alt),show.legend=FALSE) + scale_fill_gradientn(colors = grey(20:100/100)) +
  new_scale("fill") +
  geom_sf(data = mo.sf, fill=NA, lwd = 0.2) +
  coord_sf() +
  theme(legend.title=element_text(size=14), axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.ticks=element_blank(), 
axis.title.x=element_blank(), axis.title.y=element_blank()) +  
  labs(color = "PCoA1") +
  theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5), plot.background=element_rect(fill='white'))

ggsave("../output/plots/pcoa1moorea0323.cubicl.svg", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)
ggsave("../output/plots/pcoa1_moorea0323.cubicl.png", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)

########## PCoA2

pcoa2_1 <- ggplot() + theme_void() +	geom_point(data=krig4, aes(x=x, y=y, colour=pred), size=4) + 
	scale_colour_gradientn(name="",colours = cubicl(10)) +
  geom_point(data= nutsUTM.df, aes(Longitude.1, Latitude.1), size = 1.2, shape = 1, fill = "grey") +
  geom_raster(data=hill.df, aes(x=x,y=y,fill = alt),show.legend=FALSE) + scale_fill_gradientn(colors = grey(20:100/100)) +
  new_scale("fill") +
  geom_sf(data = mo.sf, fill=NA, lwd = 0.2) +
  coord_sf() +
  theme(legend.title=element_text(size=14), axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), 
axis.title.y=element_blank()) +  
  labs(color = "PCoA2") +
  theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5), plot.background=element_rect(fill='white'))

ggsave("../output/plots/pcoa2_moorea0323.cubicl.svg", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)
ggsave("../output/plots/pcoa2_moorea0323.cubicl.png", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)

pcoa2_2 <- ggplot() + theme_void() +	geom_point(data=krig4, aes(x=x, y=y, colour=pred), size=4) + 
  scale_color_viridis(begin = 0, end = 1, direction = 1, discrete = FALSE, alpha = 1, option = "D") +
  geom_point(data= nutsUTM.df, aes(Longitude.1, Latitude.1), size = 1.2, shape = 1, fill = "grey") +
  geom_raster(data=hill.df, aes(x=x,y=y,fill = alt),show.legend=FALSE) + scale_fill_gradientn(colors = grey(20:100/100)) +
  new_scale("fill") +
  geom_sf(data = mo.sf, fill=NA, lwd = 0.2) +
  coord_sf() +
  theme(legend.title=element_text(size=14), axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), 
axis.title.y=element_blank()) +  
  labs(color = "PCoA2") +
  theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5), plot.background=element_rect(fill='white'))

ggsave("../output/plots/pcoa2_moorea0323.viridis.svg", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)
ggsave("../output/plots/pcoa2_moorea0323.viridis.svg.png", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)

########## Richness

richness_1 <- ggplot() + theme_void() +	geom_point(data=krig6, aes(x=x, y=y, colour=pred), size=4) + 
  scale_colour_gradientn(name="",colours = cubicl(10)) +
  geom_point(data= nutsUTM.df, aes(Longitude.1, Latitude.1), size = 1.2, shape = 1, fill = "grey") +
  geom_raster(data=hill.df, aes(x=x,y=y,fill = alt),show.legend=FALSE) + scale_fill_gradientn(colors = grey(20:100/100)) +
  new_scale("fill") +
  geom_sf(data = mo.sf, fill=NA, lwd = 0.2) +
  coord_sf() +
  theme(legend.title=element_text(size=14), axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), 
axis.title.y=element_blank()) +  
  labs(color = "Richness") +
  theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5), plot.background=element_rect(fill='white'))



ggsave("../output/plots/richness_moorea0323.cubicl.svg", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)
ggsave("../output/plots/richness_moorea0323.cubicl.png", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)

richness_2 <- ggplot() + theme_void() +	geom_point(data=krig6, aes(x=x, y=y, colour=pred), size=4) + 
  scale_color_viridis(begin = 0, end = 1, direction = 1, discrete = FALSE, alpha = 1, option = "D") +
  geom_point(data= nutsUTM.df, aes(Longitude.1, Latitude.1), size = 1.2, shape = 1, fill = "grey") +
  geom_raster(data=hill.df, aes(x=x,y=y,fill = alt),show.legend=FALSE) + scale_fill_gradientn(colors = grey(20:100/100)) +
  new_scale("fill") +
  geom_sf(data = mo.sf, fill=NA, lwd = 0.2) +
  coord_sf() +
  labs(fill = "Richness") + 
  theme(legend.title=element_text(size=14), axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), 
axis.title.y=element_blank()) +  
  labs(color = "Richness") + 
  theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5), plot.background=element_rect(fill='white'))

ggsave("../output/plots/richness_moorea0323.viridis.svg", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)
ggsave("../output/plots/richness_moorea0323.viridis.png", bg = "transparent", limitsize = FALSE, width = 10, height = 10, dpi = 220)

