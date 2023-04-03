library(tidyverse, dplyr)
library(sp)
library(maps)
library(maptools)

setwd("~/Documents/GitHub/SARE")

pollen_inventory <- read.csv("pollen_inventory.csv")

latlong2county <- function(pollen_inventory) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per county
  counties <- map('county', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(counties$names, ":"), function(x) x[1])
  counties_sp <- map2SpatialPolygons(counties, IDs=IDs,
                                     proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pollen_inventory to a SpatialPoints object 
  pointsSP <- SpatialPoints(pollen_inventory, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, counties_sp)
  
  # Return the county names of the Polygons object containing each point
  countyNames <- sapply(counties_sp@polygons, function(x) x@ID)
  countyNames[indices]
}

counties <- data.frame(x = pollen_inventory$longitude, y = pollen_inventory$latitude)


pollen_inventory %>% 
  mutate(county = latlong2county(counties)) ->
  
  pollen_county
  
  
  
  

