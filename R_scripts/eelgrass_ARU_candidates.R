library(readr)
library(dplyr)
library(tidyverse)
library(sf)
library(raster)
library(terra)

survey_stations <- read_csv("./CHCRP_wildlife/info_passive_sampling_stations_2025_11_03.csv")

## filter csv to ARUs that are not polar bear stations, have not been removed
stations_no_polarBear_only_deployed <- survey_stations %>% 
  filter(polar_bear_station == "no") %>% 
  filter(status == "deployed") 


write_csv(stations_no_polarBear_only_deployed, "./CHCRP_wildlife/stations_no_polarBear_only_deployed.csv")

# inspect eelgrass shp 
eelgrass_polys <- read_sf("./julian_shapefiles/eelgrass_beds/Eelgrass Beds.shp")

# filter to only include present eelgrass
present_eelgrass <- eelgrass_polys %>% 
  filter(Present == "Yes") # 60 polygons

# of present eelgrass, filter to only include healthy patches
healthy_eelgrass <- present_eelgrass %>% 
  filter(Heath == "Healthy") # 8 polygons

# write as shapefiles 
# first need to project the coordinates of the ARUs

prj4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
project_points <- st_crs(prj4string)

# make aru data into an sf object
aru_sf <- st_as_sf(stations_no_polarBear_only_deployed, coords = c("lon", "lat"), crs = project_points)
st_crs(aru_sf)

# project to same crs as julian's shapefiles
projected_aru_shapefile <- st_transform(aru_sf, crs = crs(eelgrass_polys))

st_write(obj = projected_aru_shapefile, 
         dsn = "./CHCRP_wildlife", 
         layer = "stations_no_polarBear_only_deployed", 
         driver = "ESRI Shapefile",
         delete_layer = TRUE)
# success 


# write eelgrass layers as shapefiles
st_write(obj = healthy_eelgrass, 
         dsn = "./julian_shapefiles/healthy_eelgrass//", 
         layer_options = "SHPT=POLYGONZ", # specify 3D structure
         driver = "ESRI Shapefile",
         delete_layer = TRUE)


st_write(obj = present_eelgrass, 
         dsn = "./julian_shapefiles/present_eelgrass//",
         layer_options = "SHPT=POLYGONZ",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)




