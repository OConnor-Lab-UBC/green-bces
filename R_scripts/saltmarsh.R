library(dplyr)
library(tidyr)
library(readr)
library(sf)


shoreline_type <- read_sf("shapefiles/coastal_shore_type.shp")
saltmarsh <- shoreline_type %>% 
  filter(HABITAT_en == "Salt marsh")

saltmarsh_area <- sum(saltmarsh$Area/1000000)
# there are 417 km2 of saltmarsh as indicated by the map  


fresh_marsh <- shoreline_type %>% 
  filter(HABITAT_en == "Freshwater marsh")

fresh_marsh_area <- sum(fresh_marsh$Area/1000000)
# 162 km2 of freshwater marsh indicated by the map
# 579 total km2 of marsh on east coast