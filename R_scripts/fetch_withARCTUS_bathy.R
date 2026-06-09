library(dplyr)
library(tidyverse)
library(readr)
library(tidyr)
library(sf)
library(raster)
library(terra)
library(ggmap)
library(stringi)
library(purrr)
library(units)
library(tidyterra)
install.packages("devtools")
devtools::install_github("anguswg-ucsb/fetchr")
library(fetchr)

# recreating a fetch layer using high resolution (20m) bathymetry layer from ARCTUS 
#' workflow: 
#' load in bathymetry rasters
#' load in land polygons
#' load in polygon used to clip coast to extent 
#' generate points over RasterLayer of Chisasibi bathymetry
#' plot to see if looks correct 

# Speed up raster processes
rasterOptions(chunksize = 1e+05, maxmemory = 1e+09)

# start with Chisasibi DEM

# read in DEM as a raster (~ 20 m resolution)
dem_fetch_chi <- terra::rast("./JB_coastline_layers/Chissassibi_SDB_Sentinel2.tif") 

# make into a raster with raster package specifications
# dem_fetch_chi_raster <- raster(dem_fetch_chi) # formal class RasterLayer

# read in the east coast shapefile
east_coast_land <- st_read(dsn = "./shapefiles/clipped_east_coast_2021_bdry.shp")
plot(east_coast_land$geometry)

# turn into a spatvector
east_coast_polygon <- vect(east_coast_land)

# project dem to crs of coastline
dem_extent_proj <- st_as_sfc(
  st_bbox(dem_fetch_chi),
  crs = st_crs(dem_fetch_chi)
) %>%
  st_transform(st_crs(east_coast_land))

# generate points over Chisasibi raster
raster_points <- as.points(dem_fetch_chi, values = TRUE)

# filter points by value < 0
water_points <- raster_points %>% 
  tidyterra::filter(SDB_depth <= 0)

# crop water_points using land mask
# set crs
east_coast_polygon <- project(east_coast_polygon, water_points)

# run fetch function on bathymetry points for Chisasibi 
# setwd("./processed_evd/obrienjm25-REI-WaveExp-fbf17f9/Code/")
# source("./fetch_parallel.R")


polygon_layer <- st_as_sf(east_coast_polygon)
site_layer <- st_as_sf(water_points)

# there are some sites on land, so we need to get rid of those
on_land <- st_intersects(site_layer, polygon_layer) %>% 
  lengths(.) == 0

site_layer <- site_layer[on_land,]

site_layer_sample1000 <- site_layer %>% 
  slice_sample(n = 1000) %>% 
  mutate(sites = as.character(1:nrow(site_layer_sample1000))) %>% 
  select(-SDB_depth, -Standard_deviation)

source("./R_scripts/fetch_functions.R")
fetch_first_try <- fetch_parallel(polygon_layer = polygon_layer, site_var = "sites", site_layer = site_layer_sample1000, 
                                  max_dist = 300, n_directions = 8)


saveRDS(fetch_first_try, "./processed_evd/fetch/fetch_points_arctus.rds")

st_write(obj = fetch_first_try, 
         dsn = "./processed_evd/fetch", 
         layer = "fetch_points_arctus",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)

# TEST VERSION WITH CHANGED CODE #
    # pls remove me maya :D #
#fetch_first_try <- fetch_TEST(polygon_layer = polygon_layer, site_names = site_names, site_layer = first_100, 
                                  #max_dist = 300, n_directions = 8)
#fetch_first_try <- fetch_parallel_TEST(polygon_layer = polygon_layer, site_layer = first_100, 
#                              max_dist = 300, site_var = "site_names", n_directions = 8)

