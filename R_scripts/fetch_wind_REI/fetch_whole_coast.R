library(dplyr)
library(tidyverse)
library(readr)
library(tidyr)
library(sf)
library(raster)
library(terra)
library(stringi)
library(purrr)
library(units)
library(viridis)

# Speed up raster processes
rasterOptions(chunksize = 1e+05, maxmemory = 1e+09)

# read in DEM as a raster (30 m resolution)
DEM <- terra::rast("./processed_evd/environmental_layers/processed_bathy_layers/extended_30m_dem_whole_coast.tif")

# read in the cropped coast shapefile
cropped_whole_coast <- st_read(dsn = "./canada_coastline_polygons/coastlines/cropped_whole_coast.shp")
cropped_whole_coast <- vect(cropped_whole_coast)
cropped_whole_coast <- terra::aggregate(cropped_whole_coast)

# project dem to crs of coastline 
dem_extent_proj <- st_as_sfc(
  st_bbox(DEM, crs = st_crs(DEM)) %>%
  st_transform(st_crs(cropped_whole_coast)))

# generate points over the dem and sample 10000
depth_points <- terra::as.points(DEM)
points10k <- sample(depth_points, 10000)

# project points to land crs
points10k <- project(points10k, crs(cropped_whole_coast))

# turn points, cropped coast into sf objects
points10k_sf <- st_as_sf(points10k)
cropped_whole_coast_sf <- st_as_sf(cropped_whole_coast)

# check if any points are on land and discard
on_land_points <- st_intersects(points10k_sf, cropped_whole_coast_sf) %>%
  lengths(.) == 0

water_points <- points10k_sf[on_land_points,] 

# add sites column to depth_points and rename depth col
water_points <- water_points |> 
  mutate(sites = 1:nrow(water_points)) |> 
  dplyr::select(-Standard_deviation)

# save water_points as shapefile and move to fetch calculation
saveRDS(water_points, "./processed_evd/fetch/whole_coast/whole_coast_points_for_REI_calc.rds")

st_write(obj = water_points, 
         dsn = "./processed_evd/fetch/whole_coast",
         delete_layer = TRUE,
         layer = "whole_coast_points_for_REI_calc",
         driver = "ESRI Shapefile")

# rename land and point vectors to put into fetch parallel function
polygon_layer <- cropped_whole_coast_sf
site_layer <- water_points


# calculate fetch at each point
source("./R_scripts/fetch_wind_REI/fetch_functions.R") # source function
fetch_whole_coast <- fetch_parallel(polygon_layer = polygon_layer, site_var = "sites", site_layer = water_points, 
                                 max_dist = 300, n_directions = 8)

write_csv(fetch_whole_coast, "./processed_evd/fetch/whole_coast/raw_fetch_whole_coast.csv")





