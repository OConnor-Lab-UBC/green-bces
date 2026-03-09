
library(dplyr)
library(tidyverse)
library(readr)
library(tidyr)
library(sf)
library(raster)
library(INLA)
library(terra)
library(ggmap)
library(stringi)
library(purrr)
library(units)


#' testing a small subsample of the coast NW(53.8491,-79.0629), SW(53.1502,-78.9418)
#' this way I will be able sample a high-resolution of points (100) in this small
#' area, calculate effective fetch, and see how different the results are from
#' each other 
#' 
#' using this method I will determine how high of resolution we need the fetch model
#' points to be.  


# Speed up raster processes
rasterOptions(chunksize = 1e+05, maxmemory = 1e+09)

# read in DEM as a raster (~ 450 m resolution)
DEM <- terra::rast("./env_data/GEBCO_JB_bathy/gebco_2025_n55.079_s50.9963_w-82.4486_e-78.1642.asc")

dem_fetch <- raster(DEM)
dem_fetch

# read in the cropped coast shapefile
cropped_coast <- st_read(dsn = "./processed_EvD/intersection_crop.shp")

# project dem to crs of cropped coastline
dem_extent_proj <- st_as_sfc(
  st_bbox(dem_fetch),
  crs = st_crs(dem_fetch)
) %>%
  st_transform(st_crs(cropped_coast))

# create buffer
dem_buffer_300km <- st_buffer(dem_extent_proj, 300000)
plot(dem_buffer_300km)

# crop to buffer extent
buffered_coast_jb <- cropped_coast %>%
  st_intersection(dem_buffer_300km) %>%
  st_cast("MULTIPOLYGON")

plot(buffered_coast_jb)


# Create 1-km land buffer
coast_1km_buffer <- st_buffer(cropped_coast, 1000)
plot(coast_1km_buffer)

#################
# try a 100km buffer 
dem_buffer_100km <- st_buffer(dem_extent_proj, 100000)

# crop to 100km buffer extent 
buffered_coast_jb <- cropped_coast %>%
  st_intersection(dem_buffer_100km) %>%
  st_cast("MULTIPOLYGON")


# create coarser land buffer
coast_100m_buffer <- st_buffer(cropped_coast, 100)
plot(coast_100m_buffer)

# let's see how many points this generates
# Remove points further than 100 m from coast
depth_sf_50buff <- depth_sf_50[st_intersects(depth_sf_50, coast_100m_buffer) %>% 
                                 lengths > 0,] # 358 points, much better

plot(depth_sf_50buff)


# Save pts as .rds files and Shapefiles
saveRDS(depth_sf_50buff, "./processed_EvD/fetch/REI_pts_cropped_coast.rds")

st_write(obj = depth_sf_50buff, 
         dsn = "./processed_EvD/fetch", 
         layer = "REI_pts_cropped_coast",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)

#' now these points are too close to land to meaningfully calculate fetch--
#' we just wanted them to be lower resolution 

#################

# Write to shapefile
#st_write(coast_1km_buffer,
         #dsn = "./env_data",
         #layer = "coast50k_buff_1km.shp",
         #driver = "ESRI Shapefile",
         #delete_layer = TRUE)


# write dem to geotiff
#writeRaster(dem_fetch, filename = "./env_data/dem_fetch.tif", overwrite = T)

# 3) Raster to points, filter by depth, apply 1-km land buffer

# Extract raster values to tibble
depth_sf_50 <- tibble(Lon = coordinates(dem_fetch)[,1], 
                      Lat = coordinates(dem_fetch)[,2],
                      depth = raster::values(dem_fetch)) %>% 
  
  #Confine points to 50 m or shallower (bathy points are negative, so should be between 0 and -50)
  filter(., depth >= -50) %>%  
  filter(., depth <= 0) %>% # 453361 pts
  st_as_sf(., coords = c("Lon","Lat"), crs = 4326) # convert to sf object


depth_sf_50
coast_1km_buffer

# project depth_sf_50 crs
depth_sf_50 <- st_transform(depth_sf_50, st_crs(cropped_coast))

# Remove points further than 1 km from coast
depth_sf_50buff <- depth_sf_50[st_intersects(depth_sf_50, coast_1km_buffer) %>% 
                                 lengths > 0,] # 4761 points

plot(depth_sf_50buff)
# this is the one we want to work with

# check if any points are on land and discard

on_land_50 <- st_intersects(depth_sf_50buff, cropped_coast) %>%
  lengths(.) == 0

depth_sf_50buff <- depth_sf_50buff[on_land_50,] # 4266 points



# Save pts as .rds files and Shapefiles
saveRDS(depth_sf_50buff, "./processed_EvD/fetch/REI_subsample_pts.rds")

st_write(obj = depth_sf_50buff, 
         dsn = "./processed_EvD/fetch", 
         layer = "REI_subsample_pts",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)


# create a random subset of depth_sf_50buff so that it is runnable on the fetch model
# the current resolution is ~300m 
# let's see how many points are restricted to depths of 10m or less

depth_sf_10buff <- depth_sf_50buff %>% 
  filter(depth >= -10) # limited to 3008 points 
plot(depth_sf_10buff)

# now we can take a random sample of these--let's keep 1000 points
random_1000 <- depth_sf_10buff[sample(nrow(depth_sf_10buff), 1000), ]

# random_1000 will now be the site_layer input
# make variables to run fetch model  

# polygon_layer 
polygon_layer <- read_sf("./processed_EvD/intersection_crop.shp")
polygon_layer <- st_as_sf(polygon_layer)

# site_layer
site_layer <- st_as_sf(random_1000)
site_layer <- site_layer %>% 
  dplyr::select(geometry)

# project to polyon crs
site_layer <- st_transform(site_layer, st_crs(polygon_layer))

site_layer$site_var <- NA
site_var <- site_layer$site_var

# site_names
site_names <- c(1:1000) # length of site_layer
site_names <- as.character(site_names)


#one final check to make sure every variable looks good (site and polygon layers should have same CRS)
polygon_layer
site_layer
site_names

subsample_fetch <- fetch_parallel(polygon_layer, site_layer, max_dist = 300, n_directions = 8, site_var = 
                 site_var)

## effective fetch
# Source functions ----

source('./processed_EvD/fetch/obrienjm25-REI-WaveExp-fbf17f9/Code/roll_recycle_fun.R')

# Controls ----

# relative contribution of fetch vectors surrounding each heading to calculation of effective fetch

weights <- cos(c(45,33.75,22.5,11.25,0,11.25,22.5,33.75,45) * (pi/180))
# Coordinate reference system of fetch calculations (integer EPSG code)

crs_fetch <- 3347 # NAD83 / Statistics Canada Lambert


# site names column was full of NAs 
site_names_vec <- rep(c(1:1000), each = 32)
subsample_fetch$site_names <- site_names_vec
subsample_fetch$site_names <- factor(subsample_fetch$site_names)

fetch_proxies <- subsample_fetch %>% 
  group_by(geometry) %>% 
  mutate(min_fetch = min(fetch_length), sum_fetch = sum(fetch_length)) %>% 
  group_by(min_fetch, sum_fetch) %>% 
  summarise()

# calculating effective fetch

fetch_eff <- subsample_fetch %>%
  group_by(site_names) %>%
  summarise(
    fetch_eff = list({
      f <- roll.recycle(fetch_length, 9, 8, by = 4)
      f <- as.vector((weights %*% f) / sum(weights))
      tibble(
        direction = c(360, seq(45, 315, 45)),
        fetch = f
      )
    }),
    .groups = "drop"
  ) %>%
  tidyr::unnest(fetch_eff)


saveRDS(fetch_eff, "./processed_EvD/fetch/fetch_eff.rds")

st_write(obj = fetch_eff, 
         dsn = "./processed_EvD/fetch", 
         layer = "fetch_eff",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)

#' points at 300m resolution are different enough from each other along the coast that this resolution might
#' actually be too coarse


geometries <- fetch_eff_wgs84$geometry

# doing some more visualization
fetch_eff_XY <- geometries %>% 
  st_coordinates() %>% # retrieves X, Y coordinates in a matrix
  as.data.frame()

fetch_eff_comb <- cbind(fetch_eff_wgs84, fetch_eff_XY)

ggplot(fetch_eff_comb, aes(X, Y, fill= fetch)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  theme_minimal()

fetch_eff_wgs84 <- st_transform(fetch_eff, crs = 4326)



####### Calculating fetch just for the west coast (sf = ec_diff)


# read in the cropped coast shapefile
cropped_coast <- st_read(dsn = "./processed_EvD/shapefiles/ec_diff.shp")
plot(cropped_coast)

# project dem to crs of cropped coastline
dem_extent_proj <- st_as_sfc(
  st_bbox(dem_fetch),
  crs = st_crs(dem_fetch)
) %>%
  st_transform(st_crs(cropped_coast))

# create buffer
dem_buffer_300km <- st_buffer(dem_extent_proj, 300000)

# crop to buffer extent
buffered_coast_jb <- cropped_coast %>%
  st_intersection(dem_buffer_300km) %>%
  st_cast("MULTIPOLYGON")

plot(buffered_coast_jb)


# Create 1-km land buffer
coast_1km_buffer <- st_buffer(cropped_coast, 1000)
plot(coast_1km_buffer)

#################
# try a 100km buffer 
dem_buffer_100km <- st_buffer(dem_extent_proj, 100000)

# crop to 100km buffer extent 
buffered_coast_jb <- cropped_coast %>%
  st_intersection(dem_buffer_100km) %>%
  st_cast("MULTIPOLYGON")


# create coarser land buffer
coast_100m_buffer <- st_buffer(cropped_coast, 100)
plot(coast_100m_buffer)

#################

# Write to shapefile
#st_write(coast_1km_buffer,
#dsn = "./env_data",
#layer = "coast50k_buff_1km.shp",
#driver = "ESRI Shapefile",
#delete_layer = TRUE)


# write dem to geotiff
#writeRaster(dem_fetch, filename = "./env_data/dem_fetch.tif", overwrite = T)

# 3) Raster to points, filter by depth, apply 1-km land buffer

# Extract raster values to tibble
depth_sf_50 <- tibble(Lon = coordinates(dem_fetch)[,1], 
                      Lat = coordinates(dem_fetch)[,2],
                      depth = raster::values(dem_fetch)) %>% 
  
  #Confine points to 50 m or shallower (bathy points are negative, so should be between 0 and -50)
  filter(., depth >= -50) %>%  
  filter(., depth <= 0) %>% # 453361 pts
  st_as_sf(., coords = c("Lon","Lat"), crs = 4326) # convert to sf object


depth_sf_50
coast_1km_buffer

# Save pts as .rds files and Shapefiles
saveRDS(depth_sf_50buff, "./processed_EvD/fetch/REI_pts_cropped_coast.rds")

st_write(obj = depth_sf_50buff, 
         dsn = "./processed_EvD/fetch", 
         layer = "REI_pts_cropped_coast",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)

#' now these points are too close to land to meaningfully calculate fetch--
#' we just wanted them to be lower resolution 


# project depth_sf_50 crs
depth_sf_50 <- st_transform(depth_sf_50, st_crs(cropped_coast))

# Remove points further than 1 km from coast
depth_sf_50buff <- depth_sf_50[st_intersects(depth_sf_50, coast_1km_buffer) %>% 
                                 lengths > 0,] # 8346 points

plot(depth_sf_50buff)
# this is the one we want to work with

# check if any points are on land and discard

on_land_50 <- st_intersects(depth_sf_50buff, cropped_coast) %>%
  lengths(.) == 0

depth_sf_50buff <- depth_sf_50buff[on_land_50,] # 5423 points



# Save pts as .rds files and Shapefiles
saveRDS(depth_sf_50buff, "./processed_EvD/fetch/REI_westcoast_pts.rds")

st_write(obj = depth_sf_50buff, 
         dsn = "./processed_EvD/fetch", 
         layer = "REI_westcoast_pts",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)


# create a random subset of depth_sf_50buff so that it is runnable on the fetch model
# the current resolution is ~300m 
# let's see how many points are restricted to depths of 10m or less

depth_sf_10buff <- depth_sf_50buff %>% 
  filter(depth >= -10) # limited to 3891 points 
plot(depth_sf_10buff)

# depth_sf_10buff will now be the site_layer input
# make variables to run fetch model  

# polygon_layer 
polygon_layer <- read_sf("./processed_EvD/shapefiles/ec_diff.shp")
polygon_layer <- st_as_sf(polygon_layer)

# site_layer
site_layer <- st_as_sf(depth_sf_10buff)
site_layer <- site_layer %>% 
  dplyr::select(geometry)

# project to polyon crs
site_layer <- st_transform(site_layer, st_crs(polygon_layer))

site_layer$site_var <- NA
site_var <- site_layer$site_var

# site_names
site_names <- c(1:3891) # length of site_layer
site_names <- as.character(site_names)


#one final check to make sure every variable looks good (site and polygon layers should have same CRS)
polygon_layer
site_layer
site_names

westcoast_fetch <- fetch_parallel(polygon_layer, site_layer, max_dist = 300, n_directions = 8, site_var = 
                                    site_var)

## effective fetch
# Source functions ----

source('./processed_EvD/fetch/obrienjm25-REI-WaveExp-fbf17f9/Code/roll_recycle_fun.R')

# Controls ----

# relative contribution of fetch vectors surrounding each heading to calculation of effective fetch

weights <- cos(c(45,33.75,22.5,11.25,0,11.25,22.5,33.75,45) * (pi/180))
# Coordinate reference system of fetch calculations (integer EPSG code)

crs_fetch <- 3347 # NAD83 / Statistics Canada Lambert


# site names column was full of NAs 
site_names_vec <- rep(c(1:3891), each = 32)
westcoast_fetch$site_names <- site_names_vec
westcoast_fetch$site_names <- factor(westcoast_fetch$site_names)

fetch_proxies_west <- westcoast_fetch %>% 
  group_by(geometry) %>% 
  mutate(min_fetch = min(fetch_length), sum_fetch = sum(fetch_length)) %>% 
  group_by(min_fetch, sum_fetch) %>% 
  summarise()

# calculating effective fetch

westcoast_fetch_eff <- westcoast_fetch %>%
  group_by(site_names) %>%
  summarise(
    fetch_eff = list({
      f <- roll.recycle(fetch_length, 9, 8, by = 4)
      f <- as.vector((weights %*% f) / sum(weights))
      tibble(
        direction = c(360, seq(45, 315, 45)),
        fetch = f
      )
    }),
    .groups = "drop"
  ) %>%
  tidyr::unnest(fetch_eff)


saveRDS(westcoast_fetch_eff, "./processed_EvD/fetch/westcoast_fetch_eff.rds")

st_write(obj = westcoast_fetch_eff, 
         dsn = "./processed_EvD/fetch", 
         layer = "westcoast_fetch_eff",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)


#### Calculating fetch for the east coast (sf = east_coast_out)
  

# read in the cropped coast shapefile
cropped_eastcoast <- st_read(dsn = "./processed_EvD/east_coast_out.shp")
plot(cropped_eastcoast)

# project dem to crs of cropped coastline
dem_extent_proj <- st_as_sfc(
  st_bbox(dem_fetch),
  crs = st_crs(dem_fetch)
) %>%
  st_transform(st_crs(cropped_eastcoast))

# create buffer
dem_buffer_300km <- st_buffer(dem_extent_proj, 300000)

# crop to buffer extent
buffered_eastcoast_jb <- cropped_eastcoast %>%
  st_intersection(dem_buffer_300km) %>%
  st_cast("MULTIPOLYGON")

plot(buffered_eastcoast_jb)

# Create 1-km land buffer
coast_1km_buffer <- st_buffer(cropped_eastcoast, 1000)
plot(coast_1km_buffer)


# Write to shapefile
#st_write(coast_1km_buffer,
#dsn = "./env_data",
#layer = "coast50k_buff_1km.shp",
#driver = "ESRI Shapefile",
#delete_layer = TRUE)


# write dem to geotiff
#writeRaster(dem_fetch, filename = "./env_data/dem_fetch.tif", overwrite = T)

# 3) Raster to points, filter by depth, apply 1-km land buffer

# Extract raster values to tibble
depth_sf_50 <- tibble(Lon = coordinates(dem_fetch)[,1], 
                      Lat = coordinates(dem_fetch)[,2],
                      depth = raster::values(dem_fetch)) %>% 
  
  #Confine points to 50 m or shallower (bathy points are negative, so should be between 0 and -50)
  filter(., depth >= -50) %>%  
  filter(., depth <= 0) %>% # 453361 pts
  st_as_sf(., coords = c("Lon","Lat"), crs = 4326) # convert to sf object


depth_sf_50
coast_1km_buffer

# project depth_sf_50 crs
depth_sf_50 <- st_transform(depth_sf_50, st_crs(cropped_eastcoast))

# Remove points further than 1 km from coast
depth_sf_50buff <- depth_sf_50[st_intersects(depth_sf_50, coast_1km_buffer) %>% 
                                 lengths > 0,] # 23917 points

plot(depth_sf_50buff)
# this is the one we want to work with


# Save pts as .rds files and Shapefiles
saveRDS(depth_sf_50buff, "./processed_EvD/fetch/REI_pts_cropped_coast.rds")

st_write(obj = depth_sf_50buff, 
         dsn = "./processed_EvD/fetch", 
         layer = "REI_pts_cropped_coast",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)


# check if any points are on land and discard

on_land_50 <- st_intersects(depth_sf_50buff, cropped_eastcoast) %>%
  lengths(.) == 0

depth_sf_50buff <- depth_sf_50buff[on_land_50,] # 18812 points



# Save pts as .rds files and Shapefiles
saveRDS(depth_sf_50buff, "./processed_EvD/fetch/REI_eastcoast_pts.rds")

st_write(obj = depth_sf_50buff, 
         dsn = "./processed_EvD/fetch", 
         layer = "REI_eastcoast_pts",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)


# create a random subset of depth_sf_50buff so that it is runnable on the fetch model
# the current resolution is ~300m 
# let's see how many points are restricted to depths of 10m or less

depth_sf_10buff <- depth_sf_50buff %>% 
  filter(depth >= -10) # limited to 13540 points 
plot(depth_sf_10buff)

# depth_sf_10buff will now be the site_layer input
# make variables to run fetch model  

# polygon_layer 
polygon_layer <- read_sf("./processed_EvD/east_coast_out.shp")
polygon_layer <- st_as_sf(polygon_layer)

# site_layer
site_layer <- st_as_sf(depth_sf_10buff)
site_layer <- site_layer %>% 
  dplyr::select(geometry)

# project to polyon crs
site_layer <- st_transform(site_layer, st_crs(polygon_layer))

site_layer$site_var <- NA
site_var <- site_layer$site_var

# site_names
site_names <- c(1:13540) # length of site_layer
site_names <- as.character(site_names)


#one final check to make sure every variable looks good (site and polygon layers should have same CRS)
polygon_layer
site_layer
site_names

eastcoast_fetch <- fetch_parallel(polygon_layer, site_layer, max_dist = 300, n_directions = 8, site_var = 
                                    site_var)

## effective fetch
# Source functions ----

source('./processed_EvD/fetch/obrienjm25-REI-WaveExp-fbf17f9/Code/roll_recycle_fun.R')

# Controls ----

# relative contribution of fetch vectors surrounding each heading to calculation of effective fetch

weights <- cos(c(45,33.75,22.5,11.25,0,11.25,22.5,33.75,45) * (pi/180))
# Coordinate reference system of fetch calculations (integer EPSG code)

crs_fetch <- 3347 # NAD83 / Statistics Canada Lambert


# site names column was full of NAs 
site_names_vec <- rep(c(1:13540), each = 32)
eastcoast_fetch$site_names <- site_names_vec
eastcoast_fetch$site_names <- factor(eastcoast_fetch$site_names)

fetch_proxies_east <- eastcoast_fetch %>% 
  group_by(geometry) %>% 
  mutate(min_fetch = min(fetch_length), sum_fetch = sum(fetch_length)) %>% 
  group_by(min_fetch, sum_fetch) %>% 
  summarise()

# calculating effective fetch

eastcoast_fetch_eff <- eastcoast_fetch %>%
  group_by(site_names) %>%
  summarise(
    fetch_eff = list({
      f <- roll.recycle(fetch_length, 9, 8, by = 4)
      f <- as.vector((weights %*% f) / sum(weights))
      tibble(
        direction = c(360, seq(45, 315, 45)),
        fetch = f
      )
    }),
    .groups = "drop"
  ) %>%
  tidyr::unnest(fetch_eff)


saveRDS(eastcoast_fetch_eff, "./processed_EvD/fetch/eastcoast_fetch_eff.rds")

st_write(obj = eastcoast_fetch_eff, 
         dsn = "./processed_EvD/fetch", 
         layer = "eastcoast_fetch_eff",
         driver = "ESRI Shapefile",
         delete_layer = TRUE)





