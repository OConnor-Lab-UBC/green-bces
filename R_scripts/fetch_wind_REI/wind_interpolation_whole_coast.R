library(dplyr)
library(purrr)
library(tidyr)
library(sf)
library(terra)
library(raster)
library(fields)
library(tidyterra)

# read in shapefile 
wind_raw <- read_sf("./processed_evd/wind_data/whole_coast_rei/Copernicus_era5_summary_wide_whole_coast.shp") 

# add fetch sites to set extent 
whole_coast_fetch_points <- read_sf("./processed_evd/fetch/whole_coast/whole_coast_points_for_REI_calc.shp")
res = 30
raster_extent <- ext(whole_coast_fetch_points)
crs <- "EPSG:3347"

# create an empty raster 
empty_raster <-rast(
  extent = raster_extent,
  resolution = res,
  crs = crs,
  vals = 1
)

# extract xy from wind geometry
wind_xy <- st_coordinates(wind_raw$geometry)

# turn into a spatVector
wind_vect <- vect(wind_xy, crs = crs(wind_raw))

# crop to get rid of most wind points on land 
jb_water <- read_sf("./canada_coastline_polygons/geometries_used_to_crop_coastlines/whole_coast_geometry.shp")
jb_water <- vect(jb_water)
wind_cropped <- terra::crop(wind_vect, jb_water)

xy <- project(wind_cropped, crs) 
xy <- as.data.frame(geom(xy)) %>% 
  dplyr::select(x,y) 

# turn wind_raw into a vector so I can crop it and then turn it into a dataframe
wind_raw <- vect(wind_raw)
cropped <- terra::crop(wind_raw, jb_water)
wind_raw_df <- as.data.frame(cropped)

# first column is id
wind_raw_df <- wind_raw_df %>% 
  dplyr::select(-id) 

for (i in 1:ncol(wind_raw_df)) {
  
  print(i)
  v <- wind_raw_df[,i]
  
  tps <- Tps(xy, v)
  p <- rast(empty_raster)
  
  p <- interpolate(p, tps)
  
  # crop using land polygon 
  land_polygon <- vect("canada_coastline_polygons/coastlines/cropped_whole_coast.shp")
  land_polygon_projected <- land_polygon %>% 
    project(crs)
  
  cropped_rast <- mask(p, land_polygon_projected, inverse = TRUE)
  plot(cropped_rast)
  
  # write as raster and save
  writeRaster(cropped_rast, filename = 
                sprintf("./processed_evd/wind_data/whole_coast_rei/interp_wind_layers/%s_wind_raster_30mres_wholeCoast_2019_2026.tiff", 
                        colnames(wind_raw_df)[i]))
  
}

# claude's suggestion here 
# trying k-fold CV to see how well the tps is interpolating the relatively sparse data
set.seed(1)
folds <- sample(rep(1:5, length.out = nrow(xy)))

# try with only the 360 directionx
cv_rmse <- function(...) {
  sqrt(mean(unlist(lapply(1:5, function(k) {
    f <- fields::Tps(xy[folds != k, ], max_spd_360[folds != k], ...)
    (max_spd_360[folds == k] - predict(f, xy[folds == k, ]))^2
  }))))
}

cv_rmse() # lowest RMSE               
cv_rmse(lambda = 1e-3)
cv_rmse(df = 40)

# crop all of the files using the bathymetry layer 
bathy <- terra::rast("./processed_evd/environmental_layers/processed_bathy_layers/extended_30m_dem_whole_coast.tif")

# test crop
# mx_spd_360 <- terra::rast("./processed_evd/wind_data/whole_coast_rei/interp_wind_layers/mx_spd_90_wind_raster_30mres_wholeCoast_2019_2026.tiff")
# project to correct crs 
# nad83_bathy <- project(bathy, crs(mx_spd_360))
# crop_360_wind <- terra::crop(mx_spd_360)


