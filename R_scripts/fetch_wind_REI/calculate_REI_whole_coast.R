

#*Order of operations*
  
#1. Create a list of filenames for raster tiffs of i) wind frequency and ii) wind speed 
#2. Create raster stacks of wind frequency rasters and wind speed rasters
#3. Load in .rds of effective fetch data (df of dfs of effective fetch in 8 directions for each point)
#4. Ensure effective fetch is in the same CRS as wind data and project if necessary
#5. Extract frequency + speed from raster stacks at each fetch location
#6. Bind these into a matrix called wind_combine
#7. Calculate REI and place in df 
#8. Save REI as .csv file 


#######################################################################
##########Environmental Data Layers for SG sdm#########################
#######################################################################

#### Load Packages ####

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(sf)
library(raster)
library(sf)
library(prioritizr)
library(tidyterra)

# Read in interpolated wind data ----

# Create a list of file names (wind frequency)

wind_freq <- list.files(path = "./processed_evd/wind_data/whole_coast_rei/interp_wind_layers/",
                        pattern = "^freq.+tiff$",
                        full.names = T) 

wind_freq_stack <- rast(wind_freq)


# Raster stack of average daily maximum wind speed by direction

wind_spd <- list.files(path = "./processed_evd/wind_data/whole_coast_rei/interp_wind_layers/",
                       pattern = "^mx.+tiff$",
                       full.names = T) 

wind_spd_stack <- rast(wind_spd)


# Combine fetch data and wind data ----

# Read in effective fetch data
fetch_eff <- readRDS("./processed_evd/fetch/whole_coast/whole_coast_fetch_eff.rds") 

# this is an sf object, so need to convert geometries to SpatVector 
# while maintaining fetch tibble data to use in that form later 

# make geometry df with site names 
fetch_geometries <- fetch_eff %>% 
  select(site_names, geometry)

# make df of just fetch tibbles 
fetch_tibbles <- fetch_eff %>% 
  select(site_names, fetch_eff)

# turn geometries into a SpatVector 
fetch_geom_vect <- terra::vect(fetch_geometries) 

# check CRS alignment and project if necessary
crs(fetch_geom_vect)
crs(wind_freq_stack)

# project fetch_eff to crs of wind stacks 
wind_crs <- crs(wind_freq_stack)
fetch_geom_vect <- project(fetch_geom_vect, wind_crs) # project to wind CRS 


# will need to rename columns so later function can be applied 
freq_names <- c("site_names", "freq_135", "freq_180", "freq_225", "freq_270", "freq_315", "freq_360", "freq_45", "freq_90")
spd_names <- c("site_names", "mx_135", "mx_180", "mx_225", "mx_270", "mx_315", "mx_360", "mx_45", "mx_90")


# Extract frequency and speed values for fetch locations
extracted_freq <- terra::extract(wind_freq_stack, fetch_geom_vect, fun = NULL, ID = F, bind = TRUE) %>% 
  as.data.frame(.) 
  
# rename colnames
old_f_names <- colnames(extracted_freq)

extracted_freq <- extracted_freq |> 
  rename_at(old_f_names, ~freq_names)

# Extract speed values 
extracted_spd <- terra::extract(wind_spd_stack, fetch_geom_vect, fun = NULL, ID = F, bind = TRUE) %>% 
  as.data.frame(.) 

old_spd_names <- colnames(extracted_spd)

# rename cols
extracted_spd <- extracted_spd |>  
  rename_at(old_spd_names, ~spd_names)

# combine into a single df
wind_all <- left_join(extracted_freq, extracted_spd, by = "site_names")

wind_combine <- left_join(wind_all, fetch_tibbles, by = "site_names") %>% 
  drop_na() %>% 
  mutate(site_names = 1:nrow(.)) # rename site_names excluding na rows

# Calculate Relative Exposure Index ----
# caution: NA values in the wind dfs become zeros for relative exposure. 
# Best to remove these before proceeding (see previous line of code)

expo_by_site <- wind_combine %>%
  rowwise() %>% 
  mutate(freq_total = list(cbind(c_across(contains("freq")))), 
         max_wind = list(cbind(c_across(contains("mx")))), 
         wind_dir_freq = list(tibble(max_wind, 
                                     freq_total, 
                                     direction = c(135,180,225,270,315,360,45,90)))) %>%
  ungroup() %>% 
  mutate(., wind_dir_freq = map(wind_dir_freq, function(.data) {
    .data %>% 
      dplyr::arrange(., direction)
  }),
  fetch_eff = map(fetch_eff, function(.data){
    .data %>% 
      dplyr::arrange(., direction)
  })) %>% 
  mutate(REI = map2_dbl(wind_dir_freq, fetch_eff, # calculate relative exposure
                        ~sum((.x$max_wind*.x$freq_total*.y$fetch), na.rm = T))) %>%  
  dplyr::select(site_names, REI, geometry) # keep site, REI, and geometry 

# scale REI values on a scale of 0-1 
max(expo_by_site$REI) # 3731911

expo_by_site <- expo_by_site |> 
  mutate(scaled_REI = REI/3731911) |> 
  mutate(log_REI = log(REI))

# Save REI object

saveRDS(expo_by_site, './processed_evd/fetch/whole_coast/whole_coast_REI_10k_points.rds')

# Write to csv

st_write(obj = expo_by_site,
         dsn = './processed_evd/fetch/whole_coast/whole_coast_REI_10k_points.csv',
         layer_options = 'GEOMETRY=AS_XY',
         append = FALSE)


# write to shapefile
st_write(expo_by_site,
         dsn = "./processed_evd/fetch/whole_coast/whole_coast_REI_10k_points.shp",
         layer = "REI",
         driver = "ESRI Shapefile")


# Interpolating REI 
rei <- read_csv("./processed_evd/fetch/whole_coast/whole_coast_REI_10k_points.csv")

# read in shapefile with correct extent
whole_coast_fetch_points <- read_sf("./processed_evd/fetch/whole_coast/whole_coast_points_for_REI_calc.shp")
res = 30
raster_extent <- ext(whole_coast_fetch_points)
crs <- "EPSG:3347"

# create an empty raster 
temp_raster <-rast(
  extent = raster_extent,
  resolution = res,
  crs = crs,
  vals = 1
)

rei_xy <- rei %>% 
  select(X, Y)

rei_vect <- vect(rei_xy) # load in with correct CRS

xy <- as.data.frame(geom(rei_vect)) %>% select(x,y) # what space to interpolate over
v <- rei$scaled_REI# what value to interpolate

tps_rei <- Tps(xy, v) 
p <- rast(temp_raster)

p <- interpolate(p, tps_rei)
plot(p) 

cropped_temp_rast <- mask(p, land_polygon_projected, inverse = TRUE)
plot(cropped_temp_rast)


# write as raster and save
writeRaster(cropped_temp_rast, filename = "./processed_evd/environmental_layers/temperature_raster_20m_res_2019_leblanc.tiff")


