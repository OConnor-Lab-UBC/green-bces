
#######################################################################
####### Download and summarize surface wind components from CDS #######
#######################################################################

# Description: ----

# This script makes data requests to Copernicus' Climate Data Stores (CDS) for
# surface wind data from the ERA-5 reanalysis and summarizes maximum wind speed
# and frequency by direction over the specified time period

# Requirements: ----

# Free personal account to CDS with username and key:

# https://cds.climate.copernicus.eu/user/register

# Note: ----

# large data requests over a larger bounding box and multiple years will
# need to be split into multiple queries

# Load packages ----

library(dplyr)
library(purrr)
library(tidyr)
library(data.table)
library(stringr)
library(sf)
library(ecmwfr)
library(keyring)
library(doParallel)
library(foreach)
library(ncdf4)
# library(ncdf4.helpers)

# Housekeeping ----

# Sys.setenv(tz = 'GMT') # if using Azure VM, will need to set the system time zone

data_dir <- "./processed_evd/wind_data/" # directory within project folder to store downloaded data

my_uid <- "mayagreen6644@gmail.com" # your CDS store user ID

my_key <- "b89f66bf-6001-4ccb-970a-128d9da8d7bb" # your CDS store key

crs <- 32617 # WGS 84 / UTM zone 17N # projected coordinate reference system  (integer EPSG code)
# crs <- 26921 # (UTM21)

# Request and download wind data from CDS ----

# Add login details to your local keyring

wf_set_key(user = my_uid, # personal User ID for CDS
           key = my_key # key provided by CDS
          )

# Provide query arguments for CDS data request

years <- str_pad(2019:2020, 1) # year range
months <- str_pad(1:12, 2, "left", 0) # Jan - Dec
days <- str_pad(1:31,2,"left","0") # 31 days
hours <- str_c(0:23,"00",sep=":")%>%str_pad(5,"left","0") # 24 hrs
bbox <- " 54.64/-80.21/51.37/-78.25" # bounding box for region of interest (N, W, S, E)

# Loop through years, provide query as list, and download data

n_cores <- detectCores() - 1

cl <- makePSOCKcluster(n_cores)
clusterExport(cl, "wf_request")
registerDoParallel(cl)

foreach (i =  1:length(years)) %dopar% {
  
  request <- list(
    product_type = "reanalysis",
    data_format = "netcdf",
    variable = c("10m_u_component_of_wind", "10m_v_component_of_wind"),
    year = years[i],
    month = months,
    day = days,
    time = hours,
    area = bbox,
    dataset_short_name = "reanalysis-era5-single-levels",
    target = paste0("era5_", years[i], ".nc")
  )
  
  wf_request(user = my_uid,
             request = request,   
             transfer = TRUE,  
             path = data_dir,
             verbose = FALSE)
}

stopCluster(cl)

# Summarize wind data ----

era5 <- list.files(path = data_dir, pattern = "era5_2019", full.names = T) %>% 
  map(., nc_open)

time <- map(era5, ~ncvar_get(.x, "valid_time"))

u10 <- map(era5, ~as.vector(ncvar_get(.x, 'u10')))
v10 <- map(era5, ~as.vector(ncvar_get(.x, 'v10')))

latitude <- map(era5, ~as.vector(ncvar_get(.x, 'latitude'))) 

longitude <- map(era5, ~as.vector(ncvar_get(.x, 'longitude')))

latitude_long <-  map2(latitude, longitude, ~lapply(.x, function(x) rep(x, length(.y)))) %>% 
  map(., unlist) %>% 
  map2(., time, ~rep(.x, length(.y)))

longitude_long <-  map2(longitude, latitude, ~rep(.x, length(.y))) %>% 
  map2(., time, ~rep(.x, length(.y)))

# convert to arc seconds? 
date_time <- pmap(list(time, latitude, longitude), ~lapply(..1, function(x) rep(x, (length(..2) * length(..3))))) %>%  
  map(., unlist) %>% 
  map(., ~as.POSIXct(.x * 3600, origin = "1900-01-01 00:00:00", tz = "UTC"))

wind_data <- pmap(list(longitude_long, latitude_long, date_time, u10, v10), 
                  ~data.frame(longitude = ..1, latitude = ..2, date_time = ..3, u10 =  ..4, v10 = ..5)) %>% 
  rbindlist() %>% 
  mutate(wind_spd = sqrt((u10^2) + (v10^2)),
         wind_dir = (180 + ((180/pi)*atan2(u10, v10))) %% 360) %>% 
  # counts the number of entries in 45 degree bins around 8 compass headings (N, NE, E, etc.)
  # In other words, this is the frequency that the wind blew from each direction of interest      
  mutate(direction = cut(wind_dir, 
                         breaks = c(0,22.5,67.5,112.5,157.5,202.5,247.5,292.5,337.5,361),
                         right = TRUE,
                         labels = c(360,45,90,135,180,225,270,315,360))) %>%  
  # Nest data frame by lon and lat
  group_by(longitude, latitude) %>% 
  nest() %>% 
  # calculates frequency from each direction
  mutate(., q95 = map(data, function(.data){
    .data %>% 
      filter(wind_spd > quantile(wind_spd, probs = 0.95))
  })) %>% 
  mutate(wind_dir_cbin = map(q95, function(.data) {
    .data %>%
      dplyr::select(direction) %>% 
      table() %>% 
      data.frame() %>%
      mutate(freq_total = Freq/sum(Freq), 
             direction = as.numeric(as.character(direction)))})) %>%
  # average speed by direction
  mutate(., avg_wind = map(q95, function(.data) {
    .data %>% 
      mutate(date = as.Date(date_time)) %>% 
      separate(date, into = c('year', 'month', 'day')) %>% 
      group_by(year, month, day, direction) %>% 
      summarise(daily_max = max(wind_spd)) %>% 
      group_by(year, month, direction) %>% 
      summarise(monthly_max = mean(daily_max)) %>% 
      group_by(direction, .drop = FALSE) %>% 
      summarise(max_wind = mean(monthly_max)) %>% 
      mutate_if(is.factor, as.character) %>%
      mutate_if(is.character, as.numeric) %>%
      mutate(max_wind = na_if(x = max_wind, y = NA),
             max_wind = nafill(x = max_wind, fill = 0)) %>% 
      data.frame()})) %>% 
  # Join average wind speed and frequency into one table per station
  transmute(wind_dir_freq = map2(avg_wind, wind_dir_cbin, left_join)) %>%  
  st_as_sf(., coords = c('longitude','latitude'), crs = 4326) %>% # convert to sf spatial object
  st_transform(crs) # project to WGS 84 UTM 17

# save to RDS object

saveRDS(wind_data, 'processed_evd/wind_data/copernicus_era5_with_arctus.rds')

wind_data <- wind_data %>% 
  mutate(id = rownames(wind_data))

# cast to wide format and write to shapefile

wind_unnest <- unnest(wind_data, cols = wind_dir_freq)

wind_wide <- dplyr::select(wind_unnest, -Freq) %>% 
  st_drop_geometry() %>% 
  rename(., mx_spd = max_wind, freq = freq_total) %>% 
  pivot_wider(data = .,
              id_cols = id,
              names_from = direction,
              values_from = c(mx_spd, freq)) %>% 
  left_join(., dplyr::select(wind_data, id)) %>% 
  st_set_geometry(., value = .$geometry)

st_write(wind_wide,
         dsn = './processed_evd/wind_data/',
         layer = "Copernicus_era5_summary_wide",
         driver = 'ESRI Shapefile')

