#######################################################################
########## Calculate effective fetch from raw fetch vectors ###########
#######################################################################

# Description: ----

# This script takes fetch distances in 32 directions around grid points,
# the output of Python fetch script, and calculates for each point:
#
# 1) Minimum fetch: the shortest fetch distance (proxy for distance to coast)
# 2) Sum fetch: the sum of all 32 fetch distances (proxy for wave exposure)
# 3) Effective fetch: in 8 directions (N, NE, E, SE, S, SW, W, NW); 
#    needed to calculate Relative Exposure Index (REI)*
#
# *For description of effective fetch and REI, see:
#
# Fonseca MS, Bell SS (1998) Influence of physical setting on seagrass
#   landscapes near Beaufort, North Carolina, USA. Mar Ecol Prog Ser
#   171:109-121

# Requirements: ----

# csv file(s) containing outputs of Python fetch script
# i.e. fetch distance (metres) in 32 directions around grid points
# One row per point

# Load packages ----

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(data.table)
library(sf)
library(vroom)

# Source functions ----

source("./R_scripts/fetch_wind_REI/obrienjm25-REI-WaveExp-fbf17f9/Code/roll_recycle_fun.R")

# Controls ----

# relative contribution of fetch vectors surrounding each heading to calculation of effective fetch

weights <- cos(c(45,33.75,22.5,11.25,0,11.25,22.5,33.75,45) *(pi/180))

# Path to project directory (string)

projDir <- "./processed_evd/fetch"

# Coordinate reference system of fetch calculations (integer EPSG code)

crs_fetch <- 32617 # WGS 84 / UTM zone 17N

# Combine fetch data from all points into one data table ----

# list of fetch files
fetch_csv <- list.files(path = projDir,
                        pattern = 'csv',
                        recursive = TRUE,
                        full.names = TRUE)

num_fetch_vals <- 640000/32 # should be num rows in the fetch table/32 
combined_arctus_fetch <- read_csv("./processed_evd/fetch/combined_arctus_fetch.csv")

### issue: geometry column was read in as a character column, needs to be an sf object to apply st_coordinates 

# prepare combined_arctus_fetch for input into effective fetch function 
fetch_summary <- combined_arctus_fetch %>% 
  mutate(site_names = rep(1:num_fetch_vals, each = 32),
         X = st_coordinates(geometry)[,1],
         Y = st_coordinates(geometry)[,2]) %>% 
  dplyr::rename(fetch_length_m = fetch_length, cluster = cl_group) %>% 
  select(site_names,fetch_length_m,direction,X,Y)

#### Checkpoint: Effective fetch must be in the same CRS as wind rasters. Project to correct CRS below
wind_crs <- crs(wind_rast_mxspd360) # use existing wind raster to determine CRS
fetch_summary_projected <- st_transform(fetch_summary, wind_crs) # project to wind CRS


# Calculate minimum fetch (proxy for distance to coast) and sum fetch (proxy for wave exposure)
fetch_proxies <- fetch_summary %>% 
  group_by(geometry) %>% 
  mutate(min_fetch = min(fetch_length_m),
         sum_fetch = sum(fetch_length_m))


# Calculate effective fetch ----

fetch_eff <- fetch_summary %>%
  dplyr::select(X, Y, direction, fetch_length_m) %>% 
  group_by(X, Y) %>% 
  nest() %>% 
  bind_cols(., site_names = as.character(seq_len(nrow(.)))) %>%   
  # Calculate effective fetch and place in tibble
  mutate(fetch_eff = map(data, ~roll.recycle(.$fetch_length_m, 9, 8, by = 4)),
         fetch_eff = map(fetch_eff, ~as.vector((weights %*% .)/sum(weights))),
         fetch_eff = map(fetch_eff, ~tibble(direction = c(360, seq(45, 315, 45)), fetch = .))) %>% 
  st_as_sf(coords = c("X","Y"), crs = crs_fetch)

saveRDS(fetch_eff, "./processed_evd/fetch/fetch_effective.rds")

# Error in clean_columns(as.data.frame(obj), factorsAsCharacter) : 
# list columns are only allowed with raw vector contents

st_write(fetch_eff,
         dsn = './processed_evd/fetch/',
         layer = "fetch_effective",
         driver = 'ESRI Shapefile') 
