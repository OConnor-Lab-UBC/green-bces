library(terra)

chi_raster <- rast("environmental_dat/Chissassibi_SDB_Sentinel2.tif")
wem_raster <- rast("environmental_dat/Wemendji_SDB_Sentinel2.tif")
wem_raster
chi_raster

# ranging from 0 to 3m depth
wem_bathy <- as.data.frame(wem_raster, xy = T)
chi_bathy <- as.data.frame(chi_raster, xy=T)

joint_bathy <- sprc(wem_raster, chi_raster)
merged_bathy_raster <- merge(joint_bathy)
plot(merged_bathy_raster)
