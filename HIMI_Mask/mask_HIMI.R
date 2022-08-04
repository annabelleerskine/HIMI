##################################################################################################################################
### Extract HIMI fishery region 200-2000 m for Julia and Romain to get temperature time series for size spectrum models
##################################################################################################################################


## 1) Set up ----
library(raster)
library(sf)

path <- "C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Toothfish FRDC\\Analysis\\Environmental_Layers\\HIMI\\"

# EEZ shapefile
EEZ<- st_read(paste0 (path, "Raw_data\\HIMI_EEZ_polygon.shp"))

## 2) Get bathy and crop to EEZ and mask to deth range
bathy<-brick(paste0(path, "Processed_data\\HIMI_bathy_deriv"))
bathy<-raster::subset(bathy, subset="aadc_bathy")

##Depth range
bathy[bathy< -2000] <-NA
bathy[bathy > - 200] <-NA

#EEZ
bathy_eez<-crop(bathy, EEZ)
bathy_eez<-mask(bathy_eez, EEZ)

#save
writeRaster(bathy_eez, filename ="C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Toothfish FRDC\\Size_spectra\\Bathy_mask")


## 3) Extract time series of SST using this mask
sst<-brick(paste0(path, "Processed_data\\sst_yr_mth"))
sst<-crop(sst, extent(67.05, 79.45,-56.55, -49.95 ))

bathy_match<-crop(bathy_eez, extent(67.05, 79.45,-56.55, -49.95 ))
sst<-mask (sst, bathy_match)
writeRaster(sst, filename ="C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Toothfish FRDC\\Size_spectra\\SST")


## 4) Extract time series of BRAN SST using this mask
bran<-brick(paste0(path, "Processed_data\\BRAN20_temp_sfc"))
bran<-crop(bran, extent(67.05, 79.45,-56.55, -49.95 ))
bran<-mask (bran, bathy_match)
writeRaster(bran, filename ="C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Toothfish FRDC\\Size_spectra\\BRAN_surface_temp")


