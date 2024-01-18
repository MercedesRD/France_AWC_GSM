
### Onjective: Create a mask for each GSM layer based on estimated soil depth by Lacoste et al. (2016)

### Load soil depth in mm
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"
soil_depth_mm <- raster(paste0(HomeDir,"Input/soil_depth_mm.tif"))
dir.create(paste0(HomeDir, "13-soil_depth_masks"))

setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Soil_depth")
### Mask 
mask_0_5 <- calc(x = soil_depth_mm, fun = function(x)  {x[x <= 0 | is.nan(x)] <- NA; return(x) } , filename="mask_0_5.tif",type="GTiff", overwrite=T);plot(mask_0_5)
mask_5_15 <- calc(x = soil_depth_mm, fun = function(x)  {x[x <= 50 | is.nan(x)] <- NA; return(x) } , filename="mask_5_15.tif",type="GTiff", overwrite=T);plot(mask_5_15)
mask_15_30 <- calc(x = soil_depth_mm, fun = function(x)  {x[x <= 150 | is.nan(x)] <- NA; return(x) } , filename="mask_15_30.tif",type="GTiff", overwrite=T);plot(mask_15_30)
mask_30_60 <- calc(x = soil_depth_mm, fun = function(x)  {x[x <= 300 | is.nan(x)] <- NA; return(x) } , filename="mask_30_60.tif",type="GTiff", overwrite=T);plot(mask_30_60)
mask_60_100 <- calc(x = soil_depth_mm, fun = function(x)  {x[x <= 600 | is.nan(x)] <- NA; return(x) } , filename="mask_60_100.tif",type="GTiff", overwrite=T);plot(mask_60_100)
mask_100_200 <- calc(x = soil_depth_mm, fun = function(x)  {x[x <= 1000 | is.nan(x)] <- NA; return(x) } , filename="mask_100_200.tif",type="GTiff", overwrite=T);plot(mask_100_200)
save(mask_0_5, mask_5_15, mask_15_30, mask_30_60, mask_60_100, mask_100_200, file="soil_depth_mask.RData")
