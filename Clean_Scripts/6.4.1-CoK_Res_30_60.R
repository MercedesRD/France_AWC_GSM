#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 30-60 cm
###
###  Add krigged residuals (from Hocine) to Cubist predictions
### with point predictions it seems very difficult. 
### Therefore, I will just resample the 90x90 grid predictions to the dimensions and extent of my cubist rasters

###  Author: Mercedes Roman Dobarco
###  Date: 18/12/2017

####### Load packages
library(sp)
library(maptools)
library(rgeos)
library(gstat)
library(rgdal)
library(geoR)
library(lattice)
library(spatstat)
library(automap)
library(ggplot2)
library(Hmisc)
library(plyr)
library(soiltexture)
library(Cubist)
library(foreach)
library(doParallel)
library(ithir)
library(raster)
library(snow)
library(doSNOW)
library(gdalUtils)

### Load my Cubist predictions
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/5.4-Cubist_preds")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/5.4-Cubist_preds")
silt.alr.cub.30_60 <- raster("silt.alr.cub.30_60.tif" )
argile.alr.cub.30_60 <- raster( "argile.alr.cub.30_60.tif")
coarse.cub.30_60 <- raster( "coarse.cub.30_60.tif")

NAvalue(argile.alr.cub.30_60)
NAvalue(silt.alr.cub.30_60)
NAvalue(coarse.cub.30_60)

# ### Clay-alr residuals --------------------------------------------------

### change directory to the import files
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
### Read the ascii file 90m x 90m grid
argile.alr.ckR.30_60 <- read.asciigrid("R_Clay_30_60.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# argile.alr.ckR.30_60 <- as(argile.alr.ckR.30_60, "SpatialPixelsDataFrame")
### Transform to raster
argile.alr.ckR.30_60 <- raster(argile.alr.ckR.30_60, layer=1, values=TRUE)
writeRaster(x = argile.alr.ckR.30_60, filename = "argile.alr.ckR.30_60.tif", format = "GTiff", overwrite = T )
plot(argile.alr.ckR.30_60)
NAvalue(argile.alr.ckR.30_60)


### Transform to my desired grid with gdalwarp (or ArcGIS if it does not work)
dir.create("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")

# gdal_setInstallation(ignore.full_scan=FALSE)
# 
# gdalwarp(srcfile = "D:/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/argile.alr.ckR.30_60.tif",
#          dstfile = "D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/argile.alr.ckR.30_60.tif",
#          s_srs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#          t_srs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#          te = c(99226,6049647,1099999,7110524),
#          te_srs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#          tr = c(89.99756, 89.99635),
#          r = "bilinear",
#          overwrite=TRUE,
#          verbose=TRUE)
# 
# ### In LINUX (virtual machine) I have installed version 1.11.2 of GDAL
# 
# gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/argile.alr.ckR.30_60.tif",
#          dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/argile.alr.ckR.30_60.tif",
#          s_srs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#          t_srs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#          te = c(99226,6049647,1099999,7110524),
#          tr = res(argile.alr.cub.30_60),
#          r = "bilinear",
#          overwrite=TRUE,
#          verbose=TRUE)


gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/argile.alr.ckR.30_60.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/argile.alr.ckR.30_60.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(argile.alr.cub.30_60),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

### remove the large raster stored in memory
rm(argile.alr.ckR.30_60)

argile.alr.ckR.30_60.r  <- raster("argile.alr.ckR.30_60.tif")
plot(argile.alr.ckR.30_60.r ) ### It worked!

argile_30_60 <- stack(argile.alr.cub.30_60 ,argile.alr.ckR.30_60.r)
plot(argile_30_60)


#  ### Clay-alr STD -----------------------------------------------------------

### now the STD of argile.alr

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
### Read the ascii file 90m x 90m grid
argile.alr.ckSTD.30_60 <- read.asciigrid("Std_R_Clay_30_60.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# argile.alr.ckR.30_60 <- as(argile.alr.ckR.30_60, "SpatialPixelsDataFrame")
### Transform to raster
argile.alr.ckSTD.30_60 <- raster(argile.alr.ckSTD.30_60, layer=1, values=TRUE)
writeRaster(x = argile.alr.ckSTD.30_60, filename = "argile.alr.ckSTD.30_60.tif", format = "GTiff", overwrite = T )
plot(argile.alr.ckSTD.30_60)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/argile.alr.ckSTD.30_60.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/argile.alr.ckSTD.30_60.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(argile.alr.cub.30_60),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

argile.alr.ckSTD.30_60.r  <- raster("argile.alr.ckSTD.30_60.tif")

argile_30_60 <- stack(argile_30_60 , argile.alr.ckSTD.30_60.r)
plot(argile_30_60)



# ### Silt-alr residuals --------------------------------------------------

### change directory to the import files
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
### Read the ascii file 90m x 90m grid
limon.alr.ckR.30_60 <- read.asciigrid("R_Silt_30_60.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# limon.alr.ckR.30_60 <- as(limon.alr.ckR.30_60, "SpatialPixelsDataFrame")
### Transform to raster
limon.alr.ckR.30_60 <- raster(limon.alr.ckR.30_60, layer=1, values=TRUE)
writeRaster(x = limon.alr.ckR.30_60, filename = "limon.alr.ckR.30_60.tif", format = "GTiff", overwrite = T )
plot(limon.alr.ckR.30_60)
NAvalue(limon.alr.ckR.30_60)

#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")


gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/limon.alr.ckR.30_60.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/limon.alr.ckR.30_60.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.30_60),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

### remove the large raster stored in memory
rm(limon.alr.ckR.30_60)

limon.alr.ckR.30_60.r  <- raster("limon.alr.ckR.30_60.tif")
plot(limon.alr.ckR.30_60.r ) ### It worked!

limon_30_60 <- stack(silt.alr.cub.30_60 ,limon.alr.ckR.30_60.r)
plot(limon_30_60)


# ### Silt-alr STD --------------------------------------------------------

### now the STD of limon.alr

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
### Read the ascii file 90m x 90m grid
limon.alr.ckSTD.30_60 <- read.asciigrid("Std_R_Silt_30_60.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# limon.alr.ckR.30_60 <- as(limon.alr.ckR.30_60, "SpatialPixelsDataFrame")
### Transform to raster
limon.alr.ckSTD.30_60 <- raster(limon.alr.ckSTD.30_60, layer=1, values=TRUE)
writeRaster(x = limon.alr.ckSTD.30_60, filename = "limon.alr.ckSTD.30_60.tif", format = "GTiff", overwrite = T )
plot(limon.alr.ckSTD.30_60)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/limon.alr.ckSTD.30_60.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/limon.alr.ckSTD.30_60.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.30_60),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

limon.alr.ckSTD.30_60.r  <- raster("limon.alr.ckSTD.30_60.tif")

limon_30_60 <- stack(limon_30_60 , limon.alr.ckSTD.30_60.r)
plot(limon_30_60)


### end of the script (untill I get coarse elements residuals)

# ### coarse elements --------------------------------------------------------

### now the Coarse elements

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
### Read the ascii file 90m x 90m grid
coarse.kR.30_60 <- read.asciigrid("R_EG_30_60.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
### Transform to raster
coarse.kR.30_60 <- raster(coarse.kR.30_60, layer=1, values=TRUE)
writeRaster(x = coarse.kR.30_60, filename = "coarse.kR.30_60.tif", format = "GTiff", overwrite = T )
plot(coarse.kR.30_60)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/coarse.kR.30_60.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/coarse.kR.30_60.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.30_60),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

coarse.kR.30_60.r  <- raster("coarse.kR.30_60.tif")


# ### Coarse  STD --------------------------------------------------------

### now the STD of coarse elements

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60")
### Read the ascii file 90m x 90m grid
coarse.kSTD.30_60 <- read.asciigrid("Std_R_EG_30_60.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
### Transform to raster
coarse.kSTD.30_60 <- raster(coarse.kSTD.30_60, layer=1, values=TRUE)
writeRaster(x = coarse.kSTD.30_60, filename = "coarse.kSTD.30_60.tif", format = "GTiff", overwrite = T )
plot(coarse.kSTD.30_60)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_30_60/coarse.kSTD.30_60.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60/coarse.kSTD.30_60.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.30_60),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

coarse.kSTD.30_60.r  <- raster("coarse.kSTD.30_60.tif")

coarse_30_60 <- stack(coarse.kR.30_60.r , coarse.kSTD.30_60.r)
plot(coarse_30_60)
