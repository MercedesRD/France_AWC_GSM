#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 60-100 cm
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
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/5.5-Cubist_preds")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/5.5-Cubist_preds")
silt.alr.cub.60_100 <- raster("silt.alr.cub.60_100.tif" )
argile.alr.cub.60_100 <- raster( "argile.alr.cub.60_100.tif")
coarse.cub.60_100 <- raster( "coarse.cub.60_100.tif")

NAvalue(argile.alr.cub.60_100)
NAvalue(silt.alr.cub.60_100)
NAvalue(coarse.cub.60_100)

# ### Clay-alr residuals --------------------------------------------------

### change directory to the import files
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
### Read the ascii file 90m x 90m grid
argile.alr.ckR.60_100 <- read.asciigrid("R_Clay_60_100.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# argile.alr.ckR.60_100 <- as(argile.alr.ckR.60_100, "SpatialPixelsDataFrame")
### Transform to raster
argile.alr.ckR.60_100 <- raster(argile.alr.ckR.60_100, layer=1, values=TRUE)
writeRaster(x = argile.alr.ckR.60_100, filename = "argile.alr.ckR.60_100.tif", format = "GTiff", overwrite = T )
plot(argile.alr.ckR.60_100)
NAvalue(argile.alr.ckR.60_100)


### Transform to my desired grid with gdalwarp (or ArcGIS if it does not work)
dir.create("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")

# gdal_setInstallation(ignore.full_scan=FALSE)
# 
# gdalwarp(srcfile = "D:/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/argile.alr.ckR.60_100.tif",
#          dstfile = "D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/argile.alr.ckR.60_100.tif",
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
# gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/argile.alr.ckR.60_100.tif",
#          dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/argile.alr.ckR.60_100.tif",
#          s_srs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#          t_srs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
#          te = c(99226,6049647,1099999,7110524),
#          tr = res(argile.alr.cub.60_100),
#          r = "bilinear",
#          overwrite=TRUE,
#          verbose=TRUE)


gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/argile.alr.ckR.60_100.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/argile.alr.ckR.60_100.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(argile.alr.cub.60_100),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

### remove the large raster stored in memory
rm(argile.alr.ckR.60_100)

argile.alr.ckR.60_100.r  <- raster("argile.alr.ckR.60_100.tif")
plot(argile.alr.ckR.60_100.r ) ### It worked!

argile_60_100 <- stack(argile.alr.cub.60_100 ,argile.alr.ckR.60_100.r)
plot(argile_60_100)


#  ### Clay-alr STD -----------------------------------------------------------

### now the STD of argile.alr

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
### Read the ascii file 90m x 90m grid
argile.alr.ckSTD.60_100 <- read.asciigrid("Std_R_Clay_60_100.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# argile.alr.ckR.60_100 <- as(argile.alr.ckR.60_100, "SpatialPixelsDataFrame")
### Transform to raster
argile.alr.ckSTD.60_100 <- raster(argile.alr.ckSTD.60_100, layer=1, values=TRUE)
writeRaster(x = argile.alr.ckSTD.60_100, filename = "argile.alr.ckSTD.60_100.tif", format = "GTiff", overwrite = T )
plot(argile.alr.ckSTD.60_100)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/argile.alr.ckSTD.60_100.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/argile.alr.ckSTD.60_100.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(argile.alr.cub.60_100),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

argile.alr.ckSTD.60_100.r  <- raster("argile.alr.ckSTD.60_100.tif")

argile_60_100 <- stack(argile_60_100 , argile.alr.ckSTD.60_100.r)
plot(argile_60_100)



# ### Silt-alr residuals --------------------------------------------------

### change directory to the import files
#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
### Read the ascii file 90m x 90m grid
limon.alr.ckR.60_100 <- read.asciigrid("R_Silt_60_100.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# limon.alr.ckR.60_100 <- as(limon.alr.ckR.60_100, "SpatialPixelsDataFrame")
### Transform to raster
limon.alr.ckR.60_100 <- raster(limon.alr.ckR.60_100, layer=1, values=TRUE)
writeRaster(x = limon.alr.ckR.60_100, filename = "limon.alr.ckR.60_100.tif", format = "GTiff", overwrite = T )
plot(limon.alr.ckR.60_100)
NAvalue(limon.alr.ckR.60_100)

#setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")
setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")


gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/limon.alr.ckR.60_100.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/limon.alr.ckR.60_100.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.60_100),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

### remove the large raster stored in memory
rm(limon.alr.ckR.60_100)

limon.alr.ckR.60_100.r  <- raster("limon.alr.ckR.60_100.tif")
plot(limon.alr.ckR.60_100.r ) ### It worked!

limon_60_100 <- stack(silt.alr.cub.60_100 ,limon.alr.ckR.60_100.r)
plot(limon_60_100)


# ### Silt-alr STD --------------------------------------------------------

### now the STD of limon.alr

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
### Read the ascii file 90m x 90m grid
limon.alr.ckSTD.60_100 <- read.asciigrid("Std_R_Silt_60_100.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
# ### Transform to SpatialPixelsDataframe
# limon.alr.ckR.60_100 <- as(limon.alr.ckR.60_100, "SpatialPixelsDataFrame")
### Transform to raster
limon.alr.ckSTD.60_100 <- raster(limon.alr.ckSTD.60_100, layer=1, values=TRUE)
writeRaster(x = limon.alr.ckSTD.60_100, filename = "limon.alr.ckSTD.60_100.tif", format = "GTiff", overwrite = T )
plot(limon.alr.ckSTD.60_100)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/limon.alr.ckSTD.60_100.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/limon.alr.ckSTD.60_100.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.60_100),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

limon.alr.ckSTD.60_100.r  <- raster("limon.alr.ckSTD.60_100.tif")

limon_60_100 <- stack(limon_60_100 , limon.alr.ckSTD.60_100.r)
plot(limon_60_100)


### end of the script (untill I get coarse elements residuals)

# ### coarse elements --------------------------------------------------------

### now the Coarse elements

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
### Read the ascii file 90m x 90m grid
coarse.kR.60_100 <- read.asciigrid("R_EG_60_100.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
### Transform to raster
coarse.kR.60_100 <- raster(coarse.kR.60_100, layer=1, values=TRUE)
writeRaster(x = coarse.kR.60_100, filename = "coarse.kR.60_100.tif", format = "GTiff", overwrite = T )
plot(coarse.kR.60_100)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/coarse.kR.60_100.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/coarse.kR.60_100.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.60_100),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

coarse.kR.60_100.r  <- raster("coarse.kR.60_100.tif")


# ### Coarse  STD --------------------------------------------------------

### now the STD of coarse elements

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100")
### Read the ascii file 90m x 90m grid
coarse.kSTD.60_100 <- read.asciigrid("Std_R_EG_60_100.asc", proj4string ="+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ) 
### Transform to raster
coarse.kSTD.60_100 <- raster(coarse.kSTD.60_100, layer=1, values=TRUE)
writeRaster(x = coarse.kSTD.60_100, filename = "coarse.kSTD.60_100.tif", format = "GTiff", overwrite = T )
plot(coarse.kSTD.60_100)

setwd("~/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")

gdalwarp(srcfile = "/home/mercedes/romandobarco//AWC/AWC_GSM_Dec2017/Import/Results_Layer_60_100/coarse.kSTD.60_100.tif",
         dstfile = "/home/mercedes/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100/coarse.kSTD.60_100.tif",
         te = c(99226,6049647,1099999,7110524),
         tr = res(silt.alr.cub.60_100),
         r = "bilinear",
         overwrite=TRUE,
         verbose=TRUE)

coarse.kSTD.60_100.r  <- raster("coarse.kSTD.60_100.tif")

coarse_60_100 <- stack(coarse.kR.60_100.r , coarse.kSTD.60_100.r)
plot(coarse_60_100)
