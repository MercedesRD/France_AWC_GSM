#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 0-5 cm
###
###  Pretreatment of cubist and co-kriging of the residuals
###  - Calculate kriging covariance
###  Author: Mercedes Roman Dobarco
###  Date: 28/03/2018


##################################################################################################################

### Clean workspace
rm(list=ls())

# Load packages -----------------------------------------------------------
### Load packages
require(raster)
require(sp)
require(rgdal)
library(gstat)
library(maptools)
library(rgeos)
library(geoR)
library(lattice)
library(ggplot2)
library(Hmisc)
library(plyr)
library(soiltexture)
library(Cubist)
library(foreach)
library(doParallel)
library(ithir)
library(snow)
library(doSNOW)

# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")
# argile.alr.limon.alr.cov.0_5 <- raster("argile.alr.limon.alr.cov.0_5.tif")
# plot(argile.alr.limon.alr.cov.0_5,breaks=c(seq(-2,16, by=1)),
#      col=rev(viridis_pal(option="D")(18)))
# 
# 
# dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.7-LocalCovariance")
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.7-LocalCovariance")
# 
# ###### Calculate Local mean
# 
# ### test raster example
# r <- raster(ncols=10, nrows=10, xmn=0)
# values(r) <- 1:ncell(r)
# plot(r)
# ### Give some NA
# r <- calc(r, fun=function(x){ ifelse(x %in% c(22,10,35,34,48,56,73,84,96,98), NA, x)} )
# plot(r)
# # circle filter for square cells
# gf <- focalWeight(r, 40, "circle");gf
# gf[gf>0] <- 1 
# u <- focal(r, w=gf, fun=mean, na.rm=TRUE, pad=TRUE); plot(u) ### adding rows around to avoid the edge effect
# rm(gf, u, r)
# 
# ########################################################################################################################################
# 
# ## Load the kriging residuals and STD
 setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.1-CoK_Res_0_5")
 argile.alr.ckR.0_5 <- raster("argile.alr.ckR.0_5.tif"); plot(argile.alr.ckR.0_5)
 limon.alr.ckR.0_5 <- raster("limon.alr.ckR.0_5.tif"); plot(limon.alr.ckR.0_5)

setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.7-LocalCovariance")
# range <- 95050  ### Half the range of the covariogram (example, 190 km)
# 
# # circle filter for square cells
# gf <- focalWeight(argile.alr.ckR.0_5, range, "circle")
# gf[gf>0] <- 1 
# 
# argile.alr.ckR.0_5.mean <- focal(argile.alr.ckR.0_5, w=gf, fun=mean, na.rm=TRUE, pad=TRUE,
#                                  filename="argile.alr.ckR.0_5.mean.tif", overwrite=TRUE); plot(argile.alr.ckR.0_5.mean)
# limon.alr.ckR.0_5.mean <- focal(limon.alr.ckR.0_5, w=gf, fun=mean, na.rm=TRUE, pad=TRUE,
#                                 filename="limon.alr.ckR.0_5.mean.tif", overwrite=TRUE); plot(limon.alr.ckR.0_5.mean)
# 

##This operation is taking SOOOOO LOOOONG, that I did it in GRASS

### to see the original files, open GRASS in your computer, set as GRASS GIS database directory
### Y:\communs_infosol\Projets\GlobalSoilMap.net\GIS

### Select GRASS location: RGF93_2154
### Select GRASS Mapset: France

### and start GRASS session
#### the region settings can be checked with the command:
# g.region -p  

# projection: 99 (Lambert Conformal Conic)
# zone:       0
# datum:      towgs84=0,0,0,0,0,0,0
# ellipsoid:  grs80
# north:      7110524
# south:      6049647
# west:       99226
# east:       1242375
# nsres:      89.99635222
# ewres:      89.99755944
# rows:       11788
# cols:       12702
# cells:      149731176

### 1. import into GRASS region FRANCE

### 2. Calculate focal MEAN, with CIRCULAR MASK. the number or neighbours is 111 for an aproximate diameter of 10 km
# r.neighbors -c --overwrite           input=argile@France output=argile.alr.ckR.mean size=111
# r.neighbors -c --overwrite --verbose input=limon@France output=limon.alr.ckR.mean size=111

### 3. Export from GRASS
# r.out.gdal --overwrite input=argile.alr.ckR.mean@France output=D:\romandobarco\AWC\AWC_GSM_Dec2017\Output\6.7-LocalCovariance\argile.alr.ckR.0_5.fmean.tif format=GTiff
# r.out.gdal --overwrite input=limon.alr.ckR.mean@France output=D:\romandobarco\AWC\AWC_GSM_Dec2017\Output\6.7-LocalCovariance\limon.alr.ckR.0_5.fmean.tif format=GTiff

### Read in R
argile.alr.ckR.0_5.fmean <- raster("argile.alr.ckR.0_5.fmean.tif"); plot(argile.alr.ckR.0_5.fmean)
limon.alr.ckR.0_5.fmean <- raster("limon.alr.ckR.0_5.fmean.tif"); plot(limon.alr.ckR.0_5.fmean)

### CROP
argile.alr.ckR.0_5.mean <- crop(argile.alr.ckR.0_5.fmean, extent(argile.alr.ckR.0_5),
                                 filename="argile.alr.ckR.0_5.mean.tif",type="GTiff", overwrite=T)
limon.alr.ckR.0_5.mean <- crop(limon.alr.ckR.0_5.fmean, extent(limon.alr.ckR.0_5),
                                filename="limon.alr.ckR.0_5.mean.tif",type="GTiff", overwrite=T)

### Transform nan into NA
# argile.alr.ckR.0_5.fmean <- calc(argile.alr.ckR.0_5.fmean, function(x) {x[is.nan(x)] <- NA; return(x) } ,
#                                            na.rm=T, inf.rm=T, filename="argile.alr.ckR.0_5.fmean.na.tif",type="GTiff", overwrite=T)
# limon.alr.ckR.0_5.fmean <- calc(limon.alr.ckR.0_5.fmean, function(x) {x[is.nan(x)] <- NA; return(x) } ,
#                                 na.rm=T,inf.rm=T, filename="limon.alr.ckR.0_5.fmean.na.tif",type="GTiff", overwrite=T)

par(mfrow=c(1,3))
plot(argile.alr.ckR.0_5.mean)
plot(limon.alr.ckR.0_5.mean)

### Create raster stack
cov.s <- stack( argile.alr.ckR.0_5, limon.alr.ckR.0_5, argile.alr.ckR.0_5.mean, limon.alr.ckR.0_5.mean)

### write function
cov.sp <- function(x){(x[1]-x[3])*(x[2]-x[4])}

beginCluster(4)
argile.alr.limon.alr.fmean.cov.0_5 <- clusterR(cov.s, calc, args = list(cov.sp),
                       filename = "argile.alr.limon.alr.fmean.cov.0_5.tif", format = "GTiff",
                       na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
par(mfrow=c(1,1))
plot(argile.alr.limon.alr.fmean.cov.0_5)


###############################################################################################################################

### Clean workspace
rm(list=ls())

## Load the kriging residuals and STD
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.2-CoK_Res_5_15")
argile.alr.ckR.5_15 <- raster("argile.alr.ckR.5_15.tif")
limon.alr.ckR.5_15 <- raster("limon.alr.ckR.5_15.tif")

### Create raster stack
granulo.5_15.s <- stack( argile.alr.ckR.5_15,limon.alr.ckR.5_15)

### extract mean
clay.mean <- mean(getValues(argile.alr.ckR.5_15), na.rm=TRUE)
limon.mean <- mean(getValues(limon.alr.ckR.5_15), na.rm=TRUE)

### Calculate covariance

### sum cubist predictions and kriged residuals
#dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")

### write function
cov.sp <- function(x){(x[1]-clay.mean)*(x[2]-limon.mean)}

beginCluster(7)
argile.alr.limon.alr.cov.5_15 <- clusterR(granulo.5_15.s, calc, args = list(cov.sp),export=c("clay.mean", "limon.mean"),
                                         filename = "argile.alr.limon.alr.cov.5_15.tif", format = "GTiff",
                                         na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.alr.limon.alr.cov.5_15)


###############################################################################################################################

### Clean workspace
rm(list=ls())

## Load the kriging residuals and STD
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-CoK_Res_15_30")
argile.alr.ckR.15_30 <- raster("argile.alr.ckR.15_30.tif")
limon.alr.ckR.15_30 <- raster("limon.alr.ckR.15_30.tif")

### Create raster stack
granulo.15_30.s <- stack( argile.alr.ckR.15_30,limon.alr.ckR.15_30)

### extract mean
clay.mean <- mean(getValues(argile.alr.ckR.15_30), na.rm=TRUE)
limon.mean <- mean(getValues(limon.alr.ckR.15_30), na.rm=TRUE)

### Calculate covariance

### sum cubist predictions and kriged residuals
#dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")

### write function
cov.sp <- function(x){(x[1]-clay.mean)*(x[2]-limon.mean)}

beginCluster(7)
argile.alr.limon.alr.cov.15_30 <- clusterR(granulo.15_30.s, calc, args = list(cov.sp),export=c("clay.mean", "limon.mean"),
                                          filename = "argile.alr.limon.alr.cov.15_30.tif", format = "GTiff",
                                          na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.alr.limon.alr.cov.15_30)


###############################################################################################################################

### Clean workspace
rm(list=ls())

## Load the kriging residuals and STD
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.4-CoK_Res_30_60")
argile.alr.ckR.30_60 <- raster("argile.alr.ckR.30_60.tif")
limon.alr.ckR.30_60 <- raster("limon.alr.ckR.30_60.tif")

### Create raster stack
granulo.30_60.s <- stack( argile.alr.ckR.30_60,limon.alr.ckR.30_60)

### extract mean
clay.mean <- mean(getValues(argile.alr.ckR.30_60), na.rm=TRUE)
limon.mean <- mean(getValues(limon.alr.ckR.30_60), na.rm=TRUE)

### Calculate covariance

### sum cubist predictions and kriged residuals
#dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")

### write function
cov.sp <- function(x){(x[1]-clay.mean)*(x[2]-limon.mean)}

beginCluster(7)
argile.alr.limon.alr.cov.30_60 <- clusterR(granulo.30_60.s, calc, args = list(cov.sp),export=c("clay.mean", "limon.mean"),
                                           filename = "argile.alr.limon.alr.cov.30_60.tif", format = "GTiff",
                                           na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.alr.limon.alr.cov.30_60)


###############################################################################################################################

### Clean workspace
rm(list=ls())

## Load the kriging residuals and STD
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-CoK_Res_60_100")
argile.alr.ckR.60_100 <- raster("argile.alr.ckR.60_100.tif")
limon.alr.ckR.60_100 <- raster("limon.alr.ckR.60_100.tif")

### Create raster stack
granulo.60_100.s <- stack( argile.alr.ckR.60_100,limon.alr.ckR.60_100)

### extract mean
clay.mean <- mean(getValues(argile.alr.ckR.60_100), na.rm=TRUE)
limon.mean <- mean(getValues(limon.alr.ckR.60_100), na.rm=TRUE)

### Calculate covariance

### sum cubist predictions and kriged residuals
#dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")

### write function
cov.sp <- function(x){(x[1]-clay.mean)*(x[2]-limon.mean)}

beginCluster(7)
argile.alr.limon.alr.cov.60_100 <- clusterR(granulo.60_100.s, calc, args = list(cov.sp),export=c("clay.mean", "limon.mean"),
                                           filename = "argile.alr.limon.alr.cov.60_100.tif", format = "GTiff",
                                           na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.alr.limon.alr.cov.60_100)

###############################################################################################################################

### Clean workspace
rm(list=ls())

## Load the kriging residuals and STD
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.6-CoK_Res_100_200")
argile.alr.ckR.100_200 <- raster("argile.alr.ckR.100_200.tif")
limon.alr.ckR.100_200 <- raster("limon.alr.ckR.100_200.tif")

### Create raster stack
granulo.100_200.s <- stack( argile.alr.ckR.100_200,limon.alr.ckR.100_200)

### extract mean
clay.mean <- mean(getValues(argile.alr.ckR.100_200), na.rm=TRUE)
limon.mean <- mean(getValues(limon.alr.ckR.100_200), na.rm=TRUE)

### Calculate covariance

### sum cubist predictions and kriged residuals
#dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")

### write function
cov.sp <- function(x){(x[1]-clay.mean)*(x[2]-limon.mean)}

beginCluster(7)
argile.alr.limon.alr.cov.100_200 <- clusterR(granulo.100_200.s, calc, args = list(cov.sp),export=c("clay.mean", "limon.mean"),
                                            filename = "argile.alr.limon.alr.cov.100_200.tif", format = "GTiff",
                                            na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.alr.limon.alr.cov.100_200)


### End of script