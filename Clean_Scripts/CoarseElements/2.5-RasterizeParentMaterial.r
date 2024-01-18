###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date 06/03/2018

### Objective: Improving the predictions of coarse elements

#### Task: Rasterize the shapefile with Parent material classification

######## 1. Import shapefile
######## 2. Rasterize at the same extent, resolution, etc , than the other variables

####### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(raster)
library(lattice)
library(spatstat)
library(automap)
library(ggplot2)
library(plyr)
library(doParallel)
library(foreach)

### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"Clean_Output/CoarseElements"))
load("2.3-CovariatesAddChelsa.RData")

setwd(paste0(HomeDir,"Input/CoarseElements/BDGSF"))
MAT_RF <- readOGR(dsn=".", layer="soil2_L93")
MAT_RF <- spTransform(MAT_RF, proj4string(covariates.chelsa))

### copy raster
r <- covariates.chelsa[[1]]
values(r) <- NA

MAT_RF@data$MAT_RF.num <- as.numeric(as.character(MAT_RF@data$MAT_RF))
MAT_RF@data$AGL_RF.num <- as.numeric(as.character(MAT_RF@data$AGL_RF))

sort(unique(MAT_RF@data$AGL_RF.num))
sort(unique(levels(MAT_RF@data$AGL_RF)))
sort(unique(MAT_RF@data$MAT_RF.num))
sort(unique(levels(MAT_RF@data$MAT_RF)))

dir.create(paste0(HomeDir,"Input/CoarseElements/BDGSF/GeoTiff"))
setwd(paste0(HomeDir,"Input/CoarseElements/BDGSF/GeoTiff"))
MAT_RF.r <- rasterize(MAT_RF, r,field=MAT_RF$MAT_RF.num, filename = "MAT_RF.tif", format = "GTiff", overwrite = T, silent=FALSE)
plot(MAT_RF.r)
AGL_RF.r <- rasterize(MAT_RF, r,field=MAT_RF$AGL_RF.num, filename = "AGL_RF.tif", format = "GTiff", overwrite = T, silent=FALSE)
plot(AGL_RF.r)

### Create a raster stack with all predictors

### we now from previous scripts that all leve;s of soil1 and ecoclim are not present in the calibration data, so we masked them out ina  different GeoTiff
### Replace some rasters that were corrected for the models
covariates4maps <- covariates.chelsa

setwd(paste0(HomeDir,"Covariates/Prediction_tiff_EG"))
covariates4maps[["clc06"]] <- raster("clc06.tif")
covariates4maps[["ecoclim"]] <- raster("ecoclim.tif")
covariates4maps[["soil1"]] <- raster("soil1.tif")

### eliminate "clc_06" from covariates4maps
covariates4maps <- dropLayer(covariates4maps, 2)

predictors <-c("bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax","etpMean",
  "etpMedian","etpMin","EVI_median_jan","EVI_median_june","exposition","graviAL2000","graviBg2000",
  "graviG2000","hli","idpr","linear_aspect","mat11","modelclc","mrrtf","mrvbf","NDVI_median_jan",
  "NDVI_median_june","roughness","rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos",
  "slopeassin","slopeastrasp","soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin",
  "CHELSA_bio_1","CHELSA_bio_10","CHELSA_bio_11","CHELSA_bio_12","CHELSA_bio_13","CHELSA_bio_14",
  "CHELSA_bio_16","CHELSA_bio_17","CHELSA_bio_18","CHELSA_bio_19","CHELSA_bio_5","CHELSA_bio_6",
  "CHELSA_bio_8","CHELSA_bio_9","MAT_RF","AGL_RF")
setdiff(names(covariates4maps), predictors)
setdiff(predictors, names(covariates4maps))

### add the last two predictors
### Change names first
MAT_RF <- MAT_RF.r
AGL_RF <- AGL_RF.r
covariates4maps <- stack(covariates4maps,MAT_RF)
covariates4maps <- stack(covariates4maps,AGL_RF)

rm(MAT_RF.r, AGL_RF.r, MAT_RF, AGL_RF)

#### Before using it for predicting, I need to check that all levels are present in the calibration data,
#### and that the variables that are factors are considered factors
#### So, some calc and some is.factor? stuff
setwd(paste0(HomeDir,"Clean_Output/CoarseElements"))
save.image("2.PredictorStack.RData")

#### End of the script