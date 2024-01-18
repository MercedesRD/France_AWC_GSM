###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

#################  Date: 20/02/2018

### Objective: Improving the predictions of coarse elements

### Task: 1. Bring in new climate covariates, CHELSA
###       2. Cropping the maps, ONLY CROPPING, in the original projection and resolution 

# All products of CHELSA are in a geographic coordinate system referenced to the WGS 84 horizontal datum,
# with the horizontal coordinates expressed in decimal degrees. The CHELSA layer extents
# (minimum and maximum latitude and longitude) are a result of the coordinate system inherited from 
# the 1-arc-second GMTED2010 data which itself inherited the grid extent from the 1-arc-second SRTM data.

############################################### Project, resample, and crop CHELSA covariates

################### CHELSA CLIMATE VARIABLES

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
library(foreach)
library(doParallel)
library(raster)
library(snow)
library(doSNOW)
library(gdalUtils)

setwd("~/romandobarco/AWC/Covariates/CHELSA")

chelsa.files <- list.files(pattern="CHELSA")

for(i in 1:length(chelsa.files)){
    print(i)
    ### Crop to a smaller extent
    gdalwarp(srcfile = paste0("/home/mercedes/romandobarco/AWC/Covariates/CHELSA/",chelsa.files[i]),
             dstfile = gsub(".tif","_crop.tif", paste0("/home/mercedes/romandobarco/AWC/Covariates/CHELSA/",chelsa.files[i])),
             te = c(-9.166806239, 39.16652726, 13.33319511, 54.83319456),
             overwrite=TRUE,
             verbose=TRUE)
}

### The next steps, in GRASS, consisted in importing into a new REGION and LOCATION, 
### change to FRANCE EPSG:2154 location,reproject the rasters, and export them as Geo Tiff

### end of the script