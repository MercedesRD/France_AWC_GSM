###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date: 01/03/2018

######## Objective: Bring in new climate covariates, CHELSA
######## 1. Crop covariates at 90 m --------

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
library(Hmisc)
library(plyr)

### set working directory

### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

# ######## 1. Crop covariates at 90 m --------

### Directory with all covariates but MODIS
covar.dir <- paste0(HomeDir,"Covariates/CHELSA/90m")
setwd(covar.dir)

#### List with files in Gtiff directory
gtif.files <- list.files(path=".", pattern="CHELSA", all.files=TRUE)

#### This is the function to load rasters
load_raster <- function (x) {
    maps <- list()
    for (rast in 1:length(x)) {  
        maps[[rast]] <- raster(x[rast])
    }
    return(maps)
}

prediction.list <- load_raster(gtif.files)

#### Rename raster files
names(prediction.list) <- c("CHELSA_bio_1","CHELSA_bio_10","CHELSA_bio_11","CHELSA_bio_12","CHELSA_bio_13","CHELSA_bio_14",
                            "CHELSA_bio_16","CHELSA_bio_17","CHELSA_bio_18","CHELSA_bio_19","CHELSA_bio_5","CHELSA_bio_6", 
                            "CHELSA_bio_8","CHELSA_bio_9")

covar.stack <- stack(prediction.list)
#sample_raster <- raster("slope.tif")

### Crop to exclude Corse
desired.extent <- extent(covar.stack, 1, 11788, 1, 11120) 

### Try to crop a stack
#   covar.stack <- crop(covar.stack, desired.extent) ### Takes TOOOOO LONG
CropDirectory <- paste0(HomeDir,"Covariates/CHELSA/crop")
dir.create(CropDirectory)
setwd(CropDirectory)

crop.list <- list()
for (i in 1: length(prediction.list)){
    print(i)
    calc(prediction.list[[i]], function(x) {x[is.nan(x)] <- NA; return(x) }) ### Mask nan values
    crop.list[[i]] <- crop(prediction.list[[i]], desired.extent,
                           filename= paste(names(prediction.list[i]), ".tif", sep=""), type="GTiff", overwrite=T )
    gc()
}

#### Rename raster files
names(crop.list) <- c("CHELSA_bio_1","CHELSA_bio_10","CHELSA_bio_11","CHELSA_bio_12","CHELSA_bio_13","CHELSA_bio_14",
                      "CHELSA_bio_16","CHELSA_bio_17","CHELSA_bio_18","CHELSA_bio_19","CHELSA_bio_5","CHELSA_bio_6", 
                      "CHELSA_bio_8","CHELSA_bio_9")


chelsa.s <- stack(crop.list)
plot(chelsa.s)

rm(load_raster, gtif.files, prediction.list, crop.list, i, covar.dir,desired.extent)

### Add these raster files to the predictor raster stacks (and for extraction of the values at calibration dataset locations)

################################################################################################################
#### end of the script