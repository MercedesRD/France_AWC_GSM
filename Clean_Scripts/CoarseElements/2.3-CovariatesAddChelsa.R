###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date: 01/03/2018

######## Objective: Bring in new climate covariates, CHELSA
######## 1. Add CHELSA covariates to the other raster stack, including clc06
######## 2. Create raster stack with candidate, predictor covariates - CHELSA

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

###############################################################################################################################################

# ######## 1. Create raster stack with candidate, predictor covari --------

### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

### Raster files are already cropped to the same extent.
### Load raster files
covarDir <- paste0(HomeDir,"Covariates/Tiff")

setwd(covarDir)
gtif.files <- list.files(path=".", pattern=".tif$", all.files=TRUE)

### Eliminate a couple variable I don't want, "pred.tif", "tmax.tif", "tmean.tif", and "class_aspect
gtif.files <- gtif.files[ gtif.files %nin% c("class_aspect.tif", "prec.tif", "tmax.tif", "tmean.tif", "tmin.tif")]

#### This is the function to load rasters
load_raster <- function (x) {
    maps <- list()
    for (rast in 1:length(x)) {  
        maps[[rast]] <- raster(x[rast])
    }
    return(maps)
}

prediction.list <- load_raster(gtif.files)

names(prediction.list) <- c("bdforet","clc06","cti", "curv_long","curv_trans","curvature",
                            "ecoclim","eros","etpMax","etpMean" ,"etpMedian","etpMin",
                            "EVI_median_jan","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                            "hli","idpr","linear_aspect","mat11","modelclc","mrrtf",
                            "mrvbf","NDVI_median_jan","NDVI_median_june","roughness","rrMax","rrMean",
                            "rrMedian","rrMin","sar" ,"scale_pos","slope","slopeascos",
                            "slopeassin","slopeastrasp","soil1","srr" , "srtm","tMeanMax",
                            "tMeanMean" ,"tMeanMedian" ,"tMeanMin" )

covariates <- stack(prediction.list)
plot(covariates)
rm(gtif.files, prediction.list)

### Add CHELSA covariates

CropDirectory <- paste0(HomeDir, "Covariates/CHELSA/crop")
setwd(CropDirectory)

gtif.files <- list.files(path=".", pattern="CHELSA", all.files=TRUE)
prediction.list <- load_raster(gtif.files)
chelsa.s <- stack(prediction.list)

### STACK TOGETHER
covariates.chelsa <- stack(covariates, chelsa.s)

### Clean and save
rm(covarDir, CropDirectory, gtif.files, prediction.list, load_raster)
setwd(paste0(HomeDir,"Clean_Output/CoarseElements"))
save.image("2.3-CovariatesAddChelsa.RData")

### during the development of the work, this was organized in a different folder:
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData")
# dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/2.3-CovariatesAddChelsa")
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/2.3-CovariatesAddChelsa")
# save.image("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/2.3-CovariatesAddChelsa/2.3-CovariatesAddChelsa.RData")

################################################################################################################
#### end of the script