#############################################################################################################################################
###  Prediction of available water capacity in France - GSM depth intervals
###
###  In this script: IGCS data on particle size distribution for France, GSM depth intervals

######## 1. Create raster stack with candidate, predictor covariates
######## 3. Extraction of covariate values at IGCS profiles coordinates

###  Author: Mercedes Roman Dobarco
###  Date: 02/05/2017


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

### Set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"Clean_Output/3.2-Covariates"))

# 1. Texture and coarse elements data -------------------------------------------------

### load the data
load(paste0(HomeDir,"Clean_Output/3.1-Join_Texture_Gravel/soil_data.RData"))
head(soil_data)

###############################################################################################################################################


# ######## 2. Create raster stack with candidate, predictor covari --------

### Raster files are already cropped to the same extent.
### Load raster files
covarDir <- paste0(HomeDir,"Covariates/Tiff")

setwd(covarDir)
gtif.files <- list.files(path=".", pattern=".tif$", all.files=TRUE)

### Eliminate a couple variable I don't want, "pred.tif", "tmax.tif", "tmean.tif", and "class_aspect
gtif.files <- gtif.files[ gtif.files %nin% c("class_aspect", "prec.tif", "tmax.tif", "tmean.tif", "tmin.tif")]

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
rm(load_raster, gtif.files, prediction.list)


##################################################################################################################################

######## 2. Extraction of covariate values at IGCS coordinates

### IGCS data on texture for the GSM depth intervals

### spatial representation of data
proj4string(covariates) <- CRS("+init=epsg:2154")

### Extract covariates at IGCS sites
igcs_granulo_gsm <- extract(covariates, soil_data.sp, method="simple", sp=TRUE)

#### As dataframe
igcs_granulo_gsm <- as.data.frame(igcs_granulo_gsm)
igcs_granulo_gsm_bck <- igcs_granulo_gsm

### select only complete cases (for covariates)
igcs_granulo_gsm <- igcs_granulo_gsm[complete.cases(igcs_granulo_gsm[,c("bdforet","clc06","cti", "curv_long","curv_trans","curvature",
                                                                "ecoclim","eros","etpMax","etpMean" ,"etpMedian","etpMin","EVI_median_jan",
                                                                "EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000","hli","idpr",
                                                                "linear_aspect","mat11","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june",
                                                                "roughness","rrMax","rrMean","rrMedian","rrMin","sar" ,"scale_pos",
                                                                "slope","slopeascos","slopeassin","slopeastrasp","soil1","srr" , "srtm",
                                                                "tMeanMax","tMeanMean" ,"tMeanMedian" ,"tMeanMin" )]),]

### Check now number of observations by depth interval
sum(!is.na(igcs_granulo_gsm$argile_0_5)) ## 36159
sum(!is.na(igcs_granulo_gsm$argile_5_15)) ## 36108
sum(!is.na(igcs_granulo_gsm$argile_15_30)) ## 35401
sum(!is.na(igcs_granulo_gsm$argile_30_60)) ## 31494
sum(!is.na(igcs_granulo_gsm$argile_60_100)) ## 24849
sum(!is.na(igcs_granulo_gsm$argile_100_200)) ## 13086

sum(!is.na(igcs_granulo_gsm$coarse_0_5)) ## 51575
sum(!is.na(igcs_granulo_gsm$coarse_5_15)) ## 53153
sum(!is.na(igcs_granulo_gsm$coarse_15_30)) ## 53115
sum(!is.na(igcs_granulo_gsm$coarse_30_60)) ## 50110
sum(!is.na(igcs_granulo_gsm$coarse_60_100)) ## 47519
sum(!is.na(igcs_granulo_gsm$coarse_100_200)) ## 44802


### Transform some variables into factors
igcs_granulo_gsm$bdforet <-factor(igcs_granulo_gsm$bdforet) 
igcs_granulo_gsm$clc06 <-factor(igcs_granulo_gsm$clc06) 
igcs_granulo_gsm$ecoclim <-factor(igcs_granulo_gsm$ecoclim) 
igcs_granulo_gsm$mat11 <-factor(igcs_granulo_gsm$mat11)  
igcs_granulo_gsm$modelclc <-factor(igcs_granulo_gsm$modelclc) 
igcs_granulo_gsm$soil1 <-factor(igcs_granulo_gsm$soil1)  

### save image
save.image("igcs_granulo_gsm.RData")
###Export csv
write.csv2(igcs_granulo_gsm, file="igcs_granulo_gsm.csv")

################################################################################################################
#### end of the script