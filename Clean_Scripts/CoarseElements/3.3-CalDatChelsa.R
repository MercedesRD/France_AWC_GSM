###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date 06/03/2018

### Objective: Improving the predictions of coarse elements

#### Task: Extract Chelsa climate data for all GSM depths calibration data

######## 1. Import raster stack
######## 2. Extract data from stack and save

####### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(raster)
library(lattice)
library(ggplot2)
library(plyr)
library(doParallel)
library(foreach)

### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
load("EG.data_0_5.MAT_RF.RData")
load("EG.data_5_15.MAT_RF.RData")
load("EG.data_15_30.MAT_RF.RData")
load("EG.data_30_60.MAT_RF.RData")
load("EG.data_60_100.MAT_RF.RData")
load("EG.data_100_200.MAT_RF.RData")
load(paste0(HomeDir,"Clean_Output/CoarseElements/2.3-CovariatesAddChelsa.RData"))

### My list of dataframes
df.input <- list(EG.data_0_5.r, EG.data_5_15.r, EG.data_15_30.r, EG.data_30_60.r, EG.data_60_100.r, EG.data_100_200.r)

### Perform the extraction in parallel to save some time

cl <- makeCluster(6)
registerDoParallel(cl)

df.output <- foreach(i = 1:length(df.input), .packages = c("sp","raster"), .export="chelsa.s") %dopar%{
    ## Copy file
    df.input.i <- df.input[[i]]
    coordinates(df.input.i) <- ~ x+y
    proj4string(df.input.i) <- "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"
    ### Extract covariates at IGCS sites
    df.out.i <- extract(chelsa.s, df.input.i, method="simple", sp=FALSE)
    df.out.i <- as.data.frame(df.out.i)
    #### bind to the dataframe of origin
    df.out <- cbind(df.input[[i]],df.out.i )
    return(df.out)
    }

stopCluster(cl)

### save image
str(df.output)
rm(i,cl)
names(df.output) <- c("EG.data_0_5.r", "EG.data_5_15.r", "EG.data_15_30.r", "EG.data_30_60.r", "EG.data_60_100.r", "EG.data_100_200.r")

setwd(paste0(HomeDir,"Input"))
save.image("CalibrationDataCoarse.RData")
##### end of the script