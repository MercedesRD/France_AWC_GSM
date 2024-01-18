###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date 05/03/2018

### Objective: Improving the predictions of coarse elements

#### Task: EXTRACT new classification for parent material at calibration data (coarse elements)

######## 1. Import shapefile
######## 2. Extract data from 2 layers


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

### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"Input/CoarseElements/BDGSF"))
MAT_RF <- readOGR(dsn=".", layer="soil2_L93")

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
####### Just do it easier, Opoen the calibration data
load("EG.data_0_5.r.RData")

##### convert to spatial
EG.data_0_5.sp <- EG.data_0_5.r
coordinates(EG.data_0_5.sp) <- ~ x+y
proj4string(EG.data_0_5.sp) <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"

### Extract the data at the claibration locations
EG.data_0_5.soil <- over(EG.data_0_5.sp, MAT_RF, returnList = FALSE)
### Bind
EG.data_0_5.r <- cbind(EG.data_0_5.r,EG.data_0_5.soil[,c("MAT_RF", "AGL_RF")])

### check the boxplots
EG.data_0_5.r$MAT_RF <- as.factor(EG.data_0_5.r$MAT_RF)
EG.data_0_5.r$AGL_RF <- as.factor(EG.data_0_5.r$AGL_RF)

### Plot
boxplot(exp(EG.data_0_5.r$coarse_0_5) ~ EG.data_0_5.r$MAT_RF, main= "Reclassification of MAT12")
boxplot(EG.data_0_5.r$coarse_0_5 ~ EG.data_0_5.r$MAT_RF)

boxplot(exp(EG.data_0_5.r$coarse_0_5) ~ EG.data_0_5.r$AGL_RF,  main= "Reclassification of AGLIM1")
boxplot(EG.data_0_5.r$coarse_0_5 ~ EG.data_0_5.r$AGL_RF)

boxplot(exp(EG.data_0_5.r$coarse_0_5) ~ EG.data_0_5.r$mat11,  main= "MAT11")
boxplot(EG.data_0_5.r$coarse_0_5 ~ EG.data_0_5.r$mat11)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements"))
dir.create(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
save(EG.data_0_5.r, file="EG.data_0_5.MAT_RF.RData")

########################################################################################################################

### now layer 5-15
####### Just do it easier, Opoen the calibration data
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
load("EG.data_5_15.r.RData")

##### convert to spatial
EG.data_5_15.sp <- EG.data_5_15.r
coordinates(EG.data_5_15.sp) <- ~ x+y
proj4string(EG.data_5_15.sp) <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"

### Extract the data at the claibration locations
EG.data_5_15.soil <- over(EG.data_5_15.sp, MAT_RF, returnList = FALSE)
### Bind
EG.data_5_15.r <- cbind(EG.data_5_15.r,EG.data_5_15.soil[,c("MAT_RF", "AGL_RF")])

### check the boxplots
EG.data_5_15.r$MAT_RF <- as.factor(EG.data_5_15.r$MAT_RF)
EG.data_5_15.r$AGL_RF <- as.factor(EG.data_5_15.r$AGL_RF)

### Plot
boxplot(exp(EG.data_5_15.r$coarse_5_15) ~ EG.data_5_15.r$MAT_RF, main= "Reclassification of MAT12")
boxplot(EG.data_5_15.r$coarse_5_15 ~ EG.data_5_15.r$MAT_RF)

boxplot(exp(EG.data_5_15.r$coarse_5_15) ~ EG.data_5_15.r$AGL_RF,  main= "Reclassification of AGLIM1")
boxplot(EG.data_5_15.r$coarse_5_15 ~ EG.data_5_15.r$AGL_RF)

boxplot(exp(EG.data_5_15.r$coarse_5_15) ~ EG.data_5_15.r$mat11,  main= "MAT11")
boxplot(EG.data_5_15.r$coarse_5_15 ~ EG.data_5_15.r$mat11)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
save(EG.data_5_15.r, file="EG.data_5_15.MAT_RF.RData")

########################################################################################################################

### now layer 15-30
####### Just do it easier, Opoen the calibration data
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
load("EG.data_15_30.r.RData")

##### convert to spatial
EG.data_15_30.sp <- EG.data_15_30.r
coordinates(EG.data_15_30.sp) <- ~ x+y
proj4string(EG.data_15_30.sp) <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"

### Extract the data at the claibration locations
EG.data_15_30.soil <- over(EG.data_15_30.sp, MAT_RF, returnList = FALSE)
### Bind
EG.data_15_30.r <- cbind(EG.data_15_30.r,EG.data_15_30.soil[,c("MAT_RF", "AGL_RF")])

### check the boxplots
EG.data_15_30.r$MAT_RF <- as.factor(EG.data_15_30.r$MAT_RF)
EG.data_15_30.r$AGL_RF <- as.factor(EG.data_15_30.r$AGL_RF)

### Plot
boxplot(exp(EG.data_15_30.r$coarse_15_30) ~ EG.data_15_30.r$MAT_RF, main= "Reclassification of MAT12")
boxplot(EG.data_15_30.r$coarse_15_30 ~ EG.data_15_30.r$MAT_RF)

boxplot(exp(EG.data_15_30.r$coarse_15_30) ~ EG.data_15_30.r$AGL_RF,  main= "Reclassification of AGLIM1")
boxplot(EG.data_15_30.r$coarse_15_30 ~ EG.data_15_30.r$AGL_RF)

boxplot(exp(EG.data_15_30.r$coarse_15_30) ~ EG.data_15_30.r$mat11,  main= "MAT11")
boxplot(EG.data_15_30.r$coarse_15_30 ~ EG.data_15_30.r$mat11)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
save(EG.data_15_30.r, file="EG.data_15_30.MAT_RF.RData")


########################################################################################################################

### now layer 30-60
####### Just do it easier, Opoen the calibration data
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
load("EG.data_30_60.r.RData")

##### convert to spatial
EG.data_30_60.sp <- EG.data_30_60.r
coordinates(EG.data_30_60.sp) <- ~ x+y
proj4string(EG.data_30_60.sp) <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"

### Extract the data at the claibration locations
EG.data_30_60.soil <- over(EG.data_30_60.sp, MAT_RF, returnList = FALSE)
### Bind
EG.data_30_60.r <- cbind(EG.data_30_60.r,EG.data_30_60.soil[,c("MAT_RF", "AGL_RF")])

### check the boxplots
EG.data_30_60.r$MAT_RF <- as.factor(EG.data_30_60.r$MAT_RF)
EG.data_30_60.r$AGL_RF <- as.factor(EG.data_30_60.r$AGL_RF)

### Plot
boxplot(exp(EG.data_30_60.r$coarse_30_60) ~ EG.data_30_60.r$MAT_RF, main= "Reclassification of MAT12")
boxplot(EG.data_30_60.r$coarse_30_60 ~ EG.data_30_60.r$MAT_RF)

boxplot(exp(EG.data_30_60.r$coarse_30_60) ~ EG.data_30_60.r$AGL_RF,  main= "Reclassification of AGLIM1")
boxplot(EG.data_30_60.r$coarse_30_60 ~ EG.data_30_60.r$AGL_RF)

boxplot(exp(EG.data_30_60.r$coarse_30_60) ~ EG.data_30_60.r$mat11,  main= "MAT11")
boxplot(EG.data_30_60.r$coarse_30_60 ~ EG.data_30_60.r$mat11)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
save(EG.data_30_60.r, file="EG.data_30_60.MAT_RF.RData")

########################################################################################################################

### now layer 30-60
####### Just do it easier, Opoen the calibration data
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
load("EG.data_30_60.r.RData")

##### convert to spatial
EG.data_30_60.sp <- EG.data_30_60.r
coordinates(EG.data_30_60.sp) <- ~ x+y
proj4string(EG.data_30_60.sp) <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"

### Extract the data at the claibration locations
EG.data_30_60.soil <- over(EG.data_30_60.sp, MAT_RF, returnList = FALSE)
### Bind
EG.data_30_60.r <- cbind(EG.data_30_60.r,EG.data_30_60.soil[,c("MAT_RF", "AGL_RF")])

### check the boxplots
EG.data_30_60.r$MAT_RF <- as.factor(EG.data_30_60.r$MAT_RF)
EG.data_30_60.r$AGL_RF <- as.factor(EG.data_30_60.r$AGL_RF)

### Plot
boxplot(exp(EG.data_30_60.r$coarse_30_60) ~ EG.data_30_60.r$MAT_RF, main= "Reclassification of MAT12")
boxplot(EG.data_30_60.r$coarse_30_60 ~ EG.data_30_60.r$MAT_RF)

boxplot(exp(EG.data_30_60.r$coarse_30_60) ~ EG.data_30_60.r$AGL_RF,  main= "Reclassification of AGLIM1")
boxplot(EG.data_30_60.r$coarse_30_60 ~ EG.data_30_60.r$AGL_RF)

boxplot(exp(EG.data_30_60.r$coarse_30_60) ~ EG.data_30_60.r$mat11,  main= "MAT11")
boxplot(EG.data_30_60.r$coarse_30_60 ~ EG.data_30_60.r$mat11)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
save(EG.data_30_60.r, file="EG.data_30_60.MAT_RF.RData")

########################################################################################################################

### now layer 60-100
####### Just do it easier, Open the calibration data
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
load("EG.data_60_100.r.RData")

##### convert to spatial
EG.data_60_100.sp <- EG.data_60_100.r
coordinates(EG.data_60_100.sp) <- ~ x+y
proj4string(EG.data_60_100.sp) <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"

### Extract the data at the claibration locations
EG.data_60_100.soil <- over(EG.data_60_100.sp, MAT_RF, returnList = FALSE)
### Bind
EG.data_60_100.r <- cbind(EG.data_60_100.r,EG.data_60_100.soil[,c("MAT_RF", "AGL_RF")])

### check the boxplots
EG.data_60_100.r$MAT_RF <- as.factor(EG.data_60_100.r$MAT_RF)
EG.data_60_100.r$AGL_RF <- as.factor(EG.data_60_100.r$AGL_RF)

### Plot
boxplot(exp(EG.data_60_100.r$coarse_60_100) ~ EG.data_60_100.r$MAT_RF, main= "Reclassification of MAT12")
boxplot(EG.data_60_100.r$coarse_60_100 ~ EG.data_60_100.r$MAT_RF)

boxplot(exp(EG.data_60_100.r$coarse_60_100) ~ EG.data_60_100.r$AGL_RF,  main= "Reclassification of AGLIM1")
boxplot(EG.data_60_100.r$coarse_60_100 ~ EG.data_60_100.r$AGL_RF)

boxplot(exp(EG.data_60_100.r$coarse_60_100) ~ EG.data_60_100.r$mat11,  main= "MAT11")
boxplot(EG.data_60_100.r$coarse_60_100 ~ EG.data_60_100.r$mat11)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
save(EG.data_60_100.r, file="EG.data_60_100.MAT_RF.RData")

########################################################################################################################

### now layer 100-200
####### Just do it easier, Opoen the calibration data
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
load("EG.data_100_200.r.RData")

##### convert to spatial
EG.data_100_200.sp <- EG.data_100_200.r
coordinates(EG.data_100_200.sp) <- ~ x+y
proj4string(EG.data_100_200.sp) <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"

### Extract the data at the claibration locations
EG.data_100_200.soil <- over(EG.data_100_200.sp, MAT_RF, returnList = FALSE)
### Bind
EG.data_100_200.r <- cbind(EG.data_100_200.r,EG.data_100_200.soil[,c("MAT_RF", "AGL_RF")])

### check the boxplots
EG.data_100_200.r$MAT_RF <- as.factor(EG.data_100_200.r$MAT_RF)
EG.data_100_200.r$AGL_RF <- as.factor(EG.data_100_200.r$AGL_RF)

### Plot
boxplot(exp(EG.data_100_200.r$coarse_100_200) ~ EG.data_100_200.r$MAT_RF, main= "Reclassification of MAT12")
boxplot(EG.data_100_200.r$coarse_100_200 ~ EG.data_100_200.r$MAT_RF)

boxplot(exp(EG.data_100_200.r$coarse_100_200) ~ EG.data_100_200.r$AGL_RF,  main= "Reclassification of AGLIM1")
boxplot(EG.data_100_200.r$coarse_100_200 ~ EG.data_100_200.r$AGL_RF)

boxplot(exp(EG.data_100_200.r$coarse_100_200) ~ EG.data_100_200.r$mat11,  main= "MAT11")
boxplot(EG.data_100_200.r$coarse_100_200 ~ EG.data_100_200.r$mat11)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.2-CalDatParentMaterial"))
save(EG.data_100_200.r, file="EG.data_100_200.MAT_RF.RData")

#### end of the script