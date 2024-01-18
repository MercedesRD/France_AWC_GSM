###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date: 26/02/2018

### Objective: Improving the predictions of coarse elements

### Task: 1. Exclude sondage from Normandy from all calibration datasets

####### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(lattice)
library(ggplot2)
library(plyr)
library(Cubist)
library(foreach)
library(doParallel)
library(raster)


# ### GSM 0-5 cm ----------------------------------------------------------

### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"clean_Output/4.Cubist"))
load("4.1-Cubist_0_5.RData")

### Eliminate the sondages from Normandy

### Read the sites
setwd(paste0(HomeDir,"Input/CoarseElements"))
id_profil_normandy <- read.csv("id_profil_sondage_normandie.csv", sep=";")
id_profil_normandy$id_profil <- as.factor(id_profil_normandy$id_profil)

##3 Exclude the Normady sondage profiles
my_profiles <- setdiff(EG.data_0_5$id_profil, id_profil_normandy$id_profil)
EG.data_0_5.r <- EG.data_0_5[EG.data_0_5$id_profil %in% my_profiles, ]

## Only want complete.cases
EG.data_0_5.r <- EG.data_0_5.r[complete.cases(EG.data_0_5.r),]
### drop levels of factors
EG.data_0_5.r[] <- lapply(EG.data_0_5.r, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements"))
dir.create("3.1-CalibrationData")
setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
save(EG.data_0_5.r, file="EG.data_0_5.r.RData")
rm(list=ls())

# ### GSM 5-15 cm ----------------------------------------------------------
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"clean_Output/4.Cubist"))
load("4.2-Cubist_5_15.RData")

### Redo without the sondage form Normandy

### Read the sites
setwd(paste0(HomeDir,"Input/CoarseElements"))
id_profil_normandy <- read.csv("id_profil_sondage_normandie.csv", sep=";")
id_profil_normandy$id_profil <- as.factor(id_profil_normandy$id_profil)

##3 Exclude the Normady sondage profiles
my_profiles <- setdiff(EG.data_5_15$id_profil, id_profil_normandy$id_profil)
EG.data_5_15.r <- EG.data_5_15[EG.data_5_15$id_profil %in% my_profiles, ]

## Only want complete.cases
EG.data_5_15.r <- EG.data_5_15.r[complete.cases(EG.data_5_15.r),]
### drop levels of factors
EG.data_5_15.r[] <- lapply(EG.data_5_15.r, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
save(EG.data_5_15.r, file="EG.data_5_15.r.RData")
rm(list=ls())

# ### GSM 15-30 cm ----------------------------------------------------------
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"clean_Output/4.Cubist"))
load("4.3-Cubist_15_30.RData")

### Redo without the sondage form Normandy

### Read the sites
setwd(paste0(HomeDir,"Input/CoarseElements"))
### Read the sites
id_profil_normandy <- read.csv("id_profil_sondage_normandie.csv", sep=";")
id_profil_normandy$id_profil <- as.factor(id_profil_normandy$id_profil)

##3 Exclude the Normady sondage profiles
my_profiles <- setdiff(EG.data_15_30$id_profil, id_profil_normandy$id_profil)
EG.data_15_30.r <- EG.data_15_30[EG.data_15_30$id_profil %in% my_profiles, ]

## Only want complete.cases
EG.data_15_30.r <- EG.data_15_30.r[complete.cases(EG.data_15_30.r),]
### drop levels of factors
EG.data_15_30.r[] <- lapply(EG.data_15_30.r, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
save(EG.data_15_30.r, file="EG.data_15_30.r.RData")
rm(list=ls())

# ### GSM 30-60 cm ----------------------------------------------------------

### Redo without the sondage form Normandy
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"clean_Output/4.Cubist"))
### Eliminate the sondages from Normandy
load("4.4-Cubist_30_60.RData")

### Read the sites
setwd(paste0(HomeDir,"Input/CoarseElements"))
id_profil_normandy <- read.csv("id_profil_sondage_normandie.csv", sep=";")
id_profil_normandy$id_profil <- as.factor(id_profil_normandy$id_profil)

##3 Exclude the Normady sondage profiles
my_profiles <- setdiff(EG.data_30_60$id_profil, id_profil_normandy$id_profil)
EG.data_30_60.r <- EG.data_30_60[EG.data_30_60$id_profil %in% my_profiles, ]

## Only want complete.cases
EG.data_30_60.r <- EG.data_30_60.r[complete.cases(EG.data_30_60.r),]
### drop levels of factors
EG.data_30_60.r[] <- lapply(EG.data_30_60.r, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
save(EG.data_30_60.r, file="EG.data_30_60.r.RData")
rm(list=ls())


# ### GSM 60-100 cm ----------------------------------------------------------
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"clean_Output/4.Cubist"))
### Eliminate the sondages from Normandy
load("4.5-Cubist_60_100.RData")

### Read the sites
setwd(paste0(HomeDir,"Input/CoarseElements"))
id_profil_normandy <- read.csv("id_profil_sondage_normandie.csv", sep=";")
id_profil_normandy$id_profil <- as.factor(id_profil_normandy$id_profil)

##3 Exclude the Normady sondage profiles
my_profiles <- setdiff(EG.data_60_100$id_profil, id_profil_normandy$id_profil)
EG.data_60_100.r <- EG.data_60_100[EG.data_60_100$id_profil %in% my_profiles, ]

## Only want complete.cases
EG.data_60_100.r <- EG.data_60_100.r[complete.cases(EG.data_60_100.r),]
### drop levels of factors
EG.data_60_100.r[] <- lapply(EG.data_60_100.r, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
save(EG.data_60_100.r, file="EG.data_60_100.r.RData")
rm(list=ls())


# ### GSM 100-200 cm ----------------------------------------------------------
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"clean_Output/4.Cubist"))

### Eliminate the sondages from Normandy
load("4.6-Cubist_100_200.RData")

### Read the sites
setwd(paste0(HomeDir,"Input/CoarseElements"))
id_profil_normandy <- read.csv("id_profil_sondage_normandie.csv", sep=";")
id_profil_normandy$id_profil <- as.factor(id_profil_normandy$id_profil)

##3 Exclude the Normady sondage profiles
my_profiles <- setdiff(EG.data_100_200$id_profil, id_profil_normandy$id_profil)
EG.data_100_200.r <- EG.data_100_200[EG.data_100_200$id_profil %in% my_profiles, ]

## Only want complete.cases
EG.data_100_200.r <- EG.data_100_200.r[complete.cases(EG.data_100_200.r),]
### drop levels of factors
EG.data_100_200.r[] <- lapply(EG.data_100_200.r, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

setwd(paste0(HomeDir,"Clean_Output/CoarseElements/3.1-CalibrationData"))
save(EG.data_100_200.r, file="EG.data_100_200.r.RData")
rm(list=ls())

##### end of the script