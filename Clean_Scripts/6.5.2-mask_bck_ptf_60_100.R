#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 60-100 cm
###
###  Pretreatment of cubist and co-kriging of the residuals
###  - sum to obtain final predictions
###  - Backtransform to obtain sand and clay
###  - Apply PTF for calculating SMFC and SMPWP
###  - Subtract the volume of coarse elements
###  Author: Mercedes Roman Dobarco
###  Date: 16/05/2018


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

### set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

# ### Create stack with cubist and co-kriged residuals --------------------

### Load raster files
### Load my Cubist predictions
setwd(paste0(HomeDir,"Clean_Output/5-Cubist_preds"))
silt.alr.cub.60_100 <- raster("silt.alr.cub.60_100.tif" )
argile.alr.cub.60_100 <- raster( "argile.alr.cub.60_100.tif")

### Load the kriging residuals and STD
setwd(paste0(HomeDir,"Clean_Output/6.5-CoK_Res_60_100"))
argile.alr.ckR.60_100 <- raster("argile.alr.ckR.60_100.tif")
argile.alr.ckSTD.60_100  <- raster("argile.alr.ckSTD.60_100.tif")
limon.alr.ckR.60_100 <- raster("limon.alr.ckR.60_100.tif")
limon.alr.ckSTD.60_100  <- raster("limon.alr.ckSTD.60_100.tif")

### Load quantile random forest predictions for coarse elements, the log and the backtransformed
setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
coarse.qrF.mean.60_100 <- raster( "coarse_60_100.preds.mean.tif")
coarse.qrF.sd.60_100 <- raster( "coarse_60_100.preds.sd.tif")
coarse.qrF.05.60_100 <- raster( "coarse_60_100.preds.05p.tif")
coarse.qrF.95.60_100 <- raster( "coarse_60_100.preds.95p.tif")

# ### Load the covariance between clay.alr and silt.slr
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.5-covariance")
# argile.alr.limon.alr.cov.60_100 <- raster("argile.alr.limon.alr.cov.60_100.tif")

### Create raster stacks
granulo.60_100.s <- stack(argile.alr.cub.60_100, silt.alr.cub.60_100,
                        argile.alr.ckR.60_100, limon.alr.ckR.60_100,
                       argile.alr.ckSTD.60_100, limon.alr.ckSTD.60_100)
coarse.60_100.s <- stack(coarse.qrF.mean.60_100, coarse.qrF.sd.60_100, coarse.qrF.05.60_100,coarse.qrF.95.60_100)

### sum cubist predictions and kriged residuals
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))

## write function
sum_clay <- function(x){x[1]+x[3]}

beginCluster(7)
argile.alr.preds.60_100 <- clusterR(granulo.60_100.s, calc, args = list(sum_clay),
                                 filename = "argile.alr.preds.60_100.tif", format = "GTiff",
                                 na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.alr.preds.60_100)

sum_silt <- function(x){x[2]+x[4]}

beginCluster(7)
limon.alr.preds.60_100 <- clusterR(granulo.60_100.s, calc, args = list(sum_silt),
                                filename = "limon.alr.preds.60_100.tif", format = "GTiff",
                                na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(limon.alr.preds.60_100)

argile.alr.preds.60_100 <- raster("argile.alr.preds.60_100.tif")
limon.alr.preds.60_100 <- raster("limon.alr.preds.60_100.tif")

### Put in the raster stack
granulo.60_100.s <- stack(granulo.60_100.s, argile.alr.preds.60_100,limon.alr.preds.60_100)

### Load Back transformation to sand, silt, and clay ------------------------------
### Perform ALR back-transformation (in %)

### Calculate CLAY - RK

### Write function
bck_clay <- function(x) {(exp(x[7])/(1+exp(x[7])+exp(x[8])))*100}

beginCluster(7)
argile.preds.RK.60_100 <- clusterR(granulo.60_100.s, calc, args = list(bck_clay),
                                filename = "argile.preds.RK.60_100.tif", format = "GTiff",
                                na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.preds.RK.60_100)


### Calculate SILT - RK

### Write function
bck_silt <- function(x) {(exp(x[8])/(1+exp(x[7])+exp(x[8])))*100}

beginCluster(7)
limon.preds.RK.60_100 <- clusterR(granulo.60_100.s, calc, args = list(bck_silt),
                               filename = "limon.preds.RK.60_100.tif", format = "GTiff",
                               na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(limon.preds.RK.60_100)


### Calculate SAND - RK

### Write function
bck_sand <- function(x) {(1/(1+exp(x[7])+exp(x[8])))*100}

beginCluster(7)
sable.preds.RK.60_100 <- clusterR(granulo.60_100.s, calc, args = list(bck_sand),
                               filename = "sable.preds.RK.60_100.tif", format = "GTiff",
                               na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(sable.preds.RK.60_100)

argile.preds.RK.60_100 <- raster("argile.preds.RK.60_100.tif")
limon.preds.RK.60_100 <- raster("limon.preds.RK.60_100.tif")
sable.preds.RK.60_100 <- raster("sable.preds.RK.60_100.tif")

### Add to the stack
granulo.60_100.s <- stack(granulo.60_100.s, argile.preds.RK.60_100,limon.preds.RK.60_100,sable.preds.RK.60_100)
names(granulo.60_100.s) <- c("argile.alr.cub.60_100" ,  "silt.alr.cub.60_100" ,    "argile.alr.ckR.60_100" ,
                           "limon.alr.ckR.60_100"  ,  "argile.alr.ckSTD.60_100", "limon.alr.ckSTD.60_100",
                           "argile.alr.preds.60_100", "limon.alr.preds.60_100" ,
                          "clay" , "silt" ,  "sand") ##3 change the name for the lm function

# ### Apply PTF for calculating SMFC and SMPWP -----------------------------------------------------------
### Load PTF
### load("D:/romandobarco/AWC/AWC_GSM/Input/PTF_SndCly_NS.RData")
load(paste0(HomeDir,"Input/ContinuousPTF.RData"))
### Remove the PTF that I don't need
rm( "lm_w20_ClSdBd"  ,      "lm_w20_ClSdBdSOC"  ,   "lm_w20_ClSdSOC"    ,   "lm_w20_sub_ClSd"   ,  
    "lm_w20_sub_ClSdBd" ,   "lm_w20_sub_ClSdBdSOC" ,"lm_w20_sub_ClSdSOC",  "lm_w20_top_ClSd"  ,    "lm_w20_top_ClSdBd" ,   "lm_w20_top_ClSdBdSOC",
    "lm_w20_top_ClSdSOC" ,  "lm_w25_ClSd"      ,    "lm_w25_ClSdBd"    ,    "lm_w25_ClSdBdSOC" ,    "lm_w25_ClSdSOC"   ,    "lm_w25_sub_ClSd"  ,   
    "lm_w25_sub_ClSdBd"   , "lm_w25_sub_ClSdBdSOC" ,"lm_w25_sub_ClSdSOC" ,  "lm_w25_top_ClSd",      "lm_w25_top_ClSdBd" ,   "lm_w25_top_ClSdBdSOC",
    "lm_w25_top_ClSdSOC"   ,     "lm_w42_ClSdBd"   ,     "lm_w42_ClSdBdSOC"  ,   "lm_w42_ClSdSOC"   ,    "lm_w42_sub_ClSd" ,   
    "lm_w42_sub_ClSdBd"   , "lm_w42_sub_ClSdBdSOC" ,"lm_w42_sub_ClSdSOC" ,  "lm_w42_top_ClSd" ,     "lm_w42_top_ClSdBd"   , "lm_w42_top_ClSdBdSOC",
    "lm_w42_top_ClSdSOC" )

### Apply non-stratified PTF for w20
tic <- Sys.time()
beginCluster(7)
SMFC_60_100 <- clusterR(granulo.60_100.s, predict, args = list(lm_w20_ClSd),
                     filename = "SMFC_60_100.tif", format = "GTiff",
                     na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()

SMFC_60_100 <- raster("SMFC_60_100.tif")
plot(SMFC_60_100)

### Apply non-stratified PTF for w42
beginCluster(7)
SMPWP_60_100 <- clusterR(granulo.60_100.s, predict, args = list(lm_w42_ClSd),
                      filename = "SMPWP_60_100.tif", format = "GTiff",
                      na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
tac <- Sys.time()
tac-tic
SMPWP_60_100 <- raster("SMPWP_60_100.tif")
plot(SMPWP_60_100)

### Calculate difference
awc.60_100.s <- stack(SMFC_60_100, SMPWP_60_100)
ru_vol_60_100 <- raster("ru_vol_60_100.tif")

difference <- function(x){x[1]-x[2]}

beginCluster(7)
ru_vol_60_100 <- clusterR(awc.60_100.s, calc, args = list(difference),
                      filename = "ru_vol_60_100.tif", format = "GTiff",
                      na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(ru_vol_60_100)

awc.60_100.s <- stack(awc.60_100.s, ru_vol_60_100)

### Multiply by layer depth
ru_mm <- function(x){100*x}

beginCluster(7)
ru_mm_60_100 <- clusterR(ru_vol_60_100, calc, args = list(ru_mm),
                       filename = "ru_mm_60_100.tif", format = "GTiff",
                       na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(ru_mm_60_100)
ru_mm_60_100 <- raster("ru_mm_60_100.tif")
awc.60_100.s <- stack(awc.60_100.s, ru_mm_60_100)

### attach the coarse elements
awc.60_100.s <- stack(awc.60_100.s, coarse.qrF.mean.60_100)

### subtract the volume of coarse elements to the AWC (volume)
ru_coarse <- function(x){x[3]*(1-(x[5]/100))*100} ## in mm, taking into account coarse elements

beginCluster(7)
awc_mm_60_100 <- clusterR(awc.60_100.s, calc, args = list(ru_coarse),
                      filename = "awc_mm_60_100.tif", format = "GTiff",
                      na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(awc_mm_60_100)
awc_mm_60_100 <- raster("awc_mm_60_100.tif")
awc.60_100.s <- stack(awc.60_100.s, awc_mm_60_100)
plot(awc.60_100.s)

### Save again
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
save.image("6.5.2.ru_60_100.RData")

#### end of script