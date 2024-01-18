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
###  - sum to btain final predictions
###  - Backtransform to obtain sand and clay
###  - Apply PTF for calculating SMFC and SMPWP
###  - Subtract the volume of coarse elements
###  Author: Mercedes Roman Dobarco
###  Date: 05/01/2018


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

# Load raster files
### Load my Cubist predictions
setwd(paste0(HomeDir,"Clean_Output/5-Cubist_preds"))
silt.alr.cub.0_5 <- raster("silt.alr.cub.0_5.tif" )
argile.alr.cub.0_5 <- raster( "argile.alr.cub.0_5.tif")

### Load the kriging residuals and STD
setwd(paste0(HomeDir,"Clean_Output/6.1-CoK_Res_0_5"))
argile.alr.ckR.0_5 <- raster("argile.alr.ckR.0_5.tif")
argile.alr.ckSTD.0_5  <- raster("argile.alr.ckSTD.0_5.tif")
limon.alr.ckR.0_5 <- raster("limon.alr.ckR.0_5.tif")
limon.alr.ckSTD.0_5  <- raster("limon.alr.ckSTD.0_5.tif")

### Load quantile random forest predictions for coarse elements, the log and the backtransformed
setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
coarse.qrF.mean.0_5 <- raster( "coarse_0_5.preds.mean.tif")
coarse.qrF.sd.0_5 <- raster( "coarse_0_5.preds.sd.tif")
coarse.qrF.05.0_5 <- raster( "coarse_0_5.preds.05p.tif")
coarse.qrF.95.0_5 <- raster( "coarse_0_5.preds.95p.tif")

# ### Load the covariance between clay.alr and silt.slr
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/6.3-covariance")
# argile.alr.limon.alr.cov.0_5 <- raster("argile.alr.limon.alr.cov.0_5.tif")
## finally We decided not to use the covariance calculated a posteriori

### Create raster stacks
granulo.0_5.s <- stack(argile.alr.cub.0_5, silt.alr.cub.0_5,
                       argile.alr.ckR.0_5, limon.alr.ckR.0_5,
                       argile.alr.ckSTD.0_5, limon.alr.ckSTD.0_5)
plot(granulo.0_5.s)                
coarse.0_5.s <- stack(coarse.qrF.mean.0_5, coarse.qrF.sd.0_5, coarse.qrF.05.0_5,coarse.qrF.95.0_5)
plot(coarse.0_5.s)

### sum cubist predictions and kriged residuals
dir.create(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))

## write function
sum_clay <- function(x){x[1]+x[3]}

beginCluster(7)
argile.alr.preds.0_5 <- clusterR(granulo.0_5.s, calc, args = list(sum_clay),
                       filename = "argile.alr.preds.0_5.tif", format = "GTiff",
                       na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.alr.preds.0_5)

sum_silt <- function(x){x[2]+x[4]}

beginCluster(7)
limon.alr.preds.0_5 <- clusterR(granulo.0_5.s, calc, args = list(sum_silt),
                                 filename = "limon.alr.preds.0_5.tif", format = "GTiff",
                                 na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(limon.alr.preds.0_5)

argile.alr.preds.0_5 <- raster("argile.alr.preds.0_5.tif")
limon.alr.preds.0_5 <- raster("limon.alr.preds.0_5.tif")

### Put in the raster stack
granulo.0_5.s <- stack(granulo.0_5.s, argile.alr.preds.0_5,limon.alr.preds.0_5)
names(granulo.0_5.s)

# ### Back transform to sand, silt, and clay ------------------------------

### Perform ALR back-transformation (in %)

### Calculate CLAY - RK

### Write function
bck_clay <- function(x) {(exp(x[7])/(1+exp(x[7])+exp(x[8])))*100}

beginCluster(7)
argile.preds.RK.0_5 <- clusterR(granulo.0_5.s, calc, args = list(bck_clay),
                                filename = "argile.preds.RK.0_5.tif", format = "GTiff",
                                na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(argile.preds.RK.0_5)


### Calculate SILT - RK

### Write function
bck_silt <- function(x) {(exp(x[8])/(1+exp(x[7])+exp(x[8])))*100}

beginCluster(7)
limon.preds.RK.0_5 <- clusterR(granulo.0_5.s, calc, args = list(bck_silt),
                                filename = "limon.preds.RK.0_5.tif", format = "GTiff",
                                na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(limon.preds.RK.0_5)


### Calculate SAND - RK

### Write function
bck_sand <- function(x) {(1/(1+exp(x[7])+exp(x[8])))*100}

beginCluster(7)
sable.preds.RK.0_5 <- clusterR(granulo.0_5.s, calc, args = list(bck_sand),
                               filename = "sable.preds.RK.0_5.tif", format = "GTiff",
                               na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(sable.preds.RK.0_5)

argile.preds.RK.0_5 <- raster("argile.preds.RK.0_5.tif")
limon.preds.RK.0_5 <- raster("limon.preds.RK.0_5.tif")
sable.preds.RK.0_5 <- raster("sable.preds.RK.0_5.tif")

### Add to the stack
granulo.0_5.s <- stack(granulo.0_5.s, argile.preds.RK.0_5,limon.preds.RK.0_5,sable.preds.RK.0_5)
names(granulo.0_5.s) <- c("argile.alr.cub.0_5" ,  "silt.alr.cub.0_5" ,    "argile.alr.ckR.0_5" ,
                          "limon.alr.ckR.0_5"  ,  "argile.alr.ckSTD.0_5", "limon.alr.ckSTD.0_5",
                           "argile.alr.preds.0_5", "limon.alr.preds.0_5" ,
                          "clay" , "silt" ,  "sand") ### change the name for the lm function


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
SMFC_0_5 <- clusterR(granulo.0_5.s, predict, args = list(lm_w20_ClSd),
                           filename = "SMFC_0_5.tif", format = "GTiff",
                           na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()

SMFC_0_5 <- raster("SMFC_0_5.tif")
plot(SMFC_0_5)

### Apply non-stratified PTF for w42
beginCluster(7)
SMPWP_0_5 <- clusterR(granulo.0_5.s, predict, args = list(lm_w42_ClSd),
                            filename = "SMPWP_0_5.tif", format = "GTiff",
                            na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
tac <- Sys.time()
tac-tic
SMPWP_0_5 <- raster("SMPWP_0_5.tif")
plot(SMPWP_0_5)


### Calculate difference
awc.0_5.s <- stack(SMFC_0_5, SMPWP_0_5)

difference <- function(x){x[1]-x[2]}

beginCluster(7)
ru_vol_0_5 <- clusterR(awc.0_5.s, calc, args = list(difference),
                      filename = "ru_vol_0_5.tif", format = "GTiff",
                      na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(ru_vol_0_5)

awc.0_5.s <- stack(awc.0_5.s, ru_vol_0_5)

### Multiply by layer depth
ru_mm <- function(x){50*x}

beginCluster(7)
ru_mm_0_5 <- clusterR(ru_vol_0_5, calc, args = list(ru_mm),
                       filename = "ru_mm_0_5.tif", format = "GTiff",
                       na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(ru_mm_0_5)

awc.0_5.s <- stack(awc.0_5.s, ru_mm_0_5)


### attach the coarse elements
awc.0_5.s <- stack(awc.0_5.s, coarse.qrF.mean.0_5)


### subtract the volume of coarse elements to the AWC (volume)
ru_coarse <- function(x){x[3]*(1-(x[5]/100))*50} ## in mm, taking into account coarse elements

beginCluster(7)
awc_mm_0_5 <- clusterR(awc.0_5.s, calc, args = list(ru_coarse),
                      filename = "awc_mm_0_5.tif", format = "GTiff",
                      na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
plot(awc_mm_0_5)

awc.0_5.s <- stack(awc.0_5.s, awc_mm_0_5)

### Save again
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
save.image("6.1.2.ru_0_5.RData")

#### end of script