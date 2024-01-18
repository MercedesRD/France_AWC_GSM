###########################################################################################################################################

### Objective: Calculate the sensitivities and variance terms for SMFC and SMPWP by GSM layer
### Author: Mercedes Roman Dobarco
### Date: 10/09/2018

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
library(compositions)
library(rgr)
library(rasterVis)
library(viridis)
library(scales)
library(foreach)

#######################################################################################################################
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

## Define the elements of the loop
depths <- c("0_5","5_15","15_30","30_60", "60_100", "100_200")
#gsm_thickness <- c(50,100,150,300,400,1000)

## Load the PTFs
load(paste0(HomeDir,"Input/ContinuousPTF.RData"))
rm(lm_w20_ClSdBd, lm_w20_ClSdBdSOC, lm_w20_ClSdSOC,
   lm_w20_sub_ClSd, lm_w20_sub_ClSdBd,lm_w20_sub_ClSdBdSOC, lm_w20_sub_ClSdSOC,
   lm_w20_top_ClSd, lm_w20_top_ClSdBd, lm_w20_top_ClSdSOC, lm_w20_top_ClSdBdSOC,
   lm_w25_ClSd, lm_w25_ClSdBd, lm_w25_ClSdSOC, lm_w25_ClSdBdSOC,
   lm_w25_sub_ClSd, lm_w25_sub_ClSdBd, lm_w25_sub_ClSdSOC, lm_w25_sub_ClSdBdSOC,
   lm_w25_top_ClSd, lm_w25_top_ClSdBd, lm_w25_top_ClSdSOC, lm_w25_top_ClSdBdSOC,
   lm_w42_ClSdBd, lm_w42_ClSdBdSOC, lm_w42_ClSdSOC,
   lm_w42_sub_ClSd, lm_w42_sub_ClSdBd,lm_w42_sub_ClSdBdSOC, lm_w42_sub_ClSdSOC,
   lm_w42_top_ClSd, lm_w42_top_ClSdBd, lm_w42_top_ClSdSOC, lm_w42_top_ClSdBdSOC)

# ptf2.0 <- lm_w20_ClSd
# ptf4.2 <- lm_w42_ClSd

### Read variance-covariance matrix of coefficients
lm_w20_ClSd_coeff_cov <-read.csv(paste0(HomeDir,"Input/lm_w20_ClSd_coeff_cov.csv"))
lm_w42_ClSd_coeff_cov <-read.csv(paste0(HomeDir,"Input/lm_w42_ClSd_coeff_cov.csv"))
lm_w20_ClSd_coeff_cov <- lm_w20_ClSd_coeff_cov[,-1]
lm_w42_ClSd_coeff_cov <- lm_w42_ClSd_coeff_cov[,-1]

### Load the SD of the input soil variables. From vector with file names and dire
clay.alr.sd.list <- c("6.1-CoK_Res_0_5/argile.alr.ckSTD.0_5.tif",
                      "6.2-CoK_Res_5_15/argile.alr.ckSTD.5_15.tif",
                      "6.3-CoK_Res_15_30/argile.alr.ckSTD.15_30.tif",
                      "6.4-CoK_Res_30_60/argile.alr.ckSTD.30_60.tif",
                      "6.5-CoK_Res_60_100/argile.alr.ckSTD.60_100.tif",
                      "6.6-CoK_Res_100_200/argile.alr.ckSTD.100_200.tif")
silt.alr.sd.list <- c("6.1-CoK_Res_0_5/limon.alr.ckSTD.0_5.tif",
                      "6.2-CoK_Res_5_15/limon.alr.ckSTD.5_15.tif",
                      "6.3-CoK_Res_15_30/limon.alr.ckSTD.15_30.tif",
                      "6.4-CoK_Res_30_60/limon.alr.ckSTD.30_60.tif",
                      "6.5-CoK_Res_60_100/limon.alr.ckSTD.60_100.tif",
                      "6.6-CoK_Res_100_200/limon.alr.ckSTD.100_200.tif")

### Define each sensitivity term, for un unespecified ptf

### Sensitivity clay

sens.clay  <- function(s,...){
    ### first term, clay.alr
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    dg_dx <-  ( 100 * exp(x) * (ptf[1][[1]][[2]] + (ptf[1][[1]][[2]]*exp(y)) - ptf[1][[1]][[3]]) ) / ((1+exp(x)+exp(y))^2)
    return(dg_dx)
}

sens.silt  <- function(s,...){
    ### first term, clay.alr
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    dg_dy <- ( -100 * exp(y) * ( (ptf[1][[1]][[2]]*exp(x)) + ptf[1][[1]][[3]] ) ) / ((1+exp(x)+exp(y))^2) 
    return(dg_dy)
}


### The sensitivity to the PTF to the coefficients are respectively 1 (intercept), clay content (beta clay), and sand content (beta sand), 
### so there is no need to recalculate them. They ahve been calculated and available at 
### "D:\romandobarco\AWC\AWC_GSM_Dec2017\Output\Recalculation_AWC"
### but we define them as functions for calculating the variance terms afterwards

sens.beta0 <- 1

sens.beta1 <- function(s,...){
        ### define clay.alr and silt.alr
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
        dg_dbeta1 <-  100 * ( exp(x) / (1+exp(x)+exp(y)))
    return(dg_dbeta1)
}

sens.beta2 <- function(s,...){
        ### define clay.alr and silt.alr
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    dg_dbeta2 <- 100 * ( 1 / (1+exp(x)+exp(y)))
    return(dg_dbeta2)
}


######################################################################################################################################

### Now the functions to calculate the variance terms

u.texture <- function(s,...){
    Cl.sd <- s[[3]]
    Sl.sd <- s[[4]]
    
    Cl.sens <- sens.clay(s)
    Sl.sens <- sens.silt(s)
    
    u.Cl <- (Cl.sd^2)*(Cl.sens^2)
    u.Sl <- (Sl.sd^2)*(Sl.sens^2)
    u.Cl.Sl <- 2*(Cl.sd*Sl.sd*rho)*(Cl.sens*Sl.sens)  ### where rho is the global correlation between clay-alr and silt-alr for that layer
    
    u.Texture <- u.Cl + u.Sl + u.Cl.Sl
    
    return(c(u.Texture, u.Cl, u.Sl, u.Cl.Sl))
}

u.ptf<- function(s,...){
    
    ### sensitivity beta0
    beta0.sens <- 1
    ### variance of beta0
    var.beta0 <- ClSd_coeff_cov[1,1]
    ubeta0 <- var.beta0 * (beta0.sens^2)
    
    ### sensitivity beta1
    beta1.sens <- sens.beta1(s)
    ### variance of beta1
    var.beta1 <- ClSd_coeff_cov[2,2]
    ubeta1 <- var.beta1 * (beta1.sens^2)
    
    ### sensitivity beta2
    beta2.sens <- sens.beta2(s)
    ### variance of beta2
    var.beta2 <- ClSd_coeff_cov[3,3]
    ubeta2 <- var.beta2 * (beta2.sens^2)
    
    ### covariance of beta0.beta1
    cov.beta0.beta1 <- ClSd_coeff_cov[1,2]
    ubeta0.beta1 <- 2 * cov.beta0.beta1 * beta0.sens * beta1.sens
    
    ### covariance of beta0.beta2
    cov.beta0.beta2 <- ClSd_coeff_cov[1,3]
    ubeta0.beta2 <- 2 * cov.beta0.beta2 * beta0.sens * beta2.sens
    
    ### covariance of beta1.beta2
    cov.beta1.beta2 <- ClSd_coeff_cov[2,3]
    ubeta1.beta2 <- 2 * cov.beta1.beta2 * beta1.sens * beta2.sens
    
    ### Add these 6
    u.PTF <- ubeta0 + ubeta1 + ubeta2 + ubeta0.beta1 + ubeta0.beta2 + ubeta1.beta2
    
    return(c(u.PTF, ubeta0, ubeta1, ubeta2, ubeta0.beta1, ubeta0.beta2, ubeta1.beta2 ))
} 

dir.create(paste0(HomeDir,"Clean_Output/12-Sensitivity_FC_PWP"))
setwd(paste0(HomeDir,"Clean_Output/12-Sensitivity_FC_PWP"))
save.image("functions.RData")

######################################################################################################################################

### Run loops to calculate these elements for the two water potentials and the six GSM layers


tic <- Sys.time()
### chopse the depth
for(h in 1:6){
    print(h) ## what depth am I working at?
    depth <- depths[[h]]
    
    ### Load the raster files for the chosen depth
    setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
    clay.alr <- raster(paste0("argile.alr.preds.",depth,".tif"))
    silt.alr <- raster(paste0("limon.alr.preds.",depth,".tif"))
    clay.alr.sd <- raster(paste0(HomeDir,"Clean_Output/",clay.alr.sd.list[[h]]))
    silt.alr.sd <- raster(paste0(HomeDir,"Clean_Output/",silt.alr.sd.list[[h]]))
    
    ### Create a working directory for this depth
    workDir <- paste0(HomeDir,"Clean_Output/12-Sensitivity_FC_PWP/var_FC_PWP.",depth)
    dir.create(workDir)
    setwd(workDir)
    awc.input <- stack(clay.alr, silt.alr, clay.alr.sd, silt.alr.sd)
    
    ### Calculate correlation between clay.alr and silt.alr
    clay.alr.v <- getValues(clay.alr)
    silt.alr.v <- getValues(silt.alr)
    
    ### calculate correlation
    cor.clay.limon.alr <- cor.test(x=clay.alr.v, y=silt.alr.v, na.rm=TRUE, method="pearson")
    rho <- cor.clay.limon.alr$estimate

    ### Load the functions and the PTFs
    load(paste0(HomeDir,"Clean_Output/12-Sensitivity_FC_PWP/functions.RData"))
    
    ####################################################################################################################
    
    ### 1. Soil moisture at field capacity
    ptf <- lm_w20_ClSd
    ClSd_coeff_cov <- lm_w20_ClSd_coeff_cov
    
    ### Sensitivity to clay.alr
    ff <- function(x) calc(x, sens.clay)
    beginCluster(7)
    sens.Clay.FC<- clusterR(awc.input, ff, export = list("sens.clay","ptf"),
                         filename = paste0("sens.clay.FC",".",depth,".tif"), format = "GTiff",
                         na.rm=T, inf.rm=T, progress = "text", overwrite = T)
    endCluster()
    
    sens.Clay.FC <- raster(paste0("sens.clay.FC",".",depth,".tif"))
    
    ### Sensitivity to silt.alr
    ff <- function(x) calc(x, sens.silt)
    beginCluster(7)
    sens.Silt.FC <- clusterR(awc.input, ff, export = list("sens.silt","ptf"),
                          filename = paste0("sens.silt.FC",".",depth,".tif"), format = "GTiff",
                          na.rm=T, inf.rm=T, progress = "text", overwrite = T)
    endCluster()
    #plot(sens.silt)
    sens.Silt.FC <- raster(paste0("sens.silt.FC",".",depth,".tif"))
    
    #####################################################################################################################
    
    ### Now the variance terms
    
    ff <- function(x) calc(x, u.texture)
    beginCluster(7)
    u.Texture <- clusterR( awc.input, ff, export = list("u.texture","sens.clay","sens.silt","ptf", "rho"),progress = "text")
    endCluster()
    gc()
    
    u.Texture <- stack(u.Texture)
    names(u.Texture) <- c("u.Texture", "u.Cl", "u.Sl", "u.Cl.Sl") 
    
    writeRaster(u.Texture, filename = c(paste0("u.Texture.FC.",depth,".tif"),
                                        paste0("u.clay.alr.FC.",depth,".tif"),
                                        paste0("u.silt.alr.FC.",depth,".tif"),
                                        paste0("u.clay.alr.silt.alr.FC.",depth,".tif")),
                format = "GTiff", bylayer=TRUE, 
                na.rm=T, inf.rm=T, overwrite = T)
    u.Texture <- stack(lapply(c(paste0("u.Texture.FC.",depth,".tif"), paste0("u.clay.alr.FC.",depth,".tif"),
                                paste0("u.silt.alr.FC.",depth,".tif"), paste0("u.clay.alr.silt.alr.FC.",depth,".tif")),
                              FUN= raster))
    
    
    #######################################################################################################################
    
    #### Variance terms associated to the PTFs coefficients
    
    ff <- function(x) calc(x, u.ptf)
    beginCluster(7)
    u.PTF.FC <- clusterR(awc.input, ff, export = list("u.ptf", "sens.beta1", "sens.beta2", "ClSd_coeff_cov"), progress = "text")
    endCluster()
    gc()
    
    u.PTF.FC <- stack(u.PTF.FC)     
    names(u.PTF.FC) <- c("u.PTF.FC", "ubeta0.FC", "ubeta1.FC", "ubeta2.FC", "ubeta0.beta1.FC", "ubeta0.beta2.FC", "ubeta1.beta2.FC")
    
    writeRaster(u.PTF.FC, filename = c(paste0("u.PTF.FC",".",depth,".tif"),
                                       paste0("u.beta0.FC",".",depth,".tif"),
                                       paste0("u.beta1.FC",".",depth,".tif"),
                                       paste0("u.beta2.FC",".",depth,".tif"),
                                       paste0("u.beta0beta1.FC",".",depth,".tif"),
                                       paste0("u.beta0beta2.FC",".",depth,".tif"),
                                       paste0("u.beta1beta2.FC",".",depth,".tif")),
                format = "GTiff", bylayer=TRUE, 
                na.rm=T, inf.rm=T, overwrite = T)
    
    u.PTF.FC <- stack(lapply(c(paste0("u.PTF.FC",".",depth,".tif"),
                               paste0("u.beta0.FC",".",depth,".tif"),
                               paste0("u.beta1.FC",".",depth,".tif"),
                               paste0("u.beta2.FC",".",depth,".tif"),
                               paste0("u.beta0beta1.FC",".",depth,".tif"),
                               paste0("u.beta0beta2.FC",".",depth,".tif"),
                               paste0("u.beta1beta2.FC",".",depth,".tif")),
                             FUN= raster))
    
    ####################################################################################################################
    
    ### 2. Soil moisture at permanent wilting point
    ptf <- lm_w42_ClSd
    ClSd_coeff_cov <- lm_w42_ClSd_coeff_cov
    
    ### Sensitivity to clay.alr
    ff <- function(x) calc(x, sens.clay)
    beginCluster(7)
    sens.Clay.PWP<- clusterR(awc.input, ff, export = list("sens.clay","ptf"),
                            filename = paste0("sens.clay.PWP",".",depth,".tif"), format = "GTiff",
                            na.rm=T, inf.rm=T, progress = "text", overwrite = T)
    endCluster()
    
    sens.Clay.PWP <- raster(paste0("sens.clay.PWP",".",depth,".tif"))
    
    ### Sensitivity to silt.alr
    ff <- function(x) calc(x, sens.silt)
    beginCluster(7)
    sens.Silt.PWP <- clusterR(awc.input, ff, export = list("sens.silt","ptf"),
                             filename = paste0("sens.silt.PWP",".",depth,".tif"), format = "GTiff",
                             na.rm=T, inf.rm=T, progress = "text", overwrite = T)
    endCluster()
    #plot(sens.silt)
    sens.Silt.PWP <- raster(paste0("sens.silt.PWP",".",depth,".tif"))
    
    #####################################################################################################################
    
    ### Now the variance terms
    
    ff <- function(x) calc(x, u.texture)
    beginCluster(7)
    u.Texture <- clusterR( awc.input, ff, export = list("u.texture","sens.clay","sens.silt","ptf", "rho"),progress = "text")
    endCluster()
    gc()
    
    u.Texture <- stack(u.Texture)
    names(u.Texture) <- c("u.Texture", "u.Cl", "u.Sl", "u.Cl.Sl") 
    
    writeRaster(u.Texture, filename = c(paste0("u.Texture.PWP.",depth,".tif"),
                                        paste0("u.clay.alr.PWP.",depth,".tif"),
                                        paste0("u.silt.alr.PWP.",depth,".tif"),
                                        paste0("u.clay.alr.silt.alr.PWP.",depth,".tif")),
                format = "GTiff", bylayer=TRUE, 
                na.rm=T, inf.rm=T, overwrite = T)
    u.Texture <- stack(lapply(c(paste0("u.Texture.PWP.",depth,".tif"), paste0("u.clay.alr.PWP.",depth,".tif"),
                                paste0("u.silt.alr.PWP.",depth,".tif"), paste0("u.clay.alr.silt.alr.PWP.",depth,".tif")),
                              FUN= raster))
    
    
    #######################################################################################################################
    
    #### Variance terms associated to the PTFs coefficients
    
    ff <- function(x) calc(x, u.ptf)
    beginCluster(7)
    u.PTF.PWP <- clusterR(awc.input, ff, export = list("u.ptf", "sens.beta1", "sens.beta2", "ClSd_coeff_cov"), progress = "text")
    endCluster()
    gc()
    
    u.PTF.PWP <- stack(u.PTF.PWP)     
    names(u.PTF.PWP) <- c("u.PTF.PWP", "ubeta0.PWP", "ubeta1.PWP", "ubeta2.PWP", "ubeta0.beta1.PWP", "ubeta0.beta2.PWP", "ubeta1.beta2.PWP")
    
    writeRaster(u.PTF.PWP, filename = c(paste0("u.PTF.PWP",".",depth,".tif"),
                                       paste0("u.beta0.PWP",".",depth,".tif"),
                                       paste0("u.beta1.PWP",".",depth,".tif"),
                                       paste0("u.beta2.PWP",".",depth,".tif"),
                                       paste0("u.beta0beta1.PWP",".",depth,".tif"),
                                       paste0("u.beta0beta2.PWP",".",depth,".tif"),
                                       paste0("u.beta1beta2.PWP",".",depth,".tif")),
                format = "GTiff", bylayer=TRUE, 
                na.rm=T, inf.rm=T, overwrite = T)
    
    u.PTF.PWP <- stack(lapply(c(paste0("u.PTF.PWP",".",depth,".tif"),
                               paste0("u.beta0.PWP",".",depth,".tif"),
                               paste0("u.beta1.PWP",".",depth,".tif"),
                               paste0("u.beta2.PWP",".",depth,".tif"),
                               paste0("u.beta0beta1.PWP",".",depth,".tif"),
                               paste0("u.beta0beta2.PWP",".",depth,".tif"),
                               paste0("u.beta1beta2.PWP",".",depth,".tif")),
                             FUN= raster))
    
    
    save.image(paste0("FC_PWP.var.decomp.",depth,".RData"))
    gc()
    rm(h, i)
}

tac <- Sys.time()
tac-tic


###### Now plot the variance terms with levelplot :)