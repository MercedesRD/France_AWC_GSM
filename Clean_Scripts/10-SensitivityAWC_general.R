#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Decomposition of elementary AWC variance into different components:
###  - Coarse elements
###  - Texture
###  - PTF FC
###  - PTF PWP
###
###  Calculating uncertainty of AWC predictions with Taylor analysis
###  Examining sensitivity and variance of diffeernt sources
###  Author: Mercedes Roman Dobarco
###  Date: 07/08/2018

#############################################################################################################################################

### Calculate 16 terms of elementary AWC

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
dir.create(paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general"))

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

ptf2.0 <- lm_w20_ClSd
ptf4.2 <- lm_w42_ClSd

### Read variance-covariance matrix of coefficients
lm_w20_ClSd_coeff_cov <-read.csv(paste0(HomeDir,"Input/lm_w20_ClSd_coeff_cov.csv"))
lm_w42_ClSd_coeff_cov <-read.csv(paste0(HomeDir,"Input/lm_w42_ClSd_coeff_cov.csv"))
lm_w20_ClSd_coeff_cov <- lm_w20_ClSd_coeff_cov[,-1]
lm_w42_ClSd_coeff_cov <- lm_w42_ClSd_coeff_cov[,-1]

depths <- c("0_5","5_15","15_30","30_60", "60_100", "100_200")

### Define each sensibility in function of 2 PTFs, and a raster stack 

sens.rocks  <- function(s,...){
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    ### R will be rock content
    R <- s [[3]]
    dg_dR <- (-1/100) * (ptf2.0[1][[1]][[1]] - ptf4.2[1][[1]][[1]] + (((ptf2.0[1][[1]][[2]]*exp(x)) + ptf2.0[1][[1]][[3]] - (ptf4.2[1][[1]][[2]] * exp(x)) - ptf4.2[1][[1]][[3]])*100/(1+exp(x)+exp(y)))) 
    return(dg_dR)
}

sens.clay  <- function(s,...){
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    ### R will be rock content
    R <- s [[3]]
    
    dg_dx <- (1-(R/100)) * (100*exp(x)/((1+exp(x)+exp(y))^2)) * ( ((ptf2.0[1][[1]][[2]]-ptf4.2[1][[1]][[2]]) * (1+exp(y))) - ptf2.0[1][[1]][[3]] + ptf4.2[1][[1]][[3]] )
    
    return(dg_dx)
}

sens.silt  <- function(s,...){
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    ### R will be rock content
    R <- s [[3]]
    
    dg_dy <- (1-(R/100))*
        ((((ptf2.0[1][[1]][[2]] * exp(x)) + ptf2.0[1][[1]][[3]] - (ptf4.2[1][[1]][[2]] * exp(x)) - ptf4.2[1][[1]][[3]])*100*(-exp(y)))/
             ((1+exp(x)+exp(y))^2))
    return(dg_dy)
}

sens.beta0.FC  <- function(s,...){
 
    ### R will be rock content
    R <- s [[3]]
    dg_dBeta0 <- (1-(R/100))
    return(dg_dBeta0)
}

sens.beta1.FC  <- function(s,...){
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    ### R will be rock content
    R <- s [[3]]
    dg_dBeta1 <- (1-(R/100))* (100*exp(x)/(1+exp(x)+exp(y)))
    return(dg_dBeta1)
}

sens.beta2.FC  <- function(s,...){
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    ### R will be rock content
    R <- s [[3]]
    dg_dBeta2 <- (1-(R/100))* (100/(1+exp(x)+exp(y)))
    return(dg_dBeta2)
}

sens.beta0.PWP  <- function(s,...){
    
    ### R will be rock content
    R <- s [[3]]
    dg_dBeta0 <- - (1-(R/100))
    return(dg_dBeta0)
}

sens.beta1.PWP  <- function(s,...){
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    ### R will be rock content
    R <- s [[3]]
    dg_dBeta1 <- - (1-(R/100)) * (100*exp(x)/(1+exp(x)+exp(y)))
    return(dg_dBeta1)
}

sens.beta2.PWP  <- function(s,...){
    ### x will be clay.alr
    x <- s[[1]]
    ### y will be silt.alr
    y <- s[[2]]
    ### R will be rock content
    R <- s [[3]]
    dg_dBeta2 <- - (1-(R/100)) * (100/(1+exp(x)+exp(y)))
    return(dg_dBeta2)
}

##### Now define the different components of the variance (16 components)

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

tic <- Sys.time()
### chopse the depth
for(h in 1:6){
    #h <- 3
    print(h) ## what depth am I working at?
    depth <- depths[[h]]
    
    ### Load the raster files for the chosen depth
    setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
    clay.alr <- raster(paste0("argile.alr.preds.",depth,".tif"))
    silt.alr <- raster(paste0("limon.alr.preds.",depth,".tif"))
    
    setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
    rocks <- raster(paste0("coarse_",depth,".preds.mean.tif"))
    
    clay.alr.sd <- raster(paste0(HomeDir,"Clean_Output/",clay.alr.sd.list[[h]]))
    silt.alr.sd <- raster(paste0(HomeDir,"Clean_Output/",silt.alr.sd.list[[h]]))
    rocks.sd <- raster( paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps/coarse_",depth,".preds.sd.tif"))
    
   ### Create a working directory for this depth
    workDir <- paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general/varAWC.",depth)
    dir.create(workDir)
    setwd(workDir)
    awc.input <- stack(clay.alr, silt.alr, rocks, clay.alr.sd, silt.alr.sd, rocks.sd)

######################################################################################################################

## generate the 6(9) maps of sensitivities
    
    ### Define each sensibility in function of 2 PTFs, and a raster stack 
    
    sens.rocks  <- function(s,...){
        ### x will be clay.alr
        x <- s[[1]]
        ### y will be silt.alr
        y <- s[[2]]
        ### R will be rock content
        R <- s [[3]]
        dg_dR <- (-1/100) * (ptf2.0[1][[1]][[1]] - ptf4.2[1][[1]][[1]] + (((ptf2.0[1][[1]][[2]]*exp(x)) + ptf2.0[1][[1]][[3]] - (ptf4.2[1][[1]][[2]] * exp(x)) - ptf4.2[1][[1]][[3]])*100/(1+exp(x)+exp(y)))) 
        return(dg_dR)
    }
    
    sens.clay  <- function(s,...){
        ### x will be clay.alr
        x <- s[[1]]
        ### y will be silt.alr
        y <- s[[2]]
        ### R will be rock content
        R <- s [[3]]
        
        dg_dx <- (1-(R/100)) * (100*exp(x)/((1+exp(x)+exp(y))^2)) * ( ((ptf2.0[1][[1]][[2]]-ptf4.2[1][[1]][[2]]) * (1+exp(y))) - ptf2.0[1][[1]][[3]] + ptf4.2[1][[1]][[3]] )
        
        return(dg_dx)
    }
    
    sens.silt  <- function(s,...){
        ### x will be clay.alr
        x <- s[[1]]
        ### y will be silt.alr
        y <- s[[2]]
        ### R will be rock content
        R <- s [[3]]
        
        dg_dy <- (1-(R/100))*
            ((((ptf2.0[1][[1]][[2]] * exp(x)) + ptf2.0[1][[1]][[3]] - (ptf4.2[1][[1]][[2]] * exp(x)) - ptf4.2[1][[1]][[3]])*100*(-exp(y)))/
                 ((1+exp(x)+exp(y))^2))
        return(dg_dy)
    }
    
    sens.beta0.FC  <- function(s,...){
        
        ### R will be rock content
        R <- s [[3]]
        dg_dBeta0 <- (1-(R/100))
        return(dg_dBeta0)
    }
    
    sens.beta1.FC  <- function(s,...){
        ### x will be clay.alr
        x <- s[[1]]
        ### y will be silt.alr
        y <- s[[2]]
        ### R will be rock content
        R <- s [[3]]
        dg_dBeta1 <- (1-(R/100))* (100*exp(x)/(1+exp(x)+exp(y)))
        return(dg_dBeta1)
    }
    
    sens.beta2.FC  <- function(s,...){
        ### x will be clay.alr
        x <- s[[1]]
        ### y will be silt.alr
        y <- s[[2]]
        ### R will be rock content
        R <- s [[3]]
        dg_dBeta2 <- (1-(R/100))* (100/(1+exp(x)+exp(y)))
        return(dg_dBeta2)
    }
    
    sens.beta0.PWP  <- function(s,...){
        
        ### R will be rock content
        R <- s [[3]]
        dg_dBeta0 <- - ((1-(R/100)))
        return(dg_dBeta0)
    }
    
    sens.beta1.PWP  <- function(s,...){
        ### x will be clay.alr
        x <- s[[1]]
        ### y will be silt.alr
        y <- s[[2]]
        ### R will be rock content
        R <- s [[3]]
        dg_dBeta1 <- - ((1-(R/100)) * (100*exp(x)/(1+exp(x)+exp(y))))
        return(dg_dBeta1)
    }
    
    sens.beta2.PWP  <- function(s,...){
        ### x will be clay.alr
        x <- s[[1]]
        ### y will be silt.alr
        y <- s[[2]]
        ### R will be rock content
        R <- s [[3]]
        dg_dBeta2 <- - ((1-(R/100)) * (100/(1+exp(x)+exp(y))))
        return(dg_dBeta2)
    }
    
    

### Sensitivity to rocks
ff <- function(x) calc(x, sens.rocks)
beginCluster(7)
sens.Rocks <- clusterR(awc.input, ff, export = list("sens.rocks","ptf2.0", "ptf4.2"),
                             filename = paste0("sens.rocks",".",depth,".tif"), format = "GTiff",
                             na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(sens.rocks)
sens.Rocks <- raster(paste0("sens.rocks",".",depth,".tif"))

### Sensitivity to clay.alr
ff <- function(x) calc(x, sens.clay)
beginCluster(7)
sens.Clay<- clusterR(awc.input, ff, export = list("sens.clay","ptf2.0", "ptf4.2"),
                                 filename = paste0("sens.clay",".",depth,".tif"), format = "GTiff",
                                 na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(sens.clay)
sens.Clay <- raster(paste0("sens.clay",".",depth,".tif"))

### Sensitivity to silt.alr
ff <- function(x) calc(x, sens.silt)
beginCluster(7)
sens.Silt <- clusterR(awc.input, ff, export = list("sens.silt","ptf2.0", "ptf4.2"),
                            filename = paste0("sens.silt",".",depth,".tif"), format = "GTiff",
                            na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(sens.silt)
sens.Silt <- raster(paste0("sens.silt",".",depth,".tif"))

### Sensitivity to coefficients
ff <- function(x) calc(x, sens.beta0.FC)
beginCluster(7)
sens.Beta0 <- clusterR(awc.input, ff, export = list("sens.beta0.FC","ptf2.0", "ptf4.2"),
                            filename = paste0("sens.beta0.FC.",depth,".tif"), format = "GTiff",
                            na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(sens.beta0)
sens.Beta0 <- raster(paste0("sens.beta0.FC.",depth,".tif"))

### Sensitivity to clay coefficient
ff <- function(x) calc(x, sens.beta1.FC)
beginCluster(7)
sens.Beta1 <- clusterR(awc.input, ff, export = list("sens.beta1.FC","ptf2.0", "ptf4.2"),
                             filename = paste0("sens.beta1.FC.",depth,".tif"), format = "GTiff",
                             na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(sens.beta1)
sens.Beta1 <- raster(paste0("sens.beta1.FC.",depth,".tif"))


### Sensitivity to sand coefficient
ff <- function(x) calc(x, sens.beta2.FC)
beginCluster(7)
sens.Beta2 <- clusterR(awc.input, ff, export = list("sens.beta2.FC","ptf2.0", "ptf4.2"),
                             filename = paste0("sens.beta2.FC.",depth,".tif"), format = "GTiff",
                             na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(sens.beta2)
sens.Beta2 <- raster(paste0("sens.beta2.FC.",depth,".tif"))

#################################################################################################################################

#### Try to calculate the components of the variance directly

#### Rocks component

u.rocks <- function(s,...){
    R.sd <- s[[6]]
    R.sens <- sens.rocks(s)
    u.Rocks <- (R.sd^2)*(R.sens^2)
    return(u.Rocks)
}

### u rocks
ff <- function(x) calc(x, u.rocks)
beginCluster(7)
u.Rocks <- clusterR(awc.input, ff, export = list("u.rocks","sens.rocks","ptf2.0", "ptf4.2"),
                          filename = paste0("u.Rocks",".",depth,".tif"), format = "GTiff",
                          na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(u.Rocks)
u.Rocks <- raster(paste0("u.Rocks",".",depth,".tif"))

###################################################################################################################################

###################################################################################################################################

#### Texture component

u.texture <- function(s,...){
    Cl.sd <- s[[4]]
    Cl.sens <- sens.clay(s)
    
    Sl.sd <- s[[5]]
    Sl.sens <- sens.silt(s)
    
    u.Cl <- (Cl.sd^2)*(Cl.sens^2)
    u.Sl <- (Sl.sd^2)*(Sl.sens^2)
    u.Cl.Sl <- 2*(Cl.sd*Sl.sd*rho)*(Cl.sens*Sl.sens)
    
    u.Texture <- u.Cl + u.Sl + u.Cl.Sl
    
    return(c(u.Texture, u.Cl, u.Sl, u.Cl.Sl))
}

### Calculate correlation between clay.alr and silt.alr
clay.alr.v <- getValues(clay.alr)
silt.alr.v <- getValues(silt.alr)

### calculate correlation
cor.clay.limon.alr <- cor.test(x=clay.alr.v, y=silt.alr.v, na.rm=TRUE, method="pearson")
rho <- cor.clay.limon.alr$estimate


ff <- function(x) calc(x, u.texture)
beginCluster(7)
u.Texture <- clusterR( awc.input, ff, export = list("u.texture","sens.clay","sens.silt","ptf2.0", "ptf4.2", "rho"),progress = "text")
endCluster()
gc()

u.Texture <- stack(u.Texture)
names(u.Texture) <- c("u.Texture", "u.Cl", "u.Sl", "u.Cl.Sl") 

writeRaster(u.Texture, filename = c(paste0("u.Texture.",depth,".tif"),
                                    paste0("u.clay.alr.",depth,".tif"),
                                    paste0("u.silt.alr.",depth,".tif"),
                                    paste0("u.clay.alr.silt.alr.",depth,".tif")),
            format = "GTiff", bylayer=TRUE, 
            na.rm=T, inf.rm=T, overwrite = T)
u.Texture <- stack(lapply(c(paste0("u.Texture.",depth,".tif"), paste0("u.clay.alr.",depth,".tif"), paste0("u.silt.alr.",depth,".tif"), paste0("u.clay.alr.silt.alr.",depth,".tif")),
                            FUN= raster))
#plot(u.Texture)
 
##################################################################################################################################

##################################################################################################################################


### PTF for SMFC

u.ptf.FC<- function(s,...){
    
    ### sensitivity beta0
    beta0.sens <- sens.beta0.FC(s)
    ### variance of beta0
    var.beta0 <- lm_w20_ClSd_coeff_cov[1,1]
    ubeta0 <- var.beta0 * (beta0.sens^2)
    
    ### sensitivity beta1
    beta1.sens <- sens.beta1.FC(s)
    ### variance of beta1
    var.beta1 <- lm_w20_ClSd_coeff_cov[2,2]
    ubeta1 <- var.beta1 * (beta1.sens^2)
    
    ### sensitivity beta2
    beta2.sens <- sens.beta2.FC(s)
    ### variance of beta2
    var.beta2 <- lm_w20_ClSd_coeff_cov[3,3]
    ubeta2 <- var.beta2 * (beta2.sens^2)
    
    ### covariance of beta0.beta1
    cov.beta0.beta1 <- lm_w20_ClSd_coeff_cov[1,2]
    ubeta0.beta1 <- 2 * cov.beta0.beta1 * beta0.sens * beta1.sens
    
    ### covariance of beta0.beta2
    cov.beta0.beta2 <- lm_w20_ClSd_coeff_cov[1,3]
    ubeta0.beta2 <- 2 * cov.beta0.beta2 * beta0.sens * beta2.sens
  
    ### covariance of beta1.beta2
    cov.beta1.beta2 <- lm_w20_ClSd_coeff_cov[2,3]
    ubeta1.beta2 <- 2 * cov.beta1.beta2 * beta1.sens * beta2.sens
    
    ### Add these 6
    u.PTF.FC <- ubeta0 + ubeta1 + ubeta2 + ubeta0.beta1 + ubeta0.beta2 + ubeta1.beta2
    
    return(c(u.PTF.FC, ubeta0, ubeta1, ubeta2, ubeta0.beta1, ubeta0.beta2, ubeta1.beta2 ))
} 

ff <- function(x) calc(x, u.ptf.FC)
beginCluster(7)
u.PTF.FC <- clusterR(awc.input, ff, export = list("u.ptf.FC", "sens.beta0.FC", "sens.beta1.FC", "sens.beta2.FC", "lm_w20_ClSd_coeff_cov"), progress = "text")
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

#plot(u.PTF.FC)

##################################################################################################################################

##################################################################################################################################

### PTF for SMPWP

u.ptf.PWP<- function(s,...){
    
    ### sensitivity beta0
    beta0.sens <- sens.beta0.PWP(s)
    ### variance of beta0
    var.beta0 <- lm_w42_ClSd_coeff_cov[1,1]
    ubeta0 <- var.beta0 * (beta0.sens^2)
    
    ### sensitivity beta1
    beta1.sens <- sens.beta1.PWP(s)
    ### variance of beta1
    var.beta1 <- lm_w42_ClSd_coeff_cov[2,2]
    ubeta1 <- var.beta1 * (beta1.sens^2)
    
    ### sensitivity beta2
    beta2.sens <- sens.beta2.PWP(s)
    ### variance of beta2
    var.beta2 <- lm_w42_ClSd_coeff_cov[3,3]
    ubeta2 <- var.beta2 * (beta2.sens^2)
    
    ### covariance of beta0.beta1
    cov.beta0.beta1 <- lm_w42_ClSd_coeff_cov[1,2]
    ubeta0.beta1 <- 2 * cov.beta0.beta1 * beta0.sens * beta1.sens
    
    ### covariance of beta0.beta2
    cov.beta0.beta2 <- lm_w42_ClSd_coeff_cov[1,3]
    ubeta0.beta2 <- 2 * cov.beta0.beta2 * beta0.sens * beta2.sens
    
    ### covariance of beta1.beta2
    cov.beta1.beta2 <- lm_w42_ClSd_coeff_cov[2,3]
    ubeta1.beta2 <- 2 * cov.beta1.beta2 * beta1.sens * beta2.sens
    
    ### Add these 6
    u.PTF.PWP <- ubeta0 + ubeta1 + ubeta2 + ubeta0.beta1 + ubeta0.beta2 + ubeta1.beta2
    
    return(c(u.PTF.PWP, ubeta0, ubeta1, ubeta2, ubeta0.beta1, ubeta0.beta2, ubeta1.beta2 ))
} 

ff <- function(x) calc(x, u.ptf.PWP)
beginCluster(7)
u.PTF.PWP <- clusterR(awc.input, ff, export = list("u.ptf.PWP", "sens.beta0.PWP", "sens.beta1.PWP", "sens.beta2.PWP", "lm_w42_ClSd_coeff_cov"), progress = "text")
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

#plot(u.PTF.PWP)

##################################################################################################################################

##################################################################################################################################

### Add the 4 sources of uncertainty
# par(mfrow=c(2,2))
# plot(u.Rocks)
# plot(u.Texture[[1]])
# plot(u.PTF.FC[[1]])
# plot(u.PTF.PWP[[1]])

u.PTFs <- overlay(u.PTF.FC[[1]], u.PTF.PWP[[1]], fun=sum,
                  filename=paste0("u.PTFs",".",depth,".tif"),
                  format="GTiff", overwrite=TRUE)

u.rast <- stack(u.Rocks, u.Texture[[1]], u.PTFs)
#plot(u.rast)

### Calculate total ELEMENTARY VARIANCE OF AWC (cm3/cm3)
sum_all <- function(x){x[1]+x[2]+x[3]}

beginCluster(7)
var.E.AWC.dec <- clusterR(u.rast, calc, args = list(sum_all),
                                filename = paste0("var.E.AWC.dec.",depth,".tif"), format = "GTiff",
                                na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()
#plot(var.E.AWC.dec)

# ### Comparison with the elementary AWC variance.....
# sd_awc_R_vol_15_30 <- raster("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/7.3.0-AWC_Taylor_15_30/sd_awc_R_vol_15_30.tif")
# var.E.AWC.1 <- calc(sd_awc_R_vol_15_30, function(x) {x^2})
# 
# par(mfrow=c(1,4))
# plot(var.E.AWC.1)
# plot(var.E.AWC.dec)
# plot((var.E.AWC.1-var.E.AWC.dec)/var.E.AWC.dec*100)

### Calculate the contribution of each source to toal variance (%)

contrib.sources <- list()

for(i in 1:3){
    contrib.sources.i <- overlay(u.rast[[i]], var.E.AWC.dec, fun = function(x,y){x/y*100} ,
                                    filename=paste0("contrib.",names(u.rast[[i]]),".tif"),
                                    format="GTiff", overwrite=TRUE)
    contrib.sources[[i]] <- contrib.sources.i
}

contrib.sources <- stack(contrib.sources)
#plot(contrib.sources)

png(filename=paste0("Contrib.sources.", depth, ".png"), width = 3000, height = 3000, res=600)
par(mfrow=c(1,1))
plotRGB(contrib.sources, r=1,g=2,b=3, scale=100, stretch="hist")
legend(x =80000, y =6508465 ,c("U rocks", "U texture", "U PTFs"),
       col=c("red","green","blue"),
       pch=rep(15,3),
       title= expression(atop("Contribution (%)",
                              "to var(AWC)")),bty="n", cex=1)
dev.off()

save.image(paste0("AWC.var.decomp.",depth,".RData"))
gc()
rm(h, i)

}

tac <- Sys.time()
tac-tic

### end of the script