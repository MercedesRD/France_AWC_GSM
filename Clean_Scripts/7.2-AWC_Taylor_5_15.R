#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 5-15 cm
###
###  calculating uncertainty of AWC predictions with Taylor analysis
###  Author: Mercedes Roman Dobarco
###  Date: 08/01/2018


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
library(Hmisc)
library(foreach)

# ### Load image from previous script, 6.2.2 --------------------

### set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
load("6.2.2.ru_5_15.RData")

### This will not work because the links to te raster files refer to the location in my computer, so I reload them again
setwd(paste0(HomeDir,"Clean_Output/5-Cubist_preds"))
silt.alr.cub.5_15 <- raster("silt.alr.cub.5_15.tif" )
argile.alr.cub.5_15 <- raster( "argile.alr.cub.5_15.tif")

setwd(paste0(HomeDir,"Clean_Output/6.2-CoK_Res_5_15"))
argile.alr.ckR.5_15 <- raster("argile.alr.ckR.5_15.tif")
argile.alr.ckSTD.5_15  <- raster("argile.alr.ckSTD.5_15.tif")
limon.alr.ckR.5_15 <- raster("limon.alr.ckR.5_15.tif")
limon.alr.ckSTD.5_15  <- raster("limon.alr.ckSTD.5_15.tif")

setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
coarse.qrF.mean.5_15 <- raster( "coarse_5_15.preds.mean.tif")
coarse.qrF.sd.5_15 <- raster( "coarse_5_15.preds.sd.tif")
coarse.qrF.05.5_15 <- raster( "coarse_5_15.preds.05p.tif")
coarse.qrF.95.5_15 <- raster( "coarse_5_15.preds.95p.tif")
coarse.5_15.s <- stack(coarse.qrF.mean.5_15, coarse.qrF.sd.5_15, coarse.qrF.05.5_15,coarse.qrF.95.5_15)

### Load raster files generated in previous script
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
argile.alr.preds.5_15 <- raster("argile.alr.preds.5_15.tif")
limon.alr.preds.5_15 <- raster("limon.alr.preds.5_15.tif")
argile.preds.RK.5_15 <- raster("argile.preds.RK.5_15.tif")
limon.preds.RK.5_15 <- raster("limon.preds.RK.5_15.tif")
sable.preds.RK.5_15 <- raster("sable.preds.RK.5_15.tif")

granulo.5_15.s <- stack(argile.alr.cub.5_15, silt.alr.cub.5_15,
                       argile.alr.ckR.5_15, limon.alr.ckR.5_15,
                       argile.alr.ckSTD.5_15, limon.alr.ckSTD.5_15,
                       argile.alr.preds.5_15,limon.alr.preds.5_15,
                       argile.preds.RK.5_15,limon.preds.RK.5_15,sable.preds.RK.5_15)

names(granulo.5_15.s) <- c("argile.alr.cub.5_15" ,  "silt.alr.cub.5_15" ,    "argile.alr.ckR.5_15" ,
                          "limon.alr.ckR.5_15"  ,  "argile.alr.ckSTD.5_15", "limon.alr.ckSTD.5_15",
                          "argile.alr.preds.5_15", "limon.alr.preds.5_15" ,
                          "clay" , "silt" ,  "sand") ### change the name for the lm function

SMFC_5_15 <- raster("SMFC_5_15.tif")
SMPWP_5_15 <- raster("SMPWP_5_15.tif")
ru_vol_5_15 <- raster("ru_vol_5_15.tif")
ru_mm_5_15 <- raster("ru_mm_5_15.tif")
awc_mm_5_15 <- raster("awc_mm_5_15.tif")
awc.5_15.s <- stack(SMFC_5_15, SMPWP_5_15,ru_vol_5_15,ru_mm_5_15,coarse.qrF.mean.5_15,awc_mm_5_15)

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


dir.create(paste0(HomeDir,"Clean_Output/7.2-AWC_Taylor_5_15"))
setwd(paste0(HomeDir,"Clean_Output/7.2-AWC_Taylor_5_15"))

# ### 1. Correlation between clay.alr and silt.alr -------------------------------

### Calculate correlation between clay.alr and silt.alr
tic <- Sys.time()
argile.alr.v <- getValues(argile.alr.preds.5_15)
silt.alr.v <- getValues(limon.alr.preds.5_15)

### calculate correlation
cor.clay.limon.alr.5_15 <- cor.test(x=argile.alr.v, y=silt.alr.v, na.rm=TRUE, method="pearson")
cor.clay.limon.alr.5_15$estimate
tac <- Sys.time()
tac-tic; rm(tic, tac)

##########################################################################################################################

### Set parameters for this layer
t_layer <- 100 ### thickness of layer in mm
### In this case, we are working for depth1
depths <- c("5_15", "5_15", "15_30", "30_60", "60_100", "100_200")
depth <- depths[[2]]


#######################################################################################################################################################

# ### 2. Calculate uncertainty of AWC --------

## my intput raster stack is granulo.5_15.s
names(granulo.5_15.s)
### The predictions for clay-alr (regression-kriging) is layer number 7
### the predictions for silt-alr (regression-kriging) is layer number 8
### The standard deviation for clay-alr residuals is layer number 5
### The standard deviation for silt-alr residuals is layer number 6
### Read variance-covariance matrix of coefficients
lm_w20_ClSd_coeff_cov <-read.csv(paste0(HomeDir,"Input/lm_w20_ClSd_coeff_cov.csv"))
lm_w42_ClSd_coeff_cov <-read.csv(paste0(HomeDir,"Input/lm_w42_ClSd_coeff_cov.csv"))
lm_w20_ClSd_coeff_cov <- lm_w20_ClSd_coeff_cov[,-1]
lm_w42_ClSd_coeff_cov <- lm_w42_ClSd_coeff_cov[,-1]

### From the word document where I explain how to calculate the variance by Taylor analysis:

### A Generic function, for any raster stack, and any PTF did not work in the first place... need to check again
### the function accepts any raster stack, but the specific PTF

# ptf[1][[1]][[2]] == lm_w20_ClSd$coefficients[[2]]

# ptf[1][[1]][[3]] == lm_w20_ClSd$coefficients[[3]]

### remember the variables in granulo.s
#  s[[1]] "argile.alr.cub.5_15" 
#  s[[2]] "silt.alr.cub.5_15" 
#  s[[3]] "argile.alr.ckR.5_15" 
#  s[[4]] "limon.alr.ckR.5_15"  
#  s[[5]] "argile.alr.ckSTD.5_15"
#  s[[6]] "limon.alr.ckSTD.5_15"
#  s[[7]] "argile.alr.preds.5_15"
#  s[[8]] "limon.alr.preds.5_15" 
#  s[[9]] "clay" 
#  s[[10]] "silt" 
#  s[[11]] "sand"

Uncertainty_input <- function(s,...){
    
    ### sum of terms in Taylor series
    
    ### first term, clay.alr
    ### x will be clay.alr
    x <- s[[7]]
    ### y will be silt.alr
    y <- s[[8]]
    ### z will be standard deviation of clay.alr
    z <- s[[5]]
    
    dg_dx <-  ( 100 * exp(x) * (ptf[1][[1]][[2]] + (ptf[1][[1]][[2]]*exp(y)) - ptf[1][[1]][[3]]) ) / ((1+exp(x)+exp(y))^2)
    u1 <- z*z*(dg_dx^2)
    
    ### Second term, silt.alr
    ### q will be standard deviation of silt.alr
    q <- s[[6]]
    
    dg_dy <- ( -100 * exp(y) * ( (ptf[1][[1]][[2]]*exp(x)) + ptf[1][[1]][[3]] ) ) / ((1+exp(x)+exp(y))^2) 
    u2 <- q*q*(dg_dy^2)
    
    ### Third term, silt.alr and clay.alr
    ### w will be covariance clay.alr.silt.alr, that I just calculated in the script 6.7
    # w <- s[[7]] 
    ### The results with the calculated covariance were really bad (negative variance, which should not be...)
    ### So, instead I use the correlation coefficient
    
    u3 <- 2*z*q*rho*dg_dx*dg_dy
    
    U_input <- u1 +u2 +u3
    
    return(U_input)
}


### Calculate uncertainty for the PTF at pF=2.0
### Write funciton for clusterR
ff <- function(x) calc(x, Uncertainty_input)
ptf <- lm_w20_ClSd
### assign the correlation extracted fomr the rasters with the regression-kriging predictions
rho <- cor.clay.limon.alr.5_15$estimate

beginCluster(7)
var_w20_input_5_15 <- clusterR(granulo.5_15.s, ff, export = list("Uncertainty_input","ptf", "rho"),
                          filename = paste0("varInput_w20_",depth,".tif"), format = "GTiff",
                          na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

var_w20_input_5_15 <- raster(paste0("varInput_w20_",depth,".tif"))


### Calculate uncertainty for the PTF at pF=4.2
ff <- function(x) calc(x, Uncertainty_input)
ptf <- lm_w42_ClSd

beginCluster(7)
var_w42_input_5_15  <- clusterR(granulo.5_15.s, ff, export = list("Uncertainty_input","ptf","rho"),
                          filename = paste0("varInput_w42_",depth,".tif"), format = "GTiff",
                          na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

var_w42_input_5_15 <- raster(paste0("varInput_w42_",depth,".tif"))

par(mfrow=c(3,2))
plot(var_w20_input_5_15)
plot(var_w42_input_5_15)

###############################################################################################################################################

# ### 3. Uncertainty of PTF -----------------------------------------------

### Start with a generic function for any PTF

# coeff.var.covar == lm_w20_ClSd_coeff_cov
# coeff.var.covar == lm_w42_ClSd_coeff_cov

#### Define my functions outside the loop

names(granulo.5_15.s)

### I can change these numbers in the function

Uncertainty_PTF <- function(s,...){
    
    ### define clay.alr and silt.alr
    ### x will be clay.alr
    x <- s[[7]]
    ### y will be silt.alr
    y <- s[[8]]
    
    ### sum of terms in Taylor series
    
    ### first term, intercept
    ### the uncertainty is constant
    ### The multiplication of 1 x variance(beta0)
    ### I extract it form the variance covariance matrix
    u_beta0 <- coeff.var.covar[1,1]
    
    ### second term, beta1 (clay)
    ### variance of beta1
    var.beta1 <- coeff.var.covar[2,2]
    dg_dbeta1 <-  100 * ( exp(x) / (1+exp(x)+exp(y)))
    u_beta1 <- var.beta1 *(dg_dbeta1^2)
    
    ###  third term, beta2 (sand)
    ### variance of beta2
    var.beta2 <- coeff.var.covar[3,3]
    dg_dbeta2 <- 100 * ( 1 / (1+exp(x)+exp(y)))
    u_beta2 <- var.beta2 *(dg_dbeta2^2)
    
    ### 4. interaction beta0 and beta1
    cov.beta0.beta1 <- coeff.var.covar[1,2]
    u_beta0.beta1 <- 2*cov.beta0.beta1*dg_dbeta1
    
    ### 5. interaction beta0 and beta2
    cov.beta0.beta2 <- coeff.var.covar[1,3]
    u_beta0.beta2 <- 2*cov.beta0.beta2*dg_dbeta2
    
    ### 6. interaction beta1 and beta2
    cov.beta1.beta2 <- coeff.var.covar[2,3]
    u_beta1.beta2 <- 2*cov.beta1.beta2*dg_dbeta2*dg_dbeta1
    
    ### sum of all
    U_ptf <- u_beta0 + u_beta1 + u_beta2 + u_beta0.beta1 + u_beta0.beta2 + u_beta1.beta2
    
    return(U_ptf)
}

### Calculate uncertainty for the PTF at pF=2.0
### Write funciton for clusterR
ff <- function(x) calc(x, Uncertainty_PTF)
coeff.var.covar <- lm_w20_ClSd_coeff_cov

beginCluster(7)
var_w20_ptf_5_15 <- clusterR(granulo.5_15.s, ff, export = list("Uncertainty_PTF","coeff.var.covar"),
                        filename = paste0("varPTF_w20_",depth,".tif"), format = "GTiff",
                        na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

var_w20_ptf_5_15 <- raster(paste0("varPTF_w20_",depth,".tif"))


### Calculate uncertainty for the PTF at pF=4.2
ff <- function(x) calc(x, Uncertainty_PTF)
coeff.var.covar <- lm_w42_ClSd_coeff_cov

beginCluster(7)
var_w42_ptf_5_15 <- clusterR(granulo.5_15.s, ff, export = list("Uncertainty_PTF","coeff.var.covar"),
                        filename = paste0("varPTF_w42_",depth,".tif"), format = "GTiff",
                        na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

gc()
var_w42_ptf_5_15 <- raster(paste0("varPTF_w42_",depth,".tif"))

plot(var_w20_ptf_5_15)
plot(var_w42_ptf_5_15)

####################################################################################################################

# ### 4. Calculate total variance BY LAYER --------------------------------

### Create raster stack with the variances
var.awc.s <- stack (var_w20_input_5_15, var_w42_input_5_15, var_w20_ptf_5_15, var_w42_ptf_5_15)

### write function
sum_unc_w20 <- function(x){x[1]+x[3]}
sum_unc_w42 <- function(x){x[2]+x[4]}


beginCluster(7)
total_var_w20.5_15 <- clusterR(var.awc.s, calc, args = list(sum_unc_w20),
                                 filename =paste0("Total_var_w20.",depth,".tif"), format = "GTiff",
                                 na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()

total_var_w20.5_15 <- raster(paste0("Total_var_w20.",depth,".tif"))


beginCluster(7)
total_var_w42.5_15 <- clusterR(var.awc.s, calc, args = list(sum_unc_w42),
                              filename =paste0("Total_var_w42.",depth,".tif"), format = "GTiff",
                              na.rm=T,inf.rm=T, progress = "text", overwrite = T)
endCluster()

total_var_w42.5_15 <- raster(paste0("Total_var_w42.",depth,".tif"))

plot(total_var_w20.5_15)
plot(total_var_w42.5_15)

# Add to the stack
var.awc.s <- stack(var.awc.s, total_var_w20.5_15, total_var_w42.5_15)

########################################################################################################################################################
### Clean the workspace a bit
rm(ru_coarse, difference, ff, ru_mm, tic, tac, rho, rho.awc, argile.alr.cub.5_15, silt.alr.cub.5_15,
   argile.alr.cub.5_15,silt.alr.cub.5_15,argile.alr.ckR.5_15,limon.alr.ckR.5_15,
   argile.alr.ckSTD.5_15 , limon.alr.ckSTD.5_15 , argile.alr.limon.alr.cov.5_15, argile.alr.preds.5_15,        
   limon.alr.preds.5_15, clay ,silt , sand,argile.preds.RK.5_15, limon.preds.RK.5_15, sable.preds.RK.5_15,
   var_w20_input_5_15, var_w42_input_5_15, var_w20_ptf_5_15, var_w42_ptf_5_15, var_awc_5_15, total_var_w42.5_15,
   total_var_w20.5_15,coarse.qrF.mean.5_15, coarse.qrF.sd.5_15, coarse.qrF.05.5_15, coarse.qrF.95.5_15,
   SMFC_5_15, SMPWP_5_15, ru_vol_5_15, ru_mm_5_15, silt.alr.v, argile.alr.v, SMFC_5_15.v, SMPWP_5_15.v,
   ptf, coeff.var.covar)

########################################################################################################################################################

# ### Variance of backtransformed clay, sand, silt ------------------------

### remember the variables in granulo.s
#  s[[1]] "argile.alr.cub.5_15" 
#  s[[2]] "silt.alr.cub.5_15" 
#  s[[3]] "argile.alr.ckR.5_15" 
#  s[[4]] "limon.alr.ckR.5_15"  
#  s[[5]] "argile.alr.ckSTD.5_15"
#  s[[6]] "limon.alr.ckSTD.5_15"
#  s[[7]] "argile.alr.preds.5_15"
#  s[[8]] "limon.alr.preds.5_15" 
#  s[[9]] "clay" 
#  s[[10]] "silt" 
#  s[[11]] "sand"

## my intput raster stack is granulo.5_15.s
names(granulo.5_15.s)

var_back_CLAY <- function(s,...){
    
    ### sum of terms in Taylor series
    
    ### first term, clay.alr
    ### x will be clay.alr
    x <- s[[7]]
    ### y will be silt.alr
    y <- s[[8]]
    ### z will be standard deviation of clay.alr
    z <- s[[5]]
    
    dg_dx <- (100*(exp(x)*(1+exp(y))))/((1+exp(x)+exp(y))^2) 
    u1 <- (z^2)*(dg_dx^2)
    
    ### Second term, silt.alr
    ### q will be standard deviation of silt.alr
    q <- s[[6]]
    
    dg_dy <-  (-exp(x)*exp(y)*100)/((1+exp(x)+exp(y))^2)
    u2 <- (q^2)*(dg_dy^2)
    
    ### Third term, silt.alr and clay.alr
    ### w will be covariance clay.alr.silt.alr 
    #w <- s[[7]] ## I prefer to use the correlation instead
    u3 <- 2*z*q*rho*dg_dx*dg_dy
    
    U_clay <- u1 +u2 +u3
    
    return(U_clay)
    
}

### Write funciton for clusterR
ff <- function(x) calc(x, var_back_CLAY)
rho <- cor.clay.limon.alr.5_15$estimate

beginCluster(7)
var_clay_5_15 <- clusterR(granulo.5_15.s, ff, export = list("var_back_CLAY","rho"),
                             filename = paste0("var_clay_",depth,".tif"), format = "GTiff",
                             na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

plot(var_clay_5_15)

########################################################################################################################################################################

### now SILT

## my intput raster stack is granulo.5_15.s
names(granulo.5_15.s)

var_back_SILT <- function(s,...){
    
    ### sum of terms in Taylor series
    
    ### first term, clay.alr
    ### x will be clay.alr
    x <- s[[7]]
    ### y will be silt.alr
    y <- s[[8]]
    ### z will be standard deviation of clay.alr
    z <- s[[5]]
    
    dg_dx <- (-exp(x)*exp(y)*100)/((1+exp(x)+exp(y))^2) 
    u1 <- (z^2)*(dg_dx^2)
    
    ### Second term, silt.alr
    ### q will be standard deviation of silt.alr
    q <- s[[6]]
    
    dg_dy <-  (100*(exp(y)*(1+exp(x))))/((1+exp(x)+exp(y))^2)
    u2 <- (q^2)*(dg_dy^2)
    
    ### Third term, silt.alr and clay.alr
    ### w will be covariance clay.alr.silt.alr 
    ## w <- s[[7]]
    u3 <- 2*z*q*rho*dg_dx*dg_dy
    
    U_silt <- u1 +u2 +u3
    
    return(U_silt)
    
}

### Write funciton for clusterR
ff <- function(x) calc(x, var_back_SILT)
rho <- cor.clay.limon.alr.5_15$estimate

beginCluster(7)
var_silt_5_15 <- clusterR(granulo.5_15.s, ff, export = list("var_back_SILT","rho"),
                         filename = paste0("var_silt_",depth,".tif"), format = "GTiff",
                         na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

plot(var_silt_5_15)


###################################################################################################################################################################


### Now SAND

## my intput raster stack is granulo.5_15.s
names(granulo.5_15.s)

var_back_SAND <- function(s,...){
    
    ### sum of terms in Taylor series
    
    ### first term, clay.alr
    ### x will be clay.alr
    x <- s[[7]]
    ### y will be silt.alr
    y <- s[[8]]
    ### z will be standard deviation of clay.alr
    z <- s[[5]]
    
    dg_dx <- (-100*exp(x))/((1+exp(x)+exp(y))^2) 
    u1 <- (z^2)*(dg_dx^2)
    
    ### Second term, silt.alr
    ### q will be standard deviation of silt.alr
    q <- s[[6]]
    
    dg_dy <- (-100*exp(y))/((1+exp(x)+exp(y))^2)
    u2 <- (q^2)*(dg_dy^2)
    
    ### Third term, silt.alr and clay.alr
    ### w will be covariance clay.alr.silt.alr 
    ## w <- s[[7]]
    u3 <- 2*z*q*rho*dg_dx*dg_dy
    
    U_sand <- u1 +u2 +u3
    
    return(U_sand)
    
}

### Write funciton for clusterR
ff <- function(x) calc(x, var_back_SAND)
rho <- cor.clay.limon.alr.5_15$estimate

beginCluster(7)
var_sand_5_15 <- clusterR(granulo.5_15.s, ff, export = list("var_back_SAND","rho"),
                         filename = paste0("var_sand_",depth,".tif"), format = "GTiff",
                         na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

rm(rho, ff)

plot(var_sand_5_15)

var_clay_5_15 <- raster(paste0("var_clay_",depth,".tif"))
var_silt_5_15 <- raster(paste0("var_silt_",depth,".tif"))
var_sand_5_15 <- raster(paste0("var_sand_",depth,".tif"))

### Save image
save.image("7.2-AWC_Taylor_5_15.RData")

###########################################################################################################################################################################

#### end of script