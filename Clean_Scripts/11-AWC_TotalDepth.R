#### Available water capacity
#### Uncertainty calculated with Taylor analysis and Monte-Carlo simulations (for depth)

### Objectives:
### 1. Calculate AWC to the total depth (maximum 2 m) and rock fragments
### 2. Calculate associated variance
### the step 2 was repeated 13/08/2018 with the new GSM layers of elementary AWC (volumetric, without multiplying by layer thickness) variance.

### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(geoR)
library(lattice)
library(spatstat)
library(automap)
library(ggplot2)
library(Hmisc)
library(plyr)
library(soiltexture)
library(compositions)
library(rgr)
library(foreach)
library(doParallel)
library(ithir)
library(raster)
library(snow)
library(doSNOW)
library(rasterVis)
library(viridis)
library(scales)
library(foreach)
library(parallel)
library(doParallel)



# ### Load raster files ---------------------------------------------------

HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
ru_vol.0_5 <- raster("ru_vol_0_5.tif")
ru_vol.5_15 <- raster("ru_vol_5_15.tif")
ru_vol.15_30 <- raster("ru_vol_15_30.tif")
ru_vol.30_60 <- raster("ru_vol_30_60.tif")
ru_vol.60_100 <- raster("ru_vol_60_100.tif")
ru_vol.100_200 <- raster("ru_vol_100_200.tif")

setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
coarse.0_5 <- raster("coarse_0_5.preds.mean.tif")
coarse.5_15 <- raster("coarse_5_15.preds.mean.tif")
coarse.15_30 <- raster("coarse_15_30.preds.mean.tif")
coarse.30_60 <- raster("coarse_30_60.preds.mean.tif")
coarse.60_100 <- raster("coarse_60_100.preds.mean.tif")
coarse.100_200 <- raster("coarse_100_200.preds.mean.tif")


# ### 1 Take into account the depth of the soil profile and rock f --------

### make raster stack with the 4 layers (var(awc_vol))
AWC_vol.s <- stack(ru_vol.0_5,coarse.0_5,ru_vol.5_15,coarse.5_15,ru_vol.15_30,coarse.15_30,
                   ru_vol.30_60,coarse.30_60,ru_vol.60_100,coarse.60_100,ru_vol.100_200,coarse.100_200)

### Soil depth
### next steps were done previously in another script
### The GBM depth prediction
# soil_depth <- raster("D:/romandobarco/AWC/AWC_v0/GSM_90m/Depth/GBM_RMQS/GBM_QuantileCorrection/predpsol_GBM_RMQS_qc.tif")
# plot(soil_depth)
# ###Crop to the same extent
# soil_depth <- crop(soil_depth,AWC_vol.s[[1]])
# ### transform to mm
# x10 <- function(x) {x*10}
# beginCluster(8)
# soil_depth_mm <- clusterR(soil_depth, calc, args = list(x10), 
#                           filename = "soil_depth_mm.tif", format = "GTiff",
#                           na.rm=T,inf.rm=T, progress = "text",overwrite = T)
# endCluster()

### Load soil depth in mm
soil_depth_mm <- raster(paste0(HomeDir,"Clean_Output/Input/soil_depth_mm.tif"))

### Add to the stack
AWC_vol.s <- stack(AWC_vol.s, soil_depth_mm)
names(AWC_vol.s)

AWC_depth_R <- function(s) {  ### s is a raster stack with AWC by depth layer, and soil depth
    ### First, let's create an empty vector and fill the values
    totalAWC <- rep(NA, length(s[1]))
    totalAWC[is.na(s[1]) & is.na(s[3]) & is.na(s[5]) & is.na(s[7]) & is.na(s[9]) & is.na(s[11]) |
                 is.na(s[2]) & is.na(s[4]) & is.na(s[6]) & is.na(s[8]) & is.na(s[10]) & is.na(s[12])|
                 is.na(s[13])] <- NA ### If there is a missing value for soil depth, or the whole AWC profile, assign NA
    totalAWC[s[13] <= 50]                   <- s[13]*s[1]*(1-(s[2]/100))                                            
    totalAWC[s[13] > 50  & s[13]  <= 150]   <- sum((50*s[1]*(1-(s[2]/100))),((s[13]-50)*s[3]*(1-(s[4]/100))), na.rm=T)
    totalAWC[s[13] > 150 & s[13]  <= 300]   <- sum((50*s[1]*(1-(s[2]/100))),(100*s[3]*(1-(s[4]/100))), ((s[13]-150)*s[5]*(1-(s[6]/100))), na.rm=TRUE)
    totalAWC[s[13] > 300 & s[13]  <= 600]   <- sum((50*s[1]*(1-(s[2]/100))),(100*s[3]*(1-(s[4]/100))), (150*s[5]*(1-(s[6]/100))), ((s[13]-300)*s[7]*(1-(s[8]/100))), na.rm=TRUE)
    totalAWC[s[13] > 600 & s[13]  <= 1000]  <- sum((50*s[1]*(1-(s[2]/100))),(100*s[3]*(1-(s[4]/100))), (150*s[5]*(1-(s[6]/100))), (300*s[7]*(1-(s[8]/100))), ((s[13]-600)*s[9]*(1-(s[10]/100))), na.rm=TRUE)
    totalAWC[s[13] > 1000 & s[13] <= 2000]  <- sum((50*s[1]*(1-(s[2]/100))),(100*s[3]*(1-(s[4]/100))), (150*s[5]*(1-(s[6]/100))), (300*s[7]*(1-(s[8]/100))), (400*s[9]*(1-(s[10]/100))),
                                                   ((s[13]-1000)*s[11]*(1-(s[12]/100))), na.rm=TRUE)
    totalAWC[s[13] > 2000]  <- sum((50*s[1]*(1-(s[2]/100))),(100*s[3]*(1-(s[4]/100))), (150*s[5]*(1-(s[6]/100))), (300*s[7]*(1-(s[8]/100))), (400*s[9]*(1-(s[10]/100))),
                                                   (1000*s[11]*(1-(s[12]/100))), na.rm=TRUE)
    totalAWC
}

### As it is, we could apply it like
dir.create(paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth"))
setwd(paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth"))

ff <- function(x) calc(x, AWC_depth_R)
beginCluster(7)
AWCr_mm_0_200 <- clusterR(AWC_vol.s, ff, export = list("AWC_depth_R"),
                               filename = "AWCr_mm_0_200.tif", format = "GTiff",
                               na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

par(mfrow=c(1,1))
plot(AWCr_mm_0_200)
plot(AWCr_mm_0_200, breaks=c(0,20,50,80,120,160,200,250,300,360,430),
     col=rev(viridis_pal(option="D")(10)))

############################################################################################################################################################

# ### 2. Now calculate the variance for the whole profile, without considering variation in depth --------------------

###################################################################################################################################################

### Load the AWC (mm) and calculate elementary AWC by layer (cm3/cm3) (after subtraction of coarse elements, but without multiplying by layer thickness)
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
awc_mm_0_5 <- raster("awc_mm_0_5.tif")
awc_mm_5_15 <- raster("awc_mm_5_15.tif")
awc_mm_15_30 <- raster("awc_mm_15_30.tif")
awc_mm_30_60 <- raster("awc_mm_30_60.tif")
awc_mm_60_100 <- raster("awc_mm_60_100.tif")
awc_mm_100_200 <- raster("awc_mm_100_200.tif")

rasters.awc.list <- list(awc_mm_0_5,awc_mm_5_15, awc_mm_15_30,
                         awc_mm_30_60, awc_mm_60_100, awc_mm_100_200 )
t <- c(50, 100, 150, 300, 400, 1000)
depths <- c("0_5", "5_15", "15_30", "30_60", "60_100","100_200")

setwd(paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth"))
for(r in 1:length(rasters.awc.list)){
    
    print(r)
    t_layer <- t[r]
    awc_vol.f <- function(x) {x/t_layer}
    
    beginCluster(8)
    clusterR(rasters.awc.list[[r]], calc, args = list(awc_vol.f), export=list("t_layer"),
             filename = paste0("awc_R_vol_",depths[r],".tif"), format = "GTiff",
             na.rm=T,inf.rm=T, progress = "text",overwrite = T)
    endCluster()
    
}

### Calculate the correlation between elementary AWC
awcR.files <- c("awc_R_vol_0_5.tif","awc_R_vol_5_15.tif","awc_R_vol_15_30.tif",
                   "awc_R_vol_30_60.tif","awc_R_vol_60_100.tif","awc_R_vol_100_200.tif")
awcR.s <- list()
for(i in 1:6){
    awcR.s[[i]] <- raster(awcR.files[i])
    
}
awcR.s <- stack(awcR.s) 

# ### add the depth
 # awcR.s <- stack(awcR.s, soil_depth_mm)

### calculate correlations between AWC by layer
cor.AWC1.2 <- cor.test(x=getValues(awcR.s[[1]]), y=getValues(awcR.s[[2]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC1.3 <- cor.test(x=getValues(awcR.s[[1]]), y=getValues(awcR.s[[3]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC1.4 <- cor.test(x=getValues(awcR.s[[1]]), y=getValues(awcR.s[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC1.5 <- cor.test(x=getValues(awcR.s[[1]]), y=getValues(awcR.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC1.6 <- cor.test(x=getValues(awcR.s[[1]]), y=getValues(awcR.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC2.3 <- cor.test(x=getValues(awcR.s[[2]]), y=getValues(awcR.s[[3]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC2.4 <- cor.test(x=getValues(awcR.s[[2]]), y=getValues(awcR.s[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC2.5 <- cor.test(x=getValues(awcR.s[[2]]), y=getValues(awcR.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC2.6 <- cor.test(x=getValues(awcR.s[[2]]), y=getValues(awcR.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC3.4 <- cor.test(x=getValues(awcR.s[[3]]), y=getValues(awcR.s[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC3.5 <- cor.test(x=getValues(awcR.s[[3]]), y=getValues(awcR.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC3.6 <- cor.test(x=getValues(awcR.s[[3]]), y=getValues(awcR.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC4.5 <- cor.test(x=getValues(awcR.s[[4]]), y=getValues(awcR.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC4.6 <- cor.test(x=getValues(awcR.s[[4]]), y=getValues(awcR.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.AWC5.6 <- cor.test(x=getValues(awcR.s[[5]]), y=getValues(awcR.s[[6]]), na.rm=TRUE, method="pearson")$estimate

save.image("11-AWC_TotalDepth.RData")
### clean the space
rm(i, t_layer, t2, t, rasters.awc.list,r, awcR.files, awc_mm_0_5, awc_mm_5_15, awc_mm_30_60, awc_mm_60_100, awc_mm_100_200,
   awc_mm_15_30, var_awc_R_mm_0_5, var_awc_R_mm_100_200, var_awc_R_mm_5_15, var_awc_R_mm_15_30, var_awc_R_mm_30_60, var_awc_R_mm_60_100)
rm(var_awc_R_mm_0_5,var_awc_R_mm_5_15, var_awc_R_mm_15_30, var_awc_R_mm_30_60, var_awc_R_mm_60_100, var_awc_R_mm_100_200, awc_vol.f)
##############################################################################################################################################################

### Load the variances of the elementary AWc into a raster stack

### New files 13/08/2018
awcR.files <-c(paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general/varAWC.0_5/var.E.AWC.dec.0_5.tif"),
               paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general/varAWC.5_15/var.E.AWC.dec.5_15.tif"),
               paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general/varAWC.15_30/var.E.AWC.dec.15_30.tif"),
               paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general/varAWC.30_60/var.E.AWC.dec.30_60.tif"),
               paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general/varAWC.60_100/var.E.AWC.dec.60_100.tif"),
               paste0(HomeDir,"Clean_Output/10-SensitivityAWC_general/varAWC.100_200/var.E.AWC.dec.100_200.tif"))

# awcR.files <- c("var_awc_R_vol_0_5.tif","var_awc_R_vol_5_15.tif","var_awc_R_vol_15_30.tif",
#                 "var_awc_R_vol_30_60.tif","var_awc_R_vol_60_100.tif","var_awc_R_vol_100_200.tif")
awcRvar.s <- list()
for(i in 1:6){
    awcRvar.s[[i]] <- raster(awcR.files[i])
    
}
awcRvar.s <- stack(awcRvar.s) 
rm(awcR.files, i)

### Add the soil depth
awcRvar.s <- stack(awcRvar.s,soil_depth_mm)
plot(awcRvar.s)
names(awcRvar.s) <-  c("var.E.AWC.0_5","var.E.AWC.5_15",  "var.E.AWC.15_30" ,  "var.E.AWC.30_60" ,  "var.E.AWC.60_100",
                       "var.E.AWC.100_200" ,"soil_depth_mm" )

AWCr_D_Var<- function(s) {  ### s is a raster stack 
    ### First, let's create an empty vector and fill the values
    varAWC <- rep(NA, length(s[[1]]))
    
    varAWC[is.na(s[[1]]) & is.na(s[[2]]) & is.na(s[[3]]) & is.na(s[[4]]) & is.na(s[[5]]) & is.na(s[[6]]) |
               is.na(s[[7]])] <- NA ### If there is a missing value for soil depth, or the whole AWC profile, assign NA
    
    varAWC[s[[7]] <= 50]                   <- (s[[7]]^2)*s[[1]]                                           
    
    varAWC[s[[7]] > 50  & s[[7]]  <= 150]   <- sum((2500*s[[1]]),(((s[[7]]-50)^2)*s[[2]]),(2*cor.AWC1.2*sqrt(s[[1]])*sqrt(s[[2]])*50*(s[[7]]-50)), na.rm=T)
    
    varAWC[s[[7]] > 150 & s[[7]]  <= 300]   <- sum((2500*s[[1]]),(10000*s[[2]]),(((s[[7]]-150)^2)*s[[3]]),
                                               (2*cor.AWC1.2*sqrt(s[[1]])*sqrt(s[[2]])*50*100),
                                               (2*cor.AWC1.3*sqrt(s[[1]])*sqrt(s[[3]])*50*(s[[7]]-150)),
                                               (2*cor.AWC2.3*sqrt(s[[2]])*sqrt(s[[3]])*100*(s[[7]]-150)), na.rm=TRUE)
    
    varAWC[s[[7]] > 300 & s[[7]]  <= 600]   <- sum((2500*s[[1]]),(10000*s[[2]]),(22500*s[[3]]),(((s[[7]]-300)^2)*s[[4]]),
                                               (2*cor.AWC1.2*sqrt(s[[1]])*sqrt(s[[2]])*50*100),
                                               (2*cor.AWC1.3*sqrt(s[[1]])*sqrt(s[[3]])*50*150),
                                               (2*cor.AWC2.3*sqrt(s[[2]])*sqrt(s[[3]])*100*150),
                                               (2*cor.AWC1.4*sqrt(s[[1]])*sqrt(s[[4]])*50*(s[[7]]-300)),
                                               (2*cor.AWC2.4*sqrt(s[[2]])*sqrt(s[[4]])*100*(s[[7]]-300)),
                                               (2*cor.AWC3.4*sqrt(s[[3]])*sqrt(s[[4]])*150*(s[[7]]-300)), na.rm=TRUE)
    
    varAWC[s[[7]] > 600 & s[[7]]  <= 1000]  <- sum((2500*s[[1]]),(10000*s[[2]]),(22500*s[[3]]),(90000*s[[4]]),(((s[[7]]-600)^2)*s[[5]]),
                                               (2*cor.AWC1.2*sqrt(s[[1]])*sqrt(s[[2]])*50*100),
                                               (2*cor.AWC1.3*sqrt(s[[1]])*sqrt(s[[3]])*50*150),
                                               (2*cor.AWC2.3*sqrt(s[[2]])*sqrt(s[[3]])*100*150),
                                               (2*cor.AWC1.4*sqrt(s[[1]])*sqrt(s[[4]])*50*300),
                                               (2*cor.AWC2.4*sqrt(s[[2]])*sqrt(s[[4]])*100*300),
                                               (2*cor.AWC3.4*sqrt(s[[3]])*sqrt(s[[4]])*150*300),
                                               (2*cor.AWC1.5*sqrt(s[[1]])*sqrt(s[[5]])*50*(s[[7]]-600)),
                                               (2*cor.AWC2.5*sqrt(s[[2]])*sqrt(s[[5]])*100*(s[[7]]-600)),
                                               (2*cor.AWC3.5*sqrt(s[[3]])*sqrt(s[[5]])*150*(s[[7]]-600)),
                                               (2*cor.AWC4.5*sqrt(s[[4]])*sqrt(s[[5]])*300*(s[[7]]-600)), na.rm=TRUE)
    
    varAWC[s[[7]] > 1000 & s[[7]] <= 2000]  <- sum((2500*s[[1]]),(10000*s[[2]]),(22500*s[[3]]),(90000*s[[4]]),(160000*s[[5]]),(((s[[7]]-1000)^2)*s[[6]]),
                                               (2*cor.AWC1.2*sqrt(s[[1]])*sqrt(s[[2]])*50*100),
                                               (2*cor.AWC1.3*sqrt(s[[1]])*sqrt(s[[3]])*50*150),
                                               (2*cor.AWC2.3*sqrt(s[[2]])*sqrt(s[[3]])*100*150),
                                               (2*cor.AWC1.4*sqrt(s[[1]])*sqrt(s[[4]])*50*300),
                                               (2*cor.AWC2.4*sqrt(s[[2]])*sqrt(s[[4]])*100*300),
                                               (2*cor.AWC3.4*sqrt(s[[3]])*sqrt(s[[4]])*150*300),
                                               (2*cor.AWC1.5*sqrt(s[[1]])*sqrt(s[[5]])*50*400),
                                               (2*cor.AWC2.5*sqrt(s[[2]])*sqrt(s[[5]])*100*400),
                                               (2*cor.AWC3.5*sqrt(s[[3]])*sqrt(s[[5]])*150*400),
                                               (2*cor.AWC4.5*sqrt(s[[4]])*sqrt(s[[5]])*300*400),
                                               (2*cor.AWC1.6*sqrt(s[[1]])*sqrt(s[[6]])*50*(s[[7]]-1000)),
                                               (2*cor.AWC2.6*sqrt(s[[2]])*sqrt(s[[6]])*100*(s[[7]]-1000)),
                                               (2*cor.AWC3.6*sqrt(s[[3]])*sqrt(s[[6]])*150*(s[[7]]-1000)),
                                               (2*cor.AWC4.6*sqrt(s[[4]])*sqrt(s[[6]])*300*(s[[7]]-1000)),
                                               (2*cor.AWC5.6*sqrt(s[[5]])*sqrt(s[[6]])*400*(s[[7]]-1000)), na.rm=TRUE)
    
    varAWC[s[[7]] > 2000]  <- sum((2500*s[[1]]),(10000*s[[2]]),(22500*s[[3]]),(90000*s[[4]]),(160000*s[[5]]),(1000000*s[[6]]),
                                (2*cor.AWC1.2*sqrt(s[[1]])*sqrt(s[[2]])*50*100),
                                (2*cor.AWC1.3*sqrt(s[[1]])*sqrt(s[[3]])*50*150),
                                (2*cor.AWC2.3*sqrt(s[[2]])*sqrt(s[[3]])*100*150),
                                (2*cor.AWC1.4*sqrt(s[[1]])*sqrt(s[[4]])*50*300),
                                (2*cor.AWC2.4*sqrt(s[[2]])*sqrt(s[[4]])*100*300),
                                (2*cor.AWC3.4*sqrt(s[[3]])*sqrt(s[[4]])*150*300),
                                (2*cor.AWC1.5*sqrt(s[[1]])*sqrt(s[[5]])*50*400),
                                (2*cor.AWC2.5*sqrt(s[[2]])*sqrt(s[[5]])*100*400),
                                (2*cor.AWC3.5*sqrt(s[[3]])*sqrt(s[[5]])*150*400),
                                (2*cor.AWC4.5*sqrt(s[[4]])*sqrt(s[[5]])*300*400),
                                (2*cor.AWC1.6*sqrt(s[[1]])*sqrt(s[[6]])*50*1000),
                                (2*cor.AWC2.6*sqrt(s[[2]])*sqrt(s[[6]])*100*1000),
                                (2*cor.AWC3.6*sqrt(s[[3]])*sqrt(s[[6]])*150*1000),
                                (2*cor.AWC4.6*sqrt(s[[4]])*sqrt(s[[6]])*300*1000),
                                (2*cor.AWC5.6*sqrt(s[[5]])*sqrt(s[[6]])*400*1000), na.rm=TRUE)
    varAWC
}


### before applying it, I need to divide in tiles, because it is too demanding in memory
### Create the folder directory
Uncertainty_dir <-  paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth/tiles_var/")
dir.create(Uncertainty_dir)
setwd(Uncertainty_dir)

### Define the size of the blocks --- At each raster row do we start and finish each crop?
bs <- blockSize(awcRvar.s, minblocks=30)
variables <- names(awcRvar.s)

### I crop all raster files (across variables stacks) in a parallel process
### with a nested %:% and %dopar% from the foreach package (parallel nested for loops)

### identify how many cores I have in the computer
detectCores()
cl <- makeCluster(7)   ### Create cluster
registerDoParallel(cl)
getDoParWorkers()

predictor_tiles <- foreach(i=1:bs$n, .packages="raster") %:% foreach(variables=1:7, .packages="raster") %dopar% {
    ### What will be the name of the directory for this block?
    tile_Dir <- paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth/tiles_var/block",i)
    ### create the dir if it does not exist
    dir.create(tile_Dir)
    ### Set the working directory where I want to store the output
    setwd(tile_Dir)
    ### The result of one of the combinations is saved into the object "a" 
    a <- crop(awcRvar.s[[variables]], extent(awcRvar.s, bs$row[[i]], bs$row[[i]]+bs$nrows[[i]], 1, 11120),
              filename=paste0(names(awcRvar.s[[variables]]),"_tile_",i,".tif"),bylayer=TRUE, format="GTiff", overwrite=TRUE)
    gc()       ### clean garbage
    return(a)  ### when it is written this way, it is important to return the object we are interested in, to include it in the list returned by foreach
}

stopCluster(cl)

sum(is.na(getValues(predictor_tiles[[86]][[1]])))
sum(!is.na(getValues(predictor_tiles[[86]][[1]]))) ### that tile is already empty
predictor_tiles <- predictor_tiles[-c(86:93)]


### Life is easy... let's do in a for loop because I am tired

### fucntion to apply the tylor analysis
ff <- function(x) calc(x, AWCr_D_Var)
### I want to store the output in a list
AWCr_mm_0_200.VAR.tiles <- list()

tic <- Sys.time()
for(i in 1:length(predictor_tiles)){
    print(i)
    ### What will be the name of the directory for this block?
    tile_Dir <- paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth/tiles_var/block",i)
    ### Set the working directory where I want to store the output
    setwd(tile_Dir)
    
    ### create my stack now
    tile.stack <- stack(predictor_tiles[[i]])
    
    beginCluster(7)
    AWCr_mm_0_200.VAR.tile <- clusterR(tile.stack, ff, 
                                       export =list("AWCr_D_Var", "cor.AWC1.2","cor.AWC1.3","cor.AWC1.4" ,"cor.AWC1.5","cor.AWC1.6",
                                                    "cor.AWC2.3","cor.AWC2.4","cor.AWC2.5","cor.AWC2.6","cor.AWC3.4","cor.AWC3.5",
                                                    "cor.AWC3.6","cor.AWC4.5","cor.AWC4.6","cor.AWC5.6"),
                                       filename = paste0("AWCr_mm_0_200.VAR_tile_",i,".tif"), format = "GTiff",
                                       na.rm=T, inf.rm=T, progress = "text", overwrite = T)
    endCluster()
    
    AWCr_mm_0_200.VAR.tiles[[i]] <- AWCr_mm_0_200.VAR.tile
    
    rm(AWCr_mm_0_200.VAR.tile, i, tile_Dir,tile.stack)
    gc()
    
}

tac <- Sys.time()
tac-tic

save.image(paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth/11-AWC_TotalDepth.RData"))
### Let's do a mosaic with all of them!

## Assign function to mosaic
AWCr_mm_0_200.VAR.tiles$fun <- mean

## Create mosaic for whole France
setwd(paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth"))
AWCr_mm_0_200.VAR <- do.call(mosaic, AWCr_mm_0_200.VAR.tiles)
writeRaster(AWCr_mm_0_200.VAR, filename="AWCr_mm_0_200.VAR.tif", format="GTiff", overwrite=TRUE)
AWCr_mm_0_200.VAR <- raster("AWCr_mm_0_200.VAR.tif")
AWCr_mm_0_200 <- raster("AWCr_mm_0_200.tif")
### Calculate the SD
AWCr_mm_0_200.sd <- calc(AWCr_mm_0_200.VAR, sqrt,filename="AWCr_mm_0_200.sd.tif", format="GTiff", overwrite=TRUE)
AWCr_mm_0_200.sd <- raster("AWCr_mm_0_200.sd.tif")

france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(x =france,CRSobj = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"  )

AWCr_mm_0_200.sd <- extend(x =AWCr_mm_0_200.sd, y=AWCr_mm_0_200)
summary(AWCr_mm_0_200.sd/AWCr_mm_0_200)
Cv.s <- stack(AWCr_mm_0_200.sd, AWCr_mm_0_200)
plot(Cv.s)

### As it is, we could apply it like
relative_error <- function(s) { s[1]/s[2]}
ff <- function(x) calc(x, relative_error)
beginCluster(8)
rel_error_AWC <- clusterR(Cv.s, ff, export = list("relative_error"),
                          filename = "rel_error_AWC.tif", format = "GTiff",
                          na.rm=T, inf.rm=T, progress = "text", overwrite = T)
endCluster()

rel_error.v <- getValues(rel_error_AWC)
rel_error.na <- rel_error.v[!is.na(rel_error.v)]
median(rel_error.na, na.rm=TRUE)
min(rel_error.na, na.rm=TRUE)
summary(rel_error.na, na.rm=TRUE)
quantile(rel_error.na,probs=0.05)
quantile(rel_error.na,probs=0.99)
hist(rel_error.na)

par(mfrow=c(1,3), oma=c(1,3,1,4), mar=c(3,4,2,4))
jpeg("Figure6.b.jpeg", width=3740, height= 1700, units="px")
par(mar=c(5,5,5,5), oma=c(5,5,5,5), mfrow=c(1,3))
plot(AWCr_mm_0_200, breaks=c(0,50,90,130,200,400),
     col=c( "#B2182B", "#EF8A62" , "#D1E5F0", "#67A9CF", "#2166AC"), legend=FALSE)
legend("bottom", title= "AWC (mm)" , ncol=3,legend= c( "< 50", "50 - 90", "90 - 130", "130 - 200", "> 200"),
       fill =c( "#B2182B", "#EF8A62" , "#D1E5F0", "#67A9CF", "#2166AC"), border="white", bty="n", cex=6)
#plot(AWCr_mm_0_200, breaks=c(0,20,40,60,80,100,120,150,200,250,300,350,400), col=rev(viridis_pal(option="D")(12)))
plot(AWCr_mm_0_200.sd,breaks=c(0,20,40,60,80,103), col=rev(viridis_pal(option="A")(5)),legend=FALSE )
legend("bottom",title= "AWC SD (mm)", ncol=3,legend= c( "0 - 20", "20 - 40", "40 - 60", "60 - 80", "80 - 102"),
       fill =rev(viridis_pal(option="A")(5)), border="white", bty="n", cex=6)
plot(rel_error_AWC,breaks=c(0,0.1,0.2,0.3, 0.5,0.75,1), col=rev(viridis_pal(option="C")(6)), legend=FALSE)
legend("bottom",title= "CV (%)", ncol=3,legend= c( "0 - 10", "10 - 20", "20 - 30", "30 - 50", "50 - 75", "75 - 100"),
       fill =rev(viridis_pal(option="A")(6)), border="white", bty="n", cex=6)
dev.off()

##Clean the workspace
rm(cl, bs, depths, Cv.s, france, rance, rel_error.v, tac, tic, variables)
save.image(paste0(HomeDir,"Clean_Output/11-AWC_TotalDepth/11-AWC_TotalDepth.RData"))
### Here, for the article is enough. AWC and var(AWC), without considering the uncertainty for DEPTH

###########################################################################################################################################################