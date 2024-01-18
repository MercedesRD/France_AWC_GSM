
#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 0-5 cm
###
###  Check which factor levels of categorical variables present in the rasters are not found in the calibration dataset and set them to NA
###  Use clusterR to predict for the continuous metropolitan France surface

###  Author: Mercedes Roman Dobarco
###  Date: 6/12/2017


####### Load packages
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
library(Cubist)
library(foreach)
library(doParallel)
library(ithir)
library(raster)
library(snow)
library(doSNOW)


### set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

### Load data with covariates and cubist models
load(paste0(HomeDir,"Clean_Output/4.Cubist/4.1-Cubist_0_5.RData"))

########################################################################################################################################

### what do I need?
## covariates
## limon.Cub.0_5
## argile.Cub.0_5

### Copy two stacks, one for granulo, the other for coarse elements
covariates_EG <- covariates
covariates_granulo <- covariates


###### Prepare raster stack with covariates

# #### Prediction on metropolitan France
# names(covariates_granulo)
# 
# ### bdforet  
# unique(levels(granulo.data_0_5$bdforet))        # [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "30" "40" "50" All classes present
# ### and for coarse elements
# unique(levels(EG.data_0_5$bdforet)) 
# ### [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "30" "40" "50"
# 
# ### clc06
# setdiff(unique(covariates_granulo[["clc06"]]) , unlist(unique(levels(granulo.data_0_5$clc06))) )
# ### [1] 123 212 244 335 412 423 521 522 523                                              ### 9 classes are missing
# setdiff(unique(covariates_EG[["clc06"]]) , unlist(unique(levels(EG.data_0_5$clc06))) )
# ### [1] 111 123 212 213 244 335 412 421 423 521 522 523                                  ### 12 classes for Eg.granulo
# 
# ### For clc06 assign NA to those factor levels that are missing in the Cubist model
# ### it is quicker to run it in parallel, with the function clusterR
# ### mask with clusterR
# ### my function to use inside clusterR(calc...)
# dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates")
# dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_granulo")
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_granulo")
# 
# mask_clc <- function(y) {y[y == 123 | y == 212 | y == 244 | y == 335 | y == 412 | y == 423 | y == 521 | y == 522 | y == 523] <- NA ; return(y)} 
# 
# beginCluster(6)
# covariates_granulo[["clc06"]] <- clusterR(covariates_granulo[["clc06"]], calc, args = list(mask_clc),  filename="clc06.tif", type="GTiff", 
#                                        na.rm=T, inf.rm=T, progress = "text", overwrite = T)
# endCluster()
# plot(covariates_granulo[["clc06"]])
# rm(mask_clc)
# 
# summary(argile.Cub.0_5)
# summary(limon.Cub.0_5)
# 
# dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_EG")
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_EG")
# 
# mask_clc <- function(y) {y[y == 111 | y == 123 | y == 212 | y == 213 |y == 244 | y == 335 | y == 412 | y == 421 | y == 423 | y == 521 | y == 522 | y == 523] <- NA ; return(y)} 
# 
# beginCluster(6)
# covariates_EG[["clc06"]] <- clusterR(covariates_EG[["clc06"]], calc, args = list(mask_clc),  filename="clc06.tif", type="GTiff", 
#                                   na.rm=T, inf.rm=T, progress = "text", overwrite = T)
# endCluster()
# plot(covariates_EG[["clc06"]])
# rm(mask_clc)
# 
# summary(coarse.Cub.0_5) ### I can take out clc_06 form the model
# ### Refit the model in script 4.1Cubist_0_5.R
# 
# ### For ecoclim
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_granulo")
# setdiff(unique(covariates_granulo[["ecoclim"]]) , unlist(unique(levels(granulo.data_0_5$ecoclim))) )
#  ### Class 3 is missing ## and there is a weird value introduced!!! 1464736101
# ## my mask function for 1464736101
# mask_ecoclim <- function(y) { y[ is.nan(y) | y == 1464736101 | y == 3 ] <- NA ; return(y)} 
# 
# beginCluster(6)
# covariates_granulo[["ecoclim"]] <- clusterR(covariates_granulo[["ecoclim"]], calc, args = list(mask_ecoclim),
#                                          filename="ecoclim.tif", type="GTiff", na.rm=T, inf.rm=T, progress = "text", overwrite=T)
# endCluster()
# plot(covariates_granulo[["ecoclim"]])
# rm(mask_ecoclim)
# 
# ### EG
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_EG")
# setdiff(unique(covariates_EG[["ecoclim"]]) , unlist(unique(levels(EG.data_0_5$ecoclim))) )
# mask_ecoclim <- function(y) { y[ is.nan(y) | y == 1464736101 | y == 3 ] <- NA ; return(y)} 
# 
# beginCluster(6)
# covariates_EG[["ecoclim"]] <- clusterR(covariates_EG[["ecoclim"]], calc, args = list(mask_ecoclim),
#                                             filename="ecoclim.tif", type="GTiff", na.rm=T, inf.rm=T, progress = "text", overwrite=T)
# endCluster()
# plot(covariates_EG[["ecoclim"]])
# rm(mask_ecoclim)
# 
# ### For mat11 
# setdiff(unique(covariates_granulo[["mat11"]]) , unlist(unique(levels(granulo.data_0_5$mat11))) )### All classes present :)
# setdiff(unique(covariates_EG[["mat11"]]) , unlist(unique(levels(EG.data_0_5$mat11))) )### All classes present :)
# 
# ### For modelclc
# setdiff(unique(covariates_granulo[["modelclc"]]) , unlist(unique(levels(granulo.data_0_5$modelclc))) )### all present in the training dataset
# setdiff(unique(covariates_EG[["modelclc"]]) , unlist(unique(levels(EG.data_0_5$modelclc))) )### all present in the training dataset
# 
# ### For soil1
# setdiff(unique(covariates_granulo[["soil1"]]) , unlist(unique(levels(granulo.data_0_5$soil1))) )### class 25 (glacier) not present in training dataset
# setdiff(unique(covariates_EG[["soil1"]]) , unlist(unique(levels(EG.data_0_5$soil1))) )### class 25 (glacier) not present in training dataset
# 
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_granulo")
# 
# mask_soil1 <- function(y) { y[y == 25] <- NA ; return(y)} 
# beginCluster(6)
# covariates_granulo[["soil1"]] <- clusterR(covariates_granulo[["soil1"]], calc, args = list(mask_soil1), filename="soil1.tif", type="GTiff",
#                                        na.rm=T, inf.rm=T, progress = "text", overwrite=T)
# endCluster()
# plot(covariates_granulo[["soil1"]])
# rm(mask_soil1)
# 
# setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Covariates/Prediction_tiff_EG")
# mask_soil1 <- function(y) { y[y == 25] <- NA ; return(y)} 
# beginCluster(6)
# covariates_EG[["soil1"]] <- clusterR(covariates_EG[["soil1"]], calc, args = list(mask_soil1), filename="soil1.tif", type="GTiff",
#                                           na.rm=T, inf.rm=T, progress = "text", overwrite=T)
# endCluster()
# plot(covariates_EG[["soil1"]])
# rm(mask_soil1)
# 

### Or directly loading them form the folder
setwd(paste0(HomeDir,"Covariates/Prediction_tiff_granulo"))
covariates_granulo[["clc06"]] <- raster("clc06.tif")
covariates_granulo[["ecoclim"]] <- raster("ecoclim.tif")
covariates_granulo[["soil1"]] <- raster("soil1.tif")

# setwd(paste0(HomeDir,"Covariates/Prediction_tiff_EG"))
# covariates_EG[["clc06"]] <- raster("clc06.tif")
# covariates_EG[["ecoclim"]] <- raster("ecoclim.tif")
# covariates_EG[["soil1"]] <- raster("soil1.tif")

# ### eliminate "clc_06" from covariates_EG
# covariates_EG <- dropLayer(covariates_EG, 2)
# names(covariates_EG)
# names(covariates_granulo)


#### Convert some raster layers into factors
for (i in c("bdforet","clc06","ecoclim","mat11","modelclc","soil1")){
    print(i)                    ### What layer are we working in?
    print(is.factor(covariates_granulo[[i]]))    ### Is it a factor?
    covariates_granulo[[i]] <- as.factor(covariates_granulo[[i]]) ### Convert into factor, and keep it in the same stack
    print(is.factor(covariates_granulo[[i]]))    ### Did it change into factor?
    gc()                        ### Garbage clean!
}

for (i in c("bdforet","ecoclim","mat11","modelclc","soil1")){
    print(i)                    ### What layer are we working in?
    print(is.factor(covariates_EG[[i]]))    ### Is it a factor?
    covariates_EG[[i]] <- as.factor(covariates_EG[[i]]) ### Convert into factor, and keep it in the same stack
    print(is.factor(covariates_EG[[i]]))    ### Did it change into factor?
    gc()                        ### Garbage clean!
}



factor_list.granulo <- list()
elements <- c("bdforet","clc06","ecoclim","mat11","modelclc","soil1")
for (kk in 1:length(elements)){
    factor_list.granulo[[kk]] <- levels(covariates_granulo[[elements[kk]]])
}
names(factor_list.granulo) <- c("bdforet","clc06","ecoclim","mat11","modelclc","soil1")

is.factor(covariates_granulo)
unique(levels(covariates_granulo))
names(covariates_granulo)


factor_list.EG <- list()
elements <- c("bdforet","ecoclim","mat11","modelclc","soil1")
for (kk in 1:length(elements)){
    factor_list.EG[[kk]] <- levels(covariates_EG[[elements[kk]]])
}
names(factor_list.EG) <- c("bdforet","ecoclim","mat11","modelclc","soil1")

is.factor(covariates_EG)
unique(levels(covariates_EG))
names(covariates_EG)

######################################################################################################################################################
###  Save this raster stack, to use with the other layers
setwd(paste0(HomeDir,"Covariates/RData"))

save(covariates_granulo, file="covariates_granulo_0_5.RData")
save(factor_list.granulo, file="factor_list.granulo_0_5.RData")
save(covariates_EG, file="covariates_EG_0_5.RData")
save(factor_list.EG, file="factor_list.EG_0_5.RData")

rm(kk, elements)

######################################################################################################################################################

### Predict cubist models
dir.create(paste0(HomeDir,"Clean_Output/5-Cubist_preds"))
set.wd(paste0(HomeDir,"Clean_Output/5-Cubist_preds"))

beginCluster(7)
system.time(silt.alr.cub.0_5 <- clusterR(covariates_granulo, predict, args = list(limon.Cub.0_5, committees = 20),
                           filename = "silt.alr.cub.0_5.tif", format = "GTiff",
                           factors=factor_list.granulo, na.rm=T, inf.rm=T, progress = "text", overwrite = T))
endCluster()


beginCluster(7)
system.time(argile.alr.cub.0_5 <- clusterR(covariates_granulo, predict, args = list(argile.Cub.0_5, committees = 20),
                            filename = "argile.alr.cub.0_5.tif", format = "GTiff",
                            factors=factor_list.granulo, na.rm=T,inf.rm=T, progress = "text", overwrite = T))
endCluster()


### Later, we predicted coarse elements with quantile regression forests instead
# beginCluster(7)
# system.time(coarse.cub.0_5 <- clusterR(covariates_EG, predict, args = list(coarse.Cub.0_5, committees = 20),
#                                            filename = "coarse.cub.0_5.tif", format = "GTiff",
#                                            factors=factor_list.EG, na.rm=T,inf.rm=T, progress = "text", overwrite = T))
# endCluster()

par(mfrow=c(1,3))
plot(silt.alr.cub.0_5)
plot(argile.alr.cub.0_5)
#plot(coarse.cub.0_5)

#########################################################################################################################################################
### tranform from raster to grid (ASCII)
silt.alr.cub.0_5 <- raster("silt.alr.cub.0_5.tif" )
argile.alr.cub.0_5 <- raster( "argile.alr.cub.0_5.tif")
#coarse.cub.0_5 <- raster("coarse.cub.0_5.tif" )

extent(silt.alr.cub.0_5)
res(silt.alr.cub.0_5)

########################################################################################################################################################

### End of the script