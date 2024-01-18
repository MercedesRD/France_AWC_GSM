#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 30-60 cm
###
###  Check which factor levels of categorical variables present in the rasters are not found in the calibration dataset and set them to NA
###  Use clusterR to predict for the continuous metropolitan France surface

###  Author: Mercedes Roman Dobarco
###  Date: 7/12/2017


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
load(paste0(HomeDir,"Clean_Output/4.Cubist/4.4-Cubist_30_60.RData"))

### Copy two stacks, one for granulo, the other for coarse elements
#covariates_EG <- covariates
covariates_granulo <- covariates

# ### check if I need to eliminate or not any variable (clc) from the coarse.Cub.30_60 model
# 
# ### clc06
# setdiff(unique(covariates_granulo[["clc06"]]) , unlist(unique(levels(granulo.data_30_60$clc06))) )
# ### [1] 123 212 244 335 412 423 521 522 523                                              ### 9 classes are missing
# setdiff(unique(covariates_EG[["clc06"]]) , unlist(unique(levels(EG.data_30_60$clc06))) )
# ### [1] 111 123 212 213 244 335 421 423 521 522 523                                  ### 11 classes for Eg.granulo
# ### [1] 111 123 212 213 244 335 421 423 521 522 523

# ### bdforet  
# unique(levels(granulo.data_30_60$bdforet))        # [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "30" "40" "50" All classes present
# ### and for coarse elements
# unique(levels(EG.data_30_60$bdforet)) 
# ### [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "30" "40" "50"
#
# setdiff(unique(covariates_granulo[["ecoclim"]]) , unlist(unique(levels(granulo.data_30_60$ecoclim))) )
# setdiff(unique(covariates_EG[["ecoclim"]]) , unlist(unique(levels(EG.data_30_60$ecoclim))) )
### Level 3, as for 0-5 cm

# ### For mat11 
# setdiff(unique(covariates_granulo[["mat11"]]) , unlist(unique(levels(granulo.data_30_60$mat11))) )### All classes present :)
# setdiff(unique(covariates_EG[["mat11"]]) , unlist(unique(levels(EG.data_30_60$mat11))) )### All classes present :)
# 
# ### For modelclc
# setdiff(unique(covariates_granulo[["modelclc"]]) , unlist(unique(levels(granulo.data_30_60$modelclc))) )### all present in the training dataset
# setdiff(unique(covariates_EG[["modelclc"]]) , unlist(unique(levels(EG.data_30_60$modelclc))) )### all present in the training dataset
# 
# ### For soil1
# setdiff(unique(covariates_granulo[["soil1"]]) , unlist(unique(levels(granulo.data_30_60$soil1))) )### class 25 (glacier) not present in training dataset
# setdiff(unique(covariates_EG[["soil1"]]) , unlist(unique(levels(EG.data_30_60$soil1))) )### class 25 (glacier) not present in training dataset


### I went back to the script 4.4-Cubist_30_60 and eliminated clc06 from the predictor variables, to ncrease the predicted surface

########################################################################################################################################

### what do I need?
## limon.Cub.30_60
## argile.Cub.30_60
## coarse.Cub.30_60


### Or directly loading them form the folder
setwd(paste0(HomeDir,"Covariates/Prediction_tiff_granulo"))
covariates_granulo[["clc06"]] <- raster("clc06.tif")
covariates_granulo[["ecoclim"]] <- raster("ecoclim.tif")
covariates_granulo[["soil1"]] <- raster("soil1.tif")


#### Convert some raster layers into factors
for (i in c("bdforet","clc06","ecoclim","mat11","modelclc","soil1")){
    print(i)                    ### What layer are we working in?
    print(is.factor(covariates_granulo[[i]]))    ### Is it a factor?
    covariates_granulo[[i]] <- as.factor(covariates_granulo[[i]]) ### Convert into factor, and keep it in the same stack
    print(is.factor(covariates_granulo[[i]]))    ### Did it change into factor?
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


######################################################################################################################################################

### Predict cubist models
set.wd(paste0(HomeDir,"Clean_Output/5-Cubist_preds"))


beginCluster(7)
system.time(silt.alr.cub.30_60 <- clusterR(covariates_granulo, predict, args = list(limon.Cub.30_60, committees = 20),
                                         filename = "silt.alr.cub.30_60.tif", format = "GTiff",
                                         factors=factor_list.granulo, na.rm=T, inf.rm=T, progress = "text", overwrite = T))
endCluster()

beginCluster(7)
system.time(argile.alr.cub.30_60 <- clusterR(covariates_granulo, predict, args = list(argile.Cub.30_60, committees = 20),
                                           filename = "argile.alr.cub.30_60.tif", format = "GTiff",
                                           factors=factor_list.granulo, na.rm=T,inf.rm=T, progress = "text", overwrite = T))
endCluster()

par(oma=c(2,2,1,1), mar=c(3,3,1,1))
par(mfrow=c(1,3))
plot(silt.alr.cub.30_60)
plot(argile.alr.cub.30_60)

########################################################################################################################################################

### End of the script