###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date 06/03/2018

### Objective: Mapping coarse elements

####### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(raster)
library(lattice)
library(ggplot2)
library(plyr)
library(ranger)
library(quantregForest)
library(doParallel)
library(foreach)


### Load the calibration datasets
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
load(paste0(HomeDir,"Input/CalibrationDataCoarse.RData"))
rm("df.input","EG.data_0_5.r","EG.data_100_200.r","EG.data_15_30.r","EG.data_30_60.r","EG.data_5_15.r","EG.data_60_100.r")
dir.create(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))


### Conditions
### Response variable in the dataframe
responses <- c("coarse_0_5", "coarse_5_15", "coarse_15_30", "coarse_30_60", "coarse_60_100", "coarse_100_200")
paste0(responses,".qrF")

# 
# ### calibrate qrF models
# system.time(
# ### Now do my nested loop
# cl <- makeCluster(7)
# registerDoParallel(cl)

### Lets seeeee
Rocks.qrF <- foreach(d = 1:6, .packages = "quantregForest", .export = c("df.output", "responses")) %do% {
    ### Select the dataframe  for this depth
    EG.data <- df.output[[d]]
    ### The response variable name
    response <- responses[[d]]
    ### Dataframe with training data
    ### but this time, let's try just with the non-transformed data
    EG.data[,c(paste0("Unt.",response))] <- exp(EG.data[,response])
    EG.data.train <- EG.data[ ,c(paste0("Unt.",response),"bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                                 "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                                 "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                 "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                                 "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")]
    ### select only complete cases and drop empty factor levels
    EG.data.train <- EG.data.train[complete.cases(EG.data.train),]
    EG.data.train[] <- lapply(EG.data.train, function(x) if (is.factor(x)) x[, drop=TRUE] else x)
    
    ### Calibrate models
    set.seed(1984)
    qrf.train <- quantregForest(EG.data.train[,2:ncol(EG.data.train)], EG.data.train[ ,1] ,
                                ntree=1000, nodesize=20, keep.inbag = T,importance=T)
    return(qrf.train)
    }
# stopCluster(cl)
# )

names(Rocks.qrF) <- paste0(responses,".qrF")

round(importance(Rocks.qrF$coarse_0_5.qrF), digits=2)   
varImpPlot(Rocks.qrF$coarse_0_5.qrF,n.var =44,type=1)
varImpPlot(Rocks.qrF$coarse_5_15.qrF,n.var =44,type=1)
varImpPlot(Rocks.qrF$coarse_15_30.qrF,n.var =44,type=1)
varImpPlot(Rocks.qrF$coarse_30_60.qrF,n.var =44,type=1)
varImpPlot(Rocks.qrF$coarse_60_100.qrF,n.var =44,type=1)
varImpPlot(Rocks.qrF$coarse_100_200.qrF,n.var =44,type=1)

### Clean the image
rm(EG.data, EG.data.train, cl, d, response, qrf.train)

### check raster stack with predictors
load(paste0(HomeDir,"Clean_Output/CoarseElements/2.PredictorStack.RData"))
names(covariates4maps)

### Drop # 12, EVI_median_jan
covariates.qrF.PM <- dropLayer(covariates4maps, 12)
### Drop 20, mat11
covariates.qrF.PM <- dropLayer(covariates.qrF.PM, 20)
### Drop CHELSA
covariates.qrF.PM <- dropLayer(covariates.qrF.PM, 43:56)
covariates <- c("bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")
names(covariates.qrF.PM)

#####################################################################################################################################################

### If the previous lines do not work 
### because the data source no longer have the right path after changing the raster files with the covariates to your computer,
### reload the raster layers like this:

### since I changed the location of the covariates files in my computer i need to load them again
setwd(paste0(HomeDir,"Covariates/Tiff/"))
### Load all files
raster.files <- list.files(pattern=".tif$")
### Eliminate some of them
raster.files <- raster.files[-2]
raster.files <- raster.files[-12]
raster.files <- raster.files[-20]
raster.list <- lapply(raster.files, raster)
raster.names <- raster.files
raster.names <- gsub(pattern = ".tif", replacement = "", x = raster.names)
names(raster.list) <- raster.names
covariates.qrF.PM <- stack(raster.list)

### Two of them are different
setwd(paste0(HomeDir,"Covariates/Prediction_tiff_EG"))
covariates.qrF.PM[["ecoclim"]] <- raster("ecoclim.tif")
covariates.qrF.PM[["soil1"]] <- raster("soil1.tif")

### Add two more layers
setwd(paste0(HomeDir,"Input/Coarse_elements/BDGSF/GeoTiff"))
MAT_RF <- raster("MAT_RF.tif") 
AGL_RF <- raster("AGL_RF.tif")
covariates.qrF.PM <- stack(covariates,MAT_RF, AGL_RF)  

#####################################################################################################################################################

#### Convert some raster layers into factors
for (i in c("bdforet","ecoclim","modelclc","soil1", "MAT_RF", "AGL_RF")){
    print(i)                    ### What layer are we working in?
    print(is.factor(covariates.qrF.PM[[i]]))    ### Is it a factor?
    covariates.qrF.PM[[i]] <- as.factor(covariates.qrF.PM[[i]]) ### Convert into factor, and keep it in the same stack
    print(is.factor(covariates.qrF.PM[[i]]))    ### Did it change into factor?
    gc()                        ### Garbage clean!
}

factor_list.PM.EG <- list()
elements <- c("bdforet","ecoclim","modelclc","soil1", "MAT_RF", "AGL_RF")
for (kk in 1:length(elements)){
    factor_list.PM.EG[[kk]] <- levels(covariates.qrF.PM[[elements[kk]]])
}
names(factor_list.PM.EG) <- c("bdforet","ecoclim","modelclc","soil1", "MAT_RF", "AGL_RF")


### save image just in case
setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
save.image("qrF.models.CoarseElements.RData")

### Map them
rm(d, i, kk, elements, predictors, AGL_RF, chelsa.s, covariates.chelsa, MAT_RF)

for(d in 1: length(responses)){
    
    ### Which depth?
    print(d)
    
    ### Mean
    beginCluster(7)
    clusterR(covariates.qrF.PM, predict, args = list(Rocks.qrF[[d]], what=mean),
             filename = paste0(responses[[d]],".preds.mean.tif"), format = "GTiff",
                               factors=factor_list.PM.EG,  na.rm=T, inf.rm=T,
                               progress = "text", overwrite = T)
    endCluster()
    
    ### SD
    beginCluster(7)
    clusterR(covariates.qrF.PM, predict, args = list(Rocks.qrF[[d]], what=sd),
             filename = paste0(responses[[d]],".preds.sd.tif"), format = "GTiff",
                               factors=factor_list.PM.EG,  na.rm=T, inf.rm=T,
                               progress = "text", overwrite = T)
    endCluster()
    
    
    ### 5 percentile
    beginCluster(7)
    clusterR(covariates.qrF.PM, predict, args = list(Rocks.qrF[[d]], what=0.05),
             filename = paste0(responses[[d]],".preds.05p.tif"), format = "GTiff",
                               factors=factor_list.PM.EG,  na.rm=T, inf.rm=T,
                               progress = "text", overwrite = T)
    endCluster()
    
    ### 95 percentile
    beginCluster(7)
    clusterR(covariates.qrF.PM, predict, args = list(Rocks.qrF[[d]], what=0.95),
             filename = paste0(responses[[d]],".preds.95p.tif"), format = "GTiff",
                               factors=factor_list.PM.EG,  na.rm=T, inf.rm=T,
                               progress = "text", overwrite = T)
    endCluster()
    
    
}

### Load all files
raster.files <- list.files(pattern=".tif$")
raster.list <- lapply(raster.files, raster)
raster.names <- raster.files
raster.names <- gsub(pattern = ".tif", replacement = "", x = raster.names)
names(raster.list) <- raster.names

### Plot all mean and SD maps
library(scales)
library(viridis)

### Make plots
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(x =france,CRSobj = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"  )
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,3,4,2))
plot(raster.list$coarse_0_5.preds.mean, breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),
     col=rev(viridis_pal(option="D")(20)),
     main="Coarse elements 0-5 cm")
lines(france)


par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,3,4,2))
plot(raster.list$coarse_0_5.preds.05p, breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),
     col=rev(viridis_pal(option="D")(20)),
     main="Coarse elements LPL 0-5 cm")
lines(france)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,3,4,2))
plot(raster.list$coarse_0_5.preds.95p, breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),
     col=rev(viridis_pal(option="D")(20)),
     main="Coarse elements UPL 0-5 cm")
lines(france)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,3,4,2))
plot(raster.list$coarse_0_5.preds.sd, breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),
     col=rev(viridis_pal(option="D")(20)),
     main="Coarse elements SD 0-5 cm")
lines(france)

### end script