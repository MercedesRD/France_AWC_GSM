#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  IGCS data on particle size distribution and coarse elements
###          
###  Depth 100-200 cm
###
###  Silt-alr and Sand-alr as response variables
###  Regression model: Cubist because it is faster than gbm and randomForest
###  Spatial distribution of silt-alr and sand-alr
###  Variogram and cross-variograms
###  LMCR of Silt-alr and Sand-alr
###  Universal co-kriging using as trend the prediction of the regression model
###  Validation statistics with RMQs composite data
###  The script must check if all factor levels are present in the training dataset

###  Author: Mercedes Roman Dobarco
###  Date: 09/05/2017

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
library(compositions)
library(rgr)
library(randomForest)
library(Cubist)
library(gbm)
library(foreach)
library(doParallel)
library(ithir)
library(raster)
library(soiltexture)

### set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
dir.create(paste0(HomeDir,"Clean_Output/4.Cubist"))
setwd(paste0(HomeDir,"Clean_Output/4.Cubist"))

### Load the data
load(paste0(HomeDir,"Clean_Output/3.2-Covariates/igcs_granulo_gsm.RData"))
rm
rm(eg_cor2, igcs_granulo_gsm_bck, texture_gsm.m)


###############################################################################################################################################
#################                         Correct texture data
###############################################################################################################################################

### In 6699 cases, the sum of the three fractions is greater than 1000
dim(igcs_granulo_gsm[(igcs_granulo_gsm$argile_100_200 + igcs_granulo_gsm$sable_100_200 + igcs_granulo_gsm$limon_100_200  )>1000 & !is.na(igcs_granulo_gsm$argile_100_200),])[1]
dim(igcs_granulo_gsm[(igcs_granulo_gsm$argile_100_200 + igcs_granulo_gsm$sable_100_200 + igcs_granulo_gsm$limon_100_200  )<0 & !is.na(igcs_granulo_gsm$argile_100_200),])[1]

### for avoiding -Inf and Inf values for the alr transform, I give infinitesimal values to sand , silt, and clay when they are 0
igcs_granulo_gsm[igcs_granulo_gsm$sable_100_200<=0 & !is.na(igcs_granulo_gsm$argile_100_200),]$sable_100_200 <- 0.001
igcs_granulo_gsm[igcs_granulo_gsm$argile_100_200<=0 & !is.na(igcs_granulo_gsm$argile_100_200),]$argile_100_200 <- 0.001
igcs_granulo_gsm[igcs_granulo_gsm$limon_100_200<=0 & !is.na(igcs_granulo_gsm$argile_100_200),]$limon_100_200 <- 0.001

### copy backup
igcs_granulo_gsm.bck <- igcs_granulo_gsm

### I normalize with the package soiltexture
tri.igcs_granulo_100_200 <- igcs_granulo_gsm[,c("argile_100_200", "limon_100_200", "sable_100_200")] ### subset the columns of interest
tri.igcs_granulo_100_200 <- tri.igcs_granulo_100_200/10 ### Transform to %
colnames(tri.igcs_granulo_100_200) <- c("CLAY", "SILT","SAND") ### Change names as soiltexture likes
norm.tri <- TT.normalise.sum(tri.igcs_granulo_100_200[!is.na(tri.igcs_granulo_100_200$CLAY), ], text.tol=0.01 ) ### Normalize to sum 100 %
norm.tri <- norm.tri*10 ## Transform to g/kg
head(norm.tri)
head(igcs_granulo_gsm[!is.na(igcs_granulo_gsm$argile_100_200),c("argile_100_200", "limon_100_200", "sable_100_200")])

### Replace new values in old dataframe
igcs_granulo_gsm[!is.na(igcs_granulo_gsm$sable_100_200),]$sable_100_200 <- norm.tri$SAND
igcs_granulo_gsm[!is.na(igcs_granulo_gsm$sable_100_200),]$argile_100_200 <- norm.tri$CLAY
igcs_granulo_gsm[!is.na(igcs_granulo_gsm$sable_100_200),]$limon_100_200 <- norm.tri$SILT

plot(igcs_granulo_gsm$argile_100_200,igcs_granulo_gsm.bck$argile_100_200 )
plot(igcs_granulo_gsm$sable_100_200,igcs_granulo_gsm.bck$sable_100_200 )
plot(igcs_granulo_gsm$limon_100_200,igcs_granulo_gsm.bck$limon_100_200 )


### chekc this did not mess with the other data
sum(!is.na(igcs_granulo_gsm$coarse_100_200)) ## 44802
sum(!is.na(igcs_granulo_gsm$argile_100_200)) ## 13086
sum(!is.na(igcs_granulo_gsm$sable_100_200)) ## 13086
sum(!is.na(igcs_granulo_gsm$limon_100_200)) ## 13086

rm(tri.igcs_granulo_100_200, norm.tri)

###########################################################################################################################################
#######################         Perform the ALR transformation
###########################################################################################################################################

# alrTransform <- function(table, clay, silt, sand, denominator) {
#     ### Transform into a matrix the subset of texture variables
#     texture.mat <- as.matrix(table[ ,c(clay,silt,sand)])
#     ### Get the dimension
#     nrows <- dim(texture.mat)[1]
#     ### Dividing by the denominator
#     texture.alr <- alr(texture.mat,denominator)
#     ### Create a data.frame with alr-transformed data
#     texture.alr.2 <- data.frame(matrix(texture.alr, nrow=nrows, byrow=F))
#     colnames(texture.alr.2) <- c(paste0("alr.",colnames(texture.alr)[[1]]),paste0("alr.",colnames(texture.alr)[[2]]))
#     ### attach to the original data.frame
#     table <- cbind(table, texture.alr.2)
#     return (table)
# }

texture.mat <- as.matrix(igcs_granulo_gsm[ ,c("argile_100_200", "limon_100_200", "sable_100_200")])
texture.alr <- alr(texture.mat,3)
nrows <- dim(texture.alr)[1]
texture.alr.2 <- data.frame(matrix(texture.alr, nrow=nrows, byrow=F))
colnames(texture.alr.2) <- c(paste0("alr.",colnames(texture.alr)[[1]]),paste0("alr.",colnames(texture.alr)[[2]]))

igcs_granulo_gsm.alr <- igcs_granulo_gsm
igcs_granulo_gsm.alr$alr.argile_100_200 <- NA
igcs_granulo_gsm.alr$alr.limon_100_200 <- NA

length(igcs_granulo_gsm.alr[!is.na(igcs_granulo_gsm$argile_100_200),]$alr.argile_100_200)
length(texture.alr.2$alr.argile_100_200)

igcs_granulo_gsm.alr[!is.na(igcs_granulo_gsm$argile_100_200),]$alr.argile_100_200 <- texture.alr.2$alr.argile_100_200
igcs_granulo_gsm.alr[!is.na(igcs_granulo_gsm$argile_100_200),]$alr.limon_100_200 <- texture.alr.2$alr.limon_100_200
rm(texture.mat,texture.alr, nrows, texture.alr.2)

# igcs_granulo_gsm.alr <- alrTransform(igcs_granulo_gsm, "argile_100_200", "limon_100_200", "sable_100_200", 3)

### check for outliers
plot(igcs_granulo_gsm.alr$alr.argile_100_200, igcs_granulo_gsm.alr$alr.limon_100_200)
#### check correlation in the feature space
cor.test(igcs_granulo_gsm.alr$alr.argile_100_200, igcs_granulo_gsm.alr$alr.limon_100_200) ### Strong correlation
#### distribution?
par(mfrow=c(1,2))
hist(igcs_granulo_gsm.alr$alr.argile_100_200, col="cyan")
hist(igcs_granulo_gsm.alr$alr.limon_100_200, col="deeppink1") ### pretty normal

### check correlation with coarse elements:
cor.test(igcs_granulo_gsm.alr$coarse_100_200, igcs_granulo_gsm.alr$alr.limon_100_200) 
cor.test(igcs_granulo_gsm.alr$coarse_100_200, igcs_granulo_gsm.alr$alr.argile_100_200) 

par(mfrow=c(1,3))
plot(igcs_granulo_gsm.alr$alr.argile_100_200, igcs_granulo_gsm.alr$alr.limon_100_200, 
     xlab="Clay-alr (100-200 cm)", ylab="Silt-alr (100-200 cm)")
plot(igcs_granulo_gsm.alr$coarse_100_200,igcs_granulo_gsm.alr$alr.argile_100_200,
     ylab="Clay-alr (100-200 cm)", xlab="Log coarse fragments (100-200 cm)")
plot(igcs_granulo_gsm.alr$coarse_100_200, igcs_granulo_gsm.alr$alr.limon_100_200,
     ylab="Silt-alr (100-200 cm)", xlab="Log coarse fragments (100-200 cm)")




##############################################################################################################################################
##############                                cubist models for silt-alr and clay-alr                                           ##############
##############                                   Regression + Universal co-kriging                                              ##############
##############################################################################################################################################

#### Select only columns of interest
#### Data for the 100-200 cm depth

## Are all levels from clc06 in the calibration dataset?
unique(getValues(covariates[["clc06"]]))
# [1]  NA 523 423 331 322 311 231 211 112 121 142 512 123 111 242 122 511 133 141 321 324 411 421 312 131 124 243 222 313 132 412 221 522 332 333 521
# [37] 334 422 335 212 323 223 244 213

### I'll have to mask some values for the predict part
#igcs_granulo_gsm$clc06 <- droplevels(igcs_granulo_gsm$clc06)

# 211 311 221 242 231 313 243 333 324 312 112 121 142 122 512 213 411 222 321 322 131 422 323 331 511 421 124 133 132 223 111 141 332 334
### I eliminate clc form the calibration dataset

granulo.data_100_200 <- igcs_granulo_gsm.alr[,c("id_profil","x","y","alr.argile_100_200","alr.limon_100_200",
                                     "bdforet","clc06","cti","curv_long","curv_trans","curvature","ecoclim",
                                     "eros","etpMax","etpMean","etpMedian","etpMin","EVI_median_jan","EVI_median_june",
                                     "exposition","graviAL2000","graviBg2000","graviG2000","hli","idpr","linear_aspect",
                                     "mat11","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                     "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope",
                                     "slopeascos","slopeassin","slopeastrasp","soil1","srr","srtm","tMeanMax",
                                     "tMeanMean","tMeanMedian","tMeanMin")]
## Only want complete.cases
granulo.data_100_200 <- granulo.data_100_200[complete.cases(granulo.data_100_200),]
### drop levels of factors
granulo.data_100_200[] <- lapply(granulo.data_100_200, function(x) if (is.factor(x)) x[, drop=TRUE] else x)
setdiff(unique(covariates[["clc06"]]) , unlist(unique(levels(granulo.data_100_200$clc06))) )
# [1] 111 123 141 212 244 335 412 423 521 522 523

### 11 classes are not in the calibration data. I eliminate clc06 from the calibration model

### Dataframe with predictors
predictors <- granulo.data_100_200[,c("bdforet","cti","curv_long","curv_trans","curvature","ecoclim",
                           "eros","etpMax","etpMean","etpMedian","etpMin","EVI_median_jan","EVI_median_june",
                           "exposition","graviAL2000","graviBg2000","graviG2000","hli","idpr","linear_aspect",
                           "mat11","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                           "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope",
                           "slopeascos","slopeassin","slopeastrasp","soil1","srr","srtm","tMeanMax",
                           "tMeanMean","tMeanMedian","tMeanMin")]

### Fit a cubist model, for clay-alr
set.seed(1948)
argile.Cub.100_200 <- cubist(x=predictors, y=granulo.data_100_200$alr.argile_100_200, cubistControl( extrapolation = 5, unbiased=TRUE), committees = 20)

set.seed(1948)
limon.Cub.100_200 <- cubist(x=predictors, y=granulo.data_100_200$alr.limon_100_200, cubistControl( extrapolation = 5, unbiased=TRUE), committees = 20)

# summary(argile.Cub.100_200)
# summary(limon.Cub.100_200)                                  

### Predict at training sites, used for universal co-kriging
granulo.data_100_200$pred.cub.alr.argile_100_200 <- predict(argile.Cub.100_200, newdata=predictors)
granulo.data_100_200$pred.cub.alr.limon_100_200  <- predict(limon.Cub.100_200, newdata=predictors)

### Calculate residuals
granulo.data_100_200$resid.cub.alr.argile_100_200 <- granulo.data_100_200$alr.argile_100_200 - granulo.data_100_200$pred.cub.alr.argile_100_200 
granulo.data_100_200$resid.cub.alr.limon_100_200  <- granulo.data_100_200$alr.limon_100_200  - granulo.data_100_200$pred.cub.alr.limon_100_200

argile.alr.cub.cal_100_200 <- goof(observed = granulo.data_100_200$alr.argile_100_200, predicted = granulo.data_100_200$pred.cub.alr.argile_100_200, plot.it = TRUE, type="DSM")
#     R2 concordance   MSE  RMSE   bias
#  0.502       0.578 1.058 1.029 -0.001
round(argile.alr.cub.cal_100_200,digits=3)
limon.alr.cub.cal_100_200 <- goof(observed = granulo.data_100_200$alr.limon_100_200, predicted = granulo.data_100_200$pred.cub.alr.limon_100_200, plot.it = TRUE, type="DSM")    
#     R2 concordance   MSE  RMSE   bias
#  0.516       0.618 0.831 0.912    0
round(limon.alr.cub.cal_100_200,digits=3)

##########################################################################################################################

rm(predictors)

save(granulo.data_100_200,EG.data_100_200, coarse.cub.cal_100_200, argile.alr.cub.cal_100_200,limon.alr.cub.cal_100_200,igcs_granulo_gsm.alr, file="cubist.100_200.RData" )
save.image("4.6-Cubist_100_200.RData")

### Export the residuals in a table CSV
names(granulo.data_100_200)
table.granulo.data_100_200 <- granulo.data_100_200[,c( "id_profil", "x" ,"y","alr.argile_100_200", "alr.limon_100_200",
                                                     "pred.cub.alr.argile_100_200" , "pred.cub.alr.limon_100_200" ,
                                                     "resid.cub.alr.argile_100_200" ,"resid.cub.alr.limon_100_200")]
write.csv2(table.granulo.data_100_200, file = "table.granulo.data_100_200.csv")

######################################################################################################################
### end of the script