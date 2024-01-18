##########################################################################################################################################

############# validation of particle size and coarse elements predictions with RMQS measurements

### Date 12/06/2018
### Author: Mercedes Roman Dobarco

####### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(geoR)
library(lattice)
library(ggplot2)
library(Hmisc)
library(plyr)
library(soiltexture)
library(compositions)
library(rgr)
library(Cubist)
library(foreach)
library(doParallel)
library(ithir)
library(raster)
library(doBy)
library(debug)
library(testthat)

# ### Prepara dataframe observed-predicted values at RMQS sites -----------

### Load RMQS data
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"Clean_Output/3.5-RMQS_merge"))
load("RMQS.validation.RData")

### Load the predictions
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
clay.0_5 <- raster("argile.preds.RK.0_5.tif")
silt.0_5 <- raster("limon.preds.RK.0_5.tif")
sand.0_5 <- raster("sable.preds.RK.0_5.tif")

clay.5_15 <- raster("argile.preds.RK.5_15.tif")
silt.5_15 <- raster("limon.preds.RK.5_15.tif")
sand.5_15 <- raster("sable.preds.RK.5_15.tif")

clay.15_30 <- raster("argile.preds.RK.15_30.tif")
silt.15_30 <- raster("limon.preds.RK.15_30.tif")
sand.15_30 <- raster("sable.preds.RK.15_30.tif")

clay.30_60 <- raster("argile.preds.RK.30_60.tif")
silt.30_60 <- raster("limon.preds.RK.30_60.tif")
sand.30_60 <- raster("sable.preds.RK.30_60.tif")

clay.60_100 <- raster("argile.preds.RK.60_100.tif")
silt.60_100 <- raster("limon.preds.RK.60_100.tif")
sand.60_100 <- raster("sable.preds.RK.60_100.tif")

clay.100_200 <- raster("argile.preds.RK.100_200.tif")
silt.100_200 <- raster("limon.preds.RK.100_200.tif")
sand.100_200 <- raster("sable.preds.RK.100_200.tif")

### Load quantile random forest predictions for coarse elements
setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
coarse.0_5 <- raster("coarse_0_5.preds.mean.tif")
coarse.5_15 <- raster("coarse_5_15.preds.mean.tif")
coarse.15_30 <- raster("coarse_15_30.preds.mean.tif")
coarse.30_60 <- raster("coarse_30_60.preds.mean.tif")
coarse.60_100 <- raster("coarse_60_100.preds.mean.tif")
coarse.100_200 <- raster("coarse_100_200.preds.mean.tif")

preds.s <- stack(clay.0_5,clay.5_15,clay.15_30, clay.30_60,clay.60_100,clay.100_200,
                 silt.0_5,silt.5_15,silt.15_30,silt.30_60,silt.60_100,silt.100_200,
                 sand.0_5,sand.5_15,sand.15_30,sand.30_60,sand.60_100,sand.100_200,
                 coarse.0_5,coarse.5_15,coarse.15_30,coarse.30_60,coarse.60_100, coarse.100_200)

rm(clay.0_5,clay.5_15,clay.15_30, clay.30_60,clay.60_100,clay.100_200,
   silt.0_5,silt.5_15,silt.15_30,silt.30_60,silt.60_100,silt.100_200,
   sand.0_5,sand.5_15,sand.15_30,sand.30_60,sand.60_100,sand.100_200,
   coarse.0_5,coarse.5_15,coarse.15_30,coarse.30_60,coarse.60_100,coarse.100_200)

dir.create(paste0(HomeDir,"Clean_Output/9-Validation.preds.RMQS"))
setwd(paste0(HomeDir,"Clean_Output/9-Validation.preds.RMQS"))
save.image("9-Validation.RMQS.RData")

#### Load also the variance of sand, silt, sand, and coarse elements
setwd(paste0(HomeDir,"Clean_Output/7.1-AWC_Taylor_0_5"))
clay.var.0_5 <- raster("var_clay_0_5.tif")
silt.var.0_5 <- raster("var_silt_0_5.tif")
sand.var.0_5 <- raster("var_sand_0_5.tif")

setwd(paste0(HomeDir,"Clean_Output/7.2-AWC_Taylor_5_15"))
clay.var.5_15 <- raster("var_clay_5_15.tif")
silt.var.5_15 <- raster("var_silt_5_15.tif")
sand.var.5_15 <- raster("var_sand_5_15.tif")

setwd(paste0(HomeDir,"Clean_Output/7.3-AWC_Taylor_15_30"))
clay.var.15_30 <- raster("var_clay_15_30.tif")
silt.var.15_30 <- raster("var_silt_15_30.tif")
sand.var.15_30 <- raster("var_sand_15_30.tif")

setwd(paste0(HomeDir,"Clean_Output/7.4-AWC_Taylor_30_60"))
clay.var.30_60 <- raster("var_clay_30_60.tif")
silt.var.30_60 <- raster("var_silt_30_60.tif")
sand.var.30_60 <- raster("var_sand_30_60.tif")

setwd(paste0(HomeDir,"Clean_Output/7.5-AWC_Taylor_60_100"))
clay.var.60_100 <- raster("var_clay_60_100.tif")
silt.var.60_100 <- raster("var_silt_60_100.tif")
sand.var.60_100 <- raster("var_sand_60_100.tif")

setwd(paste0(HomeDir,"Clean_Output/7.6-AWC_Taylor_100_200"))
clay.var.100_200 <- raster("var_clay_100_200.tif")
silt.var.100_200 <- raster("var_silt_100_200.tif")
sand.var.100_200 <- raster("var_sand_100_200.tif")

setwd(paste0(HomeDir,"Clean_Output/8-CoarseElementsMaps"))
coarse.var.0_5 <- raster("coarse_0_5.preds.sd.tif")
coarse.var.5_15 <- raster("coarse_5_15.preds.sd.tif")
coarse.var.15_30 <- raster("coarse_15_30.preds.sd.tif")
coarse.var.30_60 <- raster("coarse_30_60.preds.sd.tif")
coarse.var.60_100 <- raster("coarse_60_100.preds.sd.tif")
coarse.var.100_200 <- raster("coarse_100_200.preds.sd.tif")


### Add to the stack
preds.s <- stack(preds.s, clay.var.0_5, clay.var.5_15, clay.var.15_30, clay.var.30_60, clay.var.60_100, clay.var.100_200, 
                 silt.var.0_5, silt.var.5_15, silt.var.15_30, silt.var.30_60, silt.var.60_100, silt.var.100_200,
                 sand.var.0_5, sand.var.5_15, sand.var.15_30, sand.var.30_60, sand.var.60_100, sand.var.100_200,
                 coarse.var.0_5, coarse.var.5_15, coarse.var.15_30, coarse.var.30_60, coarse.var.60_100, coarse.var.100_200)

### rm from the workspace
rm(clay.var.0_5, clay.var.5_15, clay.var.15_30, clay.var.30_60, clay.var.60_100, clay.var.100_200, 
   silt.var.0_5, silt.var.5_15, silt.var.15_30, silt.var.30_60, silt.var.60_100, silt.var.100_200,
   sand.var.0_5, sand.var.5_15, sand.var.15_30, sand.var.30_60, sand.var.60_100, sand.var.100_200,
   coarse.var.0_5, coarse.var.5_15, coarse.var.15_30, coarse.var.30_60, coarse.var.60_100, coarse.var.100_200)

# #### Attach land use (clc_EXPERT) ---------------------------------------
occupation <- read.csv(paste0(HomeDir,"Input/RMQS_sites/occupations.csv"), sep=";", header=TRUE,na.strings = "")
id_profil_all_sites.rmqs <- read.csv(paste0(HomeDir,"Input/RMQS_sites/id_profil_tous_sites_rmqs.csv"), sep=";", header=TRUE,na.strings = "")
id_profil_sites_foret.rmqs <- read.csv(paste0(HomeDir,"Input/RMQS_sites/id_profil_sites_forestiers.csv"), sep=";", header=TRUE,na.strings = "")

# 
id_profil_all_sites.rmqs <- merge(id_profil_all_sites.rmqs, occupation, by.x="id_site", by.y="id_site_donesol")
# 
### fill each other's columns...
RMQS.data[is.na(RMQS.data$id_site.x) & !is.na(RMQS.data$id_site.y), ]$id_site.x <-
    RMQS.data[is.na(RMQS.data$id_site.x) & !is.na(RMQS.data$id_site.y), ]$id_site.y
    
RMQS.data[is.na(RMQS.data$id_site.y) & !is.na(RMQS.data$id_site.x), ]$id_site.y <-
    RMQS.data[is.na(RMQS.data$id_site.y) & !is.na(RMQS.data$id_site.x), ]$id_site.x

RMQS.data[is.na(RMQS.data$id_profil.x) & !is.na(RMQS.data$id_profil.y), ]$id_profil.x <-
    RMQS.data[is.na(RMQS.data$id_profil.x) & !is.na(RMQS.data$id_profil.y), ]$id_profil.y

RMQS.data[is.na(RMQS.data$id_profil.y) & !is.na(RMQS.data$id_profil.x), ]$id_profil.y <-
    RMQS.data[is.na(RMQS.data$id_profil.y) & !is.na(RMQS.data$id_profil.x), ]$id_profil.x
# 
# ## merge land use
RMQS.data <- merge(RMQS.data, id_profil_all_sites.rmqs, by.x="id_profil.x", by.y="id_profil", all.x=TRUE, all.y=FALSE)
# 
summary(RMQS.data$code_occup3_dexpert)
RMQS.data$land_use <- NA
# 
RMQS.data$land_use <- ifelse(RMQS.data$code_occup3_dexpert <= 199, "urban", 
                             ifelse(RMQS.data$code_occup3_dexpert <= 299, "agricultural",
                                    ifelse(RMQS.data$code_occup3_dexpert <= 399, "forest",
                                           ifelse(RMQS.data$code_occup3_dexpert <= 499, "wetlands",
                                                  "water bodies"))))
RMQS.data$land_use <- as.factor(RMQS.data$land_use )
summary(RMQS.data$land_use)

### Just see if they are in forest
RMQS.data$forest <- NA
RMQS.data$forest <- ifelse(RMQS.data$id_profil.x %in% unique(id_profil_sites_foret.rmqs$id_profil), "forest", "agriculture")
RMQS.data$forest <- as.factor(RMQS.data$forest)
summary(RMQS.data$forest)

### the results seem wrong....

# ### Spatial extraction --------------------------------------------------

### extract predictions at RMQS sites
### copy
RMQS.data.sp <- RMQS.data
RMQS.data.sp[is.na(RMQS.data.sp$x_reel),]$x_reel <- RMQS.data.sp[is.na(RMQS.data.sp$x_reel),]$x ### fill missing coordinates
RMQS.data.sp[is.na(RMQS.data.sp$y_reel),]$y_reel <- RMQS.data.sp[is.na(RMQS.data.sp$y_reel),]$y

### Transform to spatial and assign CRS
coordinates(RMQS.data.sp) <- ~x_reel +y_reel
proj4string(RMQS.data.sp) <- " +proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
### extract predictions from raster stack
preds.df <- extract(preds.s, RMQS.data.sp,method="simple", df=TRUE)

### bind
RMQS.data <- cbind(RMQS.data, preds.df)

# ### Calculate average predictions for clay ------------------------------

### Calculate average predictions for clay

### Clay
RMQS.data.clay <- RMQS.data[!is.na(RMQS.data$clay),]

### Eliminate the horizons that are deeper than 2 m
RMQS.data.clay <- RMQS.data.clay[RMQS.data.clay$profondeur_hz_sup<200,]

### Define the arguments for this loop
### Working with clay
GSM <- c("argile.preds.RK.0_5","argile.preds.RK.5_15", "argile.preds.RK.15_30", "argile.preds.RK.30_60","argile.preds.RK.60_100","argile.preds.RK.100_200")

### The name for the new valiable
# var.pred <- "clay.pred"

### Create empty variable
RMQS.data.clay$clay.pred <- NA
    
for(Hrow in 1:nrow(RMQS.data.clay)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.clay[Hrow,]$profondeur_hz_inf <- ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_inf > 200, 
                                                      yes= 200,
                                                      no= RMQS.data.clay[Hrow,]$profondeur_hz_inf )
    
    t <- RMQS.data.clay[Hrow,]$profondeur_hz_inf - RMQS.data.clay[Hrow,]$profondeur_hz_sup 
    
    GSM.1 <- RMQS.data.clay[Hrow, GSM[[1]]] 
    GSM.2 <- RMQS.data.clay[Hrow, GSM[[2]]]
    GSM.3 <- RMQS.data.clay[Hrow, GSM[[3]]]
    GSM.4 <- RMQS.data.clay[Hrow, GSM[[4]]]
    GSM.5 <- RMQS.data.clay[Hrow, GSM[[5]]]
    GSM.6 <- RMQS.data.clay[Hrow, GSM[[6]]]
    
    lim.GSM.sup <- c(0,5,15,30,60,100)
    lim.GSM.inf <- c(5,15,30,60,100,200)
    GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
    
    #### 1. the RMQS layer falls within a GSM layer
    
    ## 1.1 In the 0-5 cm
    
    if(RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[1] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[1]) {  
        ## we are lucky, and our RMQS horizon falls in the first layer
        RMQS.data.clay[Hrow,]$clay.pred <- GSM.1
        
    ## 1.2 In the 5-15 cm    
    } else if(RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[2] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[2]){
        
        RMQS.data.clay[Hrow,]$clay.pred <- GSM.2
    
    ## 1.3 In the 15-30 cm
    } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[3] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[3]){
            
        RMQS.data.clay[Hrow,]$clay.pred <- GSM.3
            
    ## 1.4 In the 30-60 cm
    } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[4] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[4]){
        
        RMQS.data.clay[Hrow,]$clay.pred <- GSM.4
        
    ## 1.5 In the 60-100 cm
    } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[5] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[5]){
        
        RMQS.data.clay[Hrow,]$clay.pred <- GSM.5
        
    ## 1.6 In the 100-200 cm
    } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[6] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[6]){
        
        RMQS.data.clay[Hrow,]$clay.pred <- GSM.6
        
    } else {
        
        ### 2. the RMQS layer covers several GSM horizons
        
        ### Define weights for each GSM layer
        weights <- c(rep(0,6))
        
        for(i in 1:6){
            
            weights[i] <- ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                     RMQS.data.clay[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                     RMQS.data.clay[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] &
                                     RMQS.data.clay[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i],
                                 yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                 no= ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                                RMQS.data.clay[Hrow,]$profondeur_hz_sup <= lim.GSM.inf[i] &
                                                RMQS.data.clay[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] & 
                                                RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[i],
                                            yes= (RMQS.data.clay[Hrow,]$profondeur_hz_inf - lim.GSM.sup[i])/t,
                                            no= ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[i] &
                                                           RMQS.data.clay[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                                           RMQS.data.clay[Hrow,]$profondeur_hz_inf >= lim.GSM.sup[i] & 
                                                           RMQS.data.clay[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i], 
                                                       yes= (lim.GSM.inf[i] - RMQS.data.clay[Hrow,]$profondeur_hz_sup)/t,
                                                       no= 0 )))
            
            ### check that the sum of weights is 1
            
        } 
        
        if (all.equal(sum(weights),1,tolerance=0.02)){
            
            RMQS.data.clay[Hrow,]$clay.pred <- sum(weights*GSM.preds)
            
        } else {
            
            print(paste (Hrow, "error with weights for", RMQS.data.clay[Hrow,]$id_tot1 ))
        }
    }
}

rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t, weights)

### check those horizon problem ### Fix one by one

### REMEMBER THAT CLAY, SILT, SAND, ARE IN %, NOT IN G/KG
plot(RMQS.data.clay$clay, RMQS.data.clay$clay.pred*10)
abline(0,1,col="blue")

### Just in case there are sites without data or predictions
RMQS.data.clay.full <- RMQS.data.clay[!is.na(RMQS.data.clay$clay.pred), ]  ### There were around 300 horizons with NA
### Because there are NA 
error.RMQS.data.clay <- RMQS.data.clay[is.na(RMQS.data.clay$clay.pred) & !is.na(RMQS.data.clay$argile.preds.RK.0_5), ]

### Check code
validation.clay <- goof(observed = RMQS.data.clay.full$clay,
     predicted = RMQS.data.clay.full$clay.pred*10, plot.it = TRUE, type = "DSM")

par(mfrow=c(1,1), las=1)
plot(RMQS.data.clay.full$clay,RMQS.data.clay.full$clay.pred*10,
     pch=19,cex=0.7,col="darkgoldenrod1", xlab=" ",   ylab=" ", main="Clay");abline(0,1)

RMQS.data.clay.full$clay.pred.hor <- RMQS.data.clay.full$clay.pred*10
clay.plot <- ggplot(RMQS.data.clay.full, aes(clay.pred.hor, clay))
clay.plot + geom_point(alpha = 1/6,colour = "black", size = 2) + theme_bw()+
    labs(x = "Predicted clay (gkg-1)", y ="Observed clay (gkg-1)")+
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))+
     geom_abline(slope=1, intercept = 0,size=1, colour="blue")
    
save.image("9-Validation.RMQS.RData")

# ### Calculate and validate, average predictions for silt ----------------

### Calculate average predictions for silt

### Silt
RMQS.data.silt <- RMQS.data[!is.na(RMQS.data$silt),]

### Eliminate the horizons that are deeper than 2 m
RMQS.data.silt <- RMQS.data.silt[RMQS.data.silt$profondeur_hz_sup<200,]

### Define the arguments for this loop
### Working with silt
GSM <- c("limon.preds.RK.0_5","limon.preds.RK.5_15", "limon.preds.RK.15_30", "limon.preds.RK.30_60","limon.preds.RK.60_100","limon.preds.RK.100_200")
### The name for the new valiable
# var.pred <- "silt.pred"

### Create empty variable
RMQS.data.silt$silt.pred <- NA

for(Hrow in 1:nrow(RMQS.data.silt)){
    
    print(Hrow)
    
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.silt[Hrow,]$profondeur_hz_inf <- ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_inf > 200, 
                                                      yes= 200,
                                                      no= RMQS.data.silt[Hrow,]$profondeur_hz_inf )
        
        t <- RMQS.data.silt[Hrow,]$profondeur_hz_inf - RMQS.data.silt[Hrow,]$profondeur_hz_sup 
        
        GSM.1 <- RMQS.data.silt[Hrow, GSM[[1]]] 
        GSM.2 <- RMQS.data.silt[Hrow, GSM[[2]]]
        GSM.3 <- RMQS.data.silt[Hrow, GSM[[3]]]
        GSM.4 <- RMQS.data.silt[Hrow, GSM[[4]]]
        GSM.5 <- RMQS.data.silt[Hrow, GSM[[5]]]
        GSM.6 <- RMQS.data.silt[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[1] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            RMQS.data.silt[Hrow,]$silt.pred <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[2] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[2]){
            
            RMQS.data.silt[Hrow,]$silt.pred <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[3] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[3]){
            
            RMQS.data.silt[Hrow,]$silt.pred <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[4] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[4]){
            
            RMQS.data.silt[Hrow,]$silt.pred <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[5] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[5]){
            
            RMQS.data.silt[Hrow,]$silt.pred <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[6] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[6]){
            
            RMQS.data.silt[Hrow,]$silt.pred <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            weights <- c(rep(0,6))
            
            for(i in 1:6){
                
                weights[i] <- ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                         RMQS.data.silt[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                         RMQS.data.silt[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] &
                                         RMQS.data.silt[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i],
                                     yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                     no= ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                                    RMQS.data.silt[Hrow,]$profondeur_hz_sup <= lim.GSM.inf[i] &
                                                    RMQS.data.silt[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] & 
                                                    RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[i],
                                                yes= (RMQS.data.silt[Hrow,]$profondeur_hz_inf - lim.GSM.sup[i])/t,
                                                no= ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[i] &
                                                               RMQS.data.silt[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                                               RMQS.data.silt[Hrow,]$profondeur_hz_inf >= lim.GSM.sup[i] & 
                                                               RMQS.data.silt[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i], 
                                                           yes= (lim.GSM.inf[i] - RMQS.data.silt[Hrow,]$profondeur_hz_sup)/t,
                                                           no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(weights),1,tolerance=0.02)){
                
                RMQS.data.silt[Hrow,]$silt.pred <- sum(weights*GSM.preds)
                
            } else {
                
                print(paste (Hrow, "error with weights for", RMQS.data.silt[Hrow,]$id_tot1 ))
            }
        }
}

rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t, weights)

### check those horizon problem ### Fix one by one

### REMEMBER THAT CLAY, SILT, SAND, ARE IN %, NOT IN G/KG
plot(RMQS.data.silt$silt, RMQS.data.silt$silt.pred*10)
abline(0,1,col="blue")

plot(RMQS.data.silt$silt,RMQS.data.silt$silt.pred*10,
     pch=19,cex=0.7,col="antiquewhite3", xlab=" ",   ylab=" ", main="Silt");abline(0,1)

### Just in case there are sites without data or predictions
RMQS.data.silt.full <- RMQS.data.silt[!is.na(RMQS.data.silt$silt.pred), ]  ### There were around 300 horizons with NA
### Because there are NA 
error.RMQS.data.silt <- RMQS.data.silt[is.na(RMQS.data.silt$silt.pred) & !is.na(RMQS.data.silt$limon.preds.RK.0_5), ]

### Check code
validation.silt <- goof(observed = RMQS.data.silt.full$silt,
                        predicted = RMQS.data.silt.full$silt.pred*10, plot.it = TRUE, type = "DSM")

par(mfrow=c(1,1), las=1)
plot(RMQS.data.silt.full$silt,RMQS.data.silt.full$silt.pred*10,
     pch=19,cex=0.7,col="darkgoldenrod1", xlab=" ",   ylab=" ", main="silt");abline(0,1)

RMQS.data.silt.full$silt.pred.hor <- RMQS.data.silt.full$silt.pred*10
silt.plot <- ggplot(RMQS.data.silt.full, aes(silt.pred.hor, silt))
silt.plot + geom_point(alpha = 1/6,colour = "black", size = 2) + theme_bw()+
    labs(x = "Predicted silt (gkg-1)", y ="Observed silt (gkg-1)")+
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))+
    geom_abline(slope=1, intercept = 0,size=1, colour="blue")



# ### Calculate and validate, average predictions for sand ----------------

### Calculate average predictions for sand

### Sand
RMQS.data.sand <- RMQS.data[!is.na(RMQS.data$sand),]

### Eliminate the horizons that are deeper than 2 m
RMQS.data.sand <- RMQS.data.sand[RMQS.data.sand$profondeur_hz_sup<200,]

### Define the arguments for this loop
### Working with sand
GSM <- c("sable.preds.RK.0_5","sable.preds.RK.5_15", "sable.preds.RK.15_30", "sable.preds.RK.30_60","sable.preds.RK.60_100","sable.preds.RK.100_200")
### The name for the new valiable
# var.pred <- "sand.pred"

### Create empty variable
RMQS.data.sand$sand.pred <- NA

for(Hrow in 1:nrow(RMQS.data.sand)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.sand[Hrow,]$profondeur_hz_inf <- ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_inf > 200, 
                                                      yes= 200,
                                                      no= RMQS.data.sand[Hrow,]$profondeur_hz_inf )
        
        t <- RMQS.data.sand[Hrow,]$profondeur_hz_inf - RMQS.data.sand[Hrow,]$profondeur_hz_sup 
        
        GSM.1 <- RMQS.data.sand[Hrow, GSM[[1]]] 
        GSM.2 <- RMQS.data.sand[Hrow, GSM[[2]]]
        GSM.3 <- RMQS.data.sand[Hrow, GSM[[3]]]
        GSM.4 <- RMQS.data.sand[Hrow, GSM[[4]]]
        GSM.5 <- RMQS.data.sand[Hrow, GSM[[5]]]
        GSM.6 <- RMQS.data.sand[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[1] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            RMQS.data.sand[Hrow,]$sand.pred <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[2] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[2]){
            
            RMQS.data.sand[Hrow,]$sand.pred <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[3] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[3]){
            
            RMQS.data.sand[Hrow,]$sand.pred <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[4] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[4]){
            
            RMQS.data.sand[Hrow,]$sand.pred <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[5] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[5]){
            
            RMQS.data.sand[Hrow,]$sand.pred <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[6] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[6]){
            
            RMQS.data.sand[Hrow,]$sand.pred <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            weights <- c(rep(0,6))
            
            for(i in 1:6){
                
                weights[i] <- ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                         RMQS.data.sand[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                         RMQS.data.sand[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] &
                                         RMQS.data.sand[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i],
                                     yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                     no= ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                                    RMQS.data.sand[Hrow,]$profondeur_hz_sup <= lim.GSM.inf[i] &
                                                    RMQS.data.sand[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] & 
                                                    RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[i],
                                                yes= (RMQS.data.sand[Hrow,]$profondeur_hz_inf - lim.GSM.sup[i])/t,
                                                no= ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[i] &
                                                               RMQS.data.sand[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                                               RMQS.data.sand[Hrow,]$profondeur_hz_inf >= lim.GSM.sup[i] & 
                                                               RMQS.data.sand[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i], 
                                                           yes= (lim.GSM.inf[i] - RMQS.data.sand[Hrow,]$profondeur_hz_sup)/t,
                                                           no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(weights),1,tolerance=0.02)){
                
                RMQS.data.sand[Hrow,]$sand.pred <- sum(weights*GSM.preds)
                
            } else {
                
                print(paste (Hrow, "error with weights for", RMQS.data.sand[Hrow,]$id_tot1 ))
            }
        }
    }

rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t, weights)

### check those horizon problem ### Fix one by one

### REMEMBER THAT CLAY, SILT, SAND, ARE IN %, NOT IN G/KG
plot(RMQS.data.sand$sand, RMQS.data.sand$sand.pred*10)
abline(0,1,col="blue")

plot(RMQS.data.sand$sand,RMQS.data.sand$sand.pred*10,
     pch=19,cex=0.7,col="burlywood3", xlab=" ",   ylab=" ", main="Sand");abline(0,1)

### Just in case there are sites without data or predictions
RMQS.data.sand.full <- RMQS.data.sand[!is.na(RMQS.data.sand$sand.pred), ]  ### There were around 300 horizons with NA
### Because there are NA 
error.RMQS.data.sand <- RMQS.data.sand[is.na(RMQS.data.sand$sand.pred) & !is.na(RMQS.data.sand$sable.preds.RK.0_5), ]

### Check code
validation.sand <- goof(observed = RMQS.data.sand.full$sand,
                        predicted = RMQS.data.sand.full$sand.pred*10, plot.it = TRUE, type = "DSM")

par(mfrow=c(1,1), las=1)
plot(RMQS.data.sand.full$sand,RMQS.data.sand.full$sand.pred*10,
     pch=19,cex=0.7,col="darkgoldenrod1", xlab=" ",   ylab=" ", main="sand");abline(0,1)

RMQS.data.sand.full$sand.pred.hor <- RMQS.data.sand.full$sand.pred*10
sand.plot <- ggplot(RMQS.data.sand.full, aes(sand.pred.hor, sand))
sand.plot + geom_point(alpha = 1/6,colour = "black", size = 2) + theme_bw()+
    labs(x = "Predicted sand (gkg-1)", y ="Observed sand (gkg-1)")+
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))+
    geom_abline(slope=1, intercept = 0,size=1, colour="blue")

save.image("9-Validation.RMQS.RData")

# ### Calculate and validate, average predictions for coarse ----------------

### Calculate average predictions for coarse

### Coarse elements
RMQS.data.coarse <- RMQS.data[!is.na(RMQS.data$abondance_eg),]

### Eliminate the horizons that are deeper than 2 m
RMQS.data.coarse <- RMQS.data.coarse[RMQS.data.coarse$prof_sup_moy<200,]

### Define the arguments for this loop
### Working with coarse
GSM <- c("coarse_0_5.preds.mean" , "coarse_5_15.preds.mean", "coarse_15_30.preds.mean","coarse_30_60.preds.mean", "coarse_60_100.preds.mean",  "coarse_100_200.preds.mean")
### The name for the new valiable
# var.pred <- "coarse.pred"

### Create empty variable
RMQS.data.coarse$abondance_eg.pred <- NA


### check whether using 

### profondeur_hz_sup or prof_sup_moy
### profondeur_hz_inf or prof_inf_moy


dim(RMQS.data.coarse[is.na(RMQS.data.coarse$profondeur_hz_sup),])
dim(RMQS.data.coarse[is.na(RMQS.data.coarse$prof_sup_moy),])
dim(RMQS.data.coarse[is.na(RMQS.data.coarse$profondeur_hz_inf),])
dim(RMQS.data.coarse[is.na(RMQS.data.coarse$prof_inf_moy),])
dim(RMQS.data.coarse[!is.na(RMQS.data.coarse$prof_inf_moy) & !is.na(RMQS.data.coarse$profondeur_hz_inf),])

RMQS.data.coarse[!is.na(RMQS.data.coarse$prof_inf_moy) & !is.na(RMQS.data.coarse$profondeur_hz_inf), c("prof_inf_moy","profondeur_hz_inf" )]
summary(RMQS.data.coarse[!is.na(RMQS.data.coarse$prof_inf_moy) & !is.na(RMQS.data.coarse$profondeur_hz_inf),]$profondeur_hz_sup -
         RMQS.data.coarse[!is.na(RMQS.data.coarse$prof_inf_moy) & !is.na(RMQS.data.coarse$profondeur_hz_inf),]$prof_sup_moy)

head(RMQS.data.coarse[is.na(RMQS.data.coarse$profondeur_hz_sup),])

### here, use 
### prof_sup_moy and prof_inf_moy


### check some errors

tic <- Sys.time()
for(Hrow in 1:nrow(RMQS.data.coarse)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.coarse[Hrow,]$prof_inf_moy <- ifelse(test=RMQS.data.coarse[Hrow,]$prof_inf_moy > 200, 
                                                        yes= 200,
                                                        no= RMQS.data.coarse[Hrow,]$prof_inf_moy )
        
        t <- RMQS.data.coarse[Hrow,]$prof_inf_moy - RMQS.data.coarse[Hrow,]$prof_sup_moy 
       
        GSM.1 <- RMQS.data.coarse[Hrow, GSM[[1]]] 
        GSM.2 <- RMQS.data.coarse[Hrow, GSM[[2]]]
        GSM.3 <- RMQS.data.coarse[Hrow, GSM[[3]]]
        GSM.4 <- RMQS.data.coarse[Hrow, GSM[[4]]]
        GSM.5 <- RMQS.data.coarse[Hrow, GSM[[5]]]
        GSM.6 <- RMQS.data.coarse[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[1] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            RMQS.data.coarse[Hrow,]$abondance_eg.pred <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[2] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[2]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[3] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[3]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[4] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[4]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[5] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[5]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[6] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[6]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            weights <- c(rep(0,6))
            
            for(i in 1:6){
                
                weights[i] <- ifelse(test=RMQS.data.coarse[Hrow,]$prof_sup_moy <= lim.GSM.sup[i] &
                                         RMQS.data.coarse[Hrow,]$prof_sup_moy < lim.GSM.inf[i] &
                                         RMQS.data.coarse[Hrow,]$prof_inf_moy > lim.GSM.sup[i] &
                                         RMQS.data.coarse[Hrow,]$prof_inf_moy >= lim.GSM.inf[i],
                                     yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                     no= ifelse(test=RMQS.data.coarse[Hrow,]$prof_sup_moy <= lim.GSM.sup[i] &
                                                    RMQS.data.coarse[Hrow,]$prof_sup_moy <= lim.GSM.inf[i] &
                                                    RMQS.data.coarse[Hrow,]$prof_inf_moy > lim.GSM.sup[i] & 
                                                    RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[i],
                                                yes= (RMQS.data.coarse[Hrow,]$prof_inf_moy - lim.GSM.sup[i])/t,
                                                no= ifelse(test=RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[i] &
                                                               RMQS.data.coarse[Hrow,]$prof_sup_moy < lim.GSM.inf[i] &
                                                               RMQS.data.coarse[Hrow,]$prof_inf_moy >= lim.GSM.sup[i] & 
                                                               RMQS.data.coarse[Hrow,]$prof_inf_moy >= lim.GSM.inf[i], 
                                                           yes= (lim.GSM.inf[i] - RMQS.data.coarse[Hrow,]$prof_sup_moy)/t,
                                                           no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(weights),1,tolerance=0.02)){
                
                RMQS.data.coarse[Hrow,]$abondance_eg.pred <- sum(weights*GSM.preds)
                
            } else {
                
                print(paste (Hrow, "error with weights for", RMQS.data.coarse[Hrow,]$id_tot1 ))
            }
        }
    }
tac <- Sys.time()
tac-tic

rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t, weights)

### check those horizon problem ### Fix one by one

### REMEMBER THAT COARSE ELEMENTS IS IN %, so no need to multiply or divide
plot(RMQS.data.coarse$abondance_eg, RMQS.data.coarse$abondance_eg.pred)
abline(0,1,col="blue")

plot(RMQS.data.coarse$abondance_eg,RMQS.data.coarse$abondance_eg.pred,
     pch=19,cex=0.7,col="lightpink2", xlab=" ",   ylab=" ", main="Coarse elements");abline(0,1)

### Just in case there are sites without data or predictions
RMQS.data.coarse.full <- RMQS.data.coarse[!is.na(RMQS.data.coarse$abondance_eg.pred), ]  
### Because there are NA 
error.RMQS.data.coarse <- RMQS.data.coarse[is.na(RMQS.data.coarse$abondance_eg.pred) & !is.na(RMQS.data.coarse$sable.preds.RK.0_5), ]
rm(error.RMQS.data.coarse, error.RMQS.data.sand, error.RMQS.data.silt, error.RMQS.data.clay)

validation.coarse <- goof(observed = RMQS.data.coarse.full$abondance_eg,
                        predicted = RMQS.data.coarse.full$abondance_eg.pred, plot.it = TRUE, type = "DSM")

par(mfrow=c(1,1), las=1)
plot(RMQS.data.coarse.full$abondance_eg,RMQS.data.coarse.full$abondance_eg.pred,
     pch=19,cex=0.7,col="darkgoldenrod1", xlab=" ",   ylab=" ", main="coarse elements");abline(0,1)

coarse.plot <- ggplot(RMQS.data.coarse.full, aes(abondance_eg.pred, abondance_eg))
coarse.plot + geom_point(alpha = 1/6,colour = "black", size = 2) + theme_bw()+
    labs(x = "Predicted coarse elements (%)", y ="Observed coarse elements (%)")+
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))+
    geom_abline(slope=1, intercept = 0,size=1, colour="blue")

#dir.create("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/8.Validation.preds.RMQS")
save.image("9-Validation.RMQS.RData")

### Clean workspace
rm(tac, tic, GSM.preds, preds.df, rmqs.eg.profil, rmqs.gr.profil)

# ### Calculate VARIANCE of weighed average predictions -------------------

### First, I need to calculate the correlation between the GSM layer predictions (puff!)
tic <- Sys.time()
cor.clay1.2 <- cor.test(x=getValues(preds.s[[1]]), y=getValues(preds.s[[2]]), na.rm=TRUE, method="pearson")$estimate
cor.clay1.3 <- cor.test(x=getValues(preds.s[[1]]), y=getValues(preds.s[[3]]), na.rm=TRUE, method="pearson")$estimate
cor.clay1.4 <- cor.test(x=getValues(preds.s[[1]]), y=getValues(preds.s[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.clay1.5 <- cor.test(x=getValues(preds.s[[1]]), y=getValues(preds.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.clay1.6 <- cor.test(x=getValues(preds.s[[1]]), y=getValues(preds.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.clay2.3 <- cor.test(x=getValues(preds.s[[2]]), y=getValues(preds.s[[3]]), na.rm=TRUE, method="pearson")$estimate
cor.clay2.4 <- cor.test(x=getValues(preds.s[[2]]), y=getValues(preds.s[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.clay2.5 <- cor.test(x=getValues(preds.s[[2]]), y=getValues(preds.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.clay2.6 <- cor.test(x=getValues(preds.s[[2]]), y=getValues(preds.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.clay3.4 <- cor.test(x=getValues(preds.s[[3]]), y=getValues(preds.s[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.clay3.5 <- cor.test(x=getValues(preds.s[[3]]), y=getValues(preds.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.clay3.6 <- cor.test(x=getValues(preds.s[[3]]), y=getValues(preds.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.clay4.5 <- cor.test(x=getValues(preds.s[[4]]), y=getValues(preds.s[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.clay4.6 <- cor.test(x=getValues(preds.s[[4]]), y=getValues(preds.s[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.clay5.6 <- cor.test(x=getValues(preds.s[[5]]), y=getValues(preds.s[[6]]), na.rm=TRUE, method="pearson")$estimate

cor.silt1.2 <- cor.test(x=getValues(preds.s[[7]]), y=getValues(preds.s[[8]]), na.rm=TRUE, method="pearson")$estimate
cor.silt1.3 <- cor.test(x=getValues(preds.s[[7]]), y=getValues(preds.s[[9]]), na.rm=TRUE, method="pearson")$estimate
cor.silt1.4 <- cor.test(x=getValues(preds.s[[7]]), y=getValues(preds.s[[10]]), na.rm=TRUE, method="pearson")$estimate
cor.silt1.5 <- cor.test(x=getValues(preds.s[[7]]), y=getValues(preds.s[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.silt1.6 <- cor.test(x=getValues(preds.s[[7]]), y=getValues(preds.s[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.silt2.3 <- cor.test(x=getValues(preds.s[[8]]), y=getValues(preds.s[[9]]), na.rm=TRUE, method="pearson")$estimate
cor.silt2.4 <- cor.test(x=getValues(preds.s[[8]]), y=getValues(preds.s[[10]]), na.rm=TRUE, method="pearson")$estimate
cor.silt2.5 <- cor.test(x=getValues(preds.s[[8]]), y=getValues(preds.s[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.silt2.6 <- cor.test(x=getValues(preds.s[[8]]), y=getValues(preds.s[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.silt3.4 <- cor.test(x=getValues(preds.s[[9]]), y=getValues(preds.s[[10]]), na.rm=TRUE, method="pearson")$estimate
cor.silt3.5 <- cor.test(x=getValues(preds.s[[9]]), y=getValues(preds.s[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.silt3.6 <- cor.test(x=getValues(preds.s[[9]]), y=getValues(preds.s[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.silt4.5 <- cor.test(x=getValues(preds.s[[10]]), y=getValues(preds.s[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.silt4.6 <- cor.test(x=getValues(preds.s[[10]]), y=getValues(preds.s[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.silt5.6 <- cor.test(x=getValues(preds.s[[11]]), y=getValues(preds.s[[12]]), na.rm=TRUE, method="pearson")$estimate

cor.sand1.2 <- cor.test(x=getValues(preds.s[[13]]), y=getValues(preds.s[[14]]), na.rm=TRUE, method="pearson")$estimate
cor.sand1.3 <- cor.test(x=getValues(preds.s[[13]]), y=getValues(preds.s[[15]]), na.rm=TRUE, method="pearson")$estimate
cor.sand1.4 <- cor.test(x=getValues(preds.s[[13]]), y=getValues(preds.s[[16]]), na.rm=TRUE, method="pearson")$estimate
cor.sand1.5 <- cor.test(x=getValues(preds.s[[13]]), y=getValues(preds.s[[17]]), na.rm=TRUE, method="pearson")$estimate
cor.sand1.6 <- cor.test(x=getValues(preds.s[[13]]), y=getValues(preds.s[[18]]), na.rm=TRUE, method="pearson")$estimate
cor.sand2.3 <- cor.test(x=getValues(preds.s[[14]]), y=getValues(preds.s[[15]]), na.rm=TRUE, method="pearson")$estimate
cor.sand2.4 <- cor.test(x=getValues(preds.s[[14]]), y=getValues(preds.s[[16]]), na.rm=TRUE, method="pearson")$estimate
cor.sand2.5 <- cor.test(x=getValues(preds.s[[14]]), y=getValues(preds.s[[17]]), na.rm=TRUE, method="pearson")$estimate
cor.sand2.6 <- cor.test(x=getValues(preds.s[[14]]), y=getValues(preds.s[[18]]), na.rm=TRUE, method="pearson")$estimate
cor.sand3.4 <- cor.test(x=getValues(preds.s[[15]]), y=getValues(preds.s[[16]]), na.rm=TRUE, method="pearson")$estimate
cor.sand3.5 <- cor.test(x=getValues(preds.s[[15]]), y=getValues(preds.s[[17]]), na.rm=TRUE, method="pearson")$estimate
cor.sand3.6 <- cor.test(x=getValues(preds.s[[15]]), y=getValues(preds.s[[18]]), na.rm=TRUE, method="pearson")$estimate
cor.sand4.5 <- cor.test(x=getValues(preds.s[[16]]), y=getValues(preds.s[[17]]), na.rm=TRUE, method="pearson")$estimate
cor.sand4.6 <- cor.test(x=getValues(preds.s[[16]]), y=getValues(preds.s[[18]]), na.rm=TRUE, method="pearson")$estimate
cor.sand5.6 <- cor.test(x=getValues(preds.s[[17]]), y=getValues(preds.s[[18]]), na.rm=TRUE, method="pearson")$estimate

cor.coarse1.2 <- cor.test(x=getValues(preds.s[[19]]), y=getValues(preds.s[[20]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse1.3 <- cor.test(x=getValues(preds.s[[19]]), y=getValues(preds.s[[21]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse1.4 <- cor.test(x=getValues(preds.s[[19]]), y=getValues(preds.s[[22]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse1.5 <- cor.test(x=getValues(preds.s[[19]]), y=getValues(preds.s[[23]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse1.6 <- cor.test(x=getValues(preds.s[[19]]), y=getValues(preds.s[[24]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse2.3 <- cor.test(x=getValues(preds.s[[20]]), y=getValues(preds.s[[21]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse2.4 <- cor.test(x=getValues(preds.s[[20]]), y=getValues(preds.s[[22]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse2.5 <- cor.test(x=getValues(preds.s[[20]]), y=getValues(preds.s[[23]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse2.6 <- cor.test(x=getValues(preds.s[[20]]), y=getValues(preds.s[[24]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse3.4 <- cor.test(x=getValues(preds.s[[21]]), y=getValues(preds.s[[22]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse3.5 <- cor.test(x=getValues(preds.s[[21]]), y=getValues(preds.s[[23]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse3.6 <- cor.test(x=getValues(preds.s[[21]]), y=getValues(preds.s[[24]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse4.5 <- cor.test(x=getValues(preds.s[[22]]), y=getValues(preds.s[[23]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse4.6 <- cor.test(x=getValues(preds.s[[22]]), y=getValues(preds.s[[24]]), na.rm=TRUE, method="pearson")$estimate
cor.coarse5.6 <- cor.test(x=getValues(preds.s[[23]]), y=getValues(preds.s[[24]]), na.rm=TRUE, method="pearson")$estimate
tac <- Sys.time()

rm(tac, tic)
save.image("9-Validation.RMQS.RData")

### Now, Calculate the variance of the average, taking into account the correlation between layers......


# ### CLAY VARIANCE (%^2) -------------------------------------------------

### Define the arguments for this loop

### Working with var.clay
GSM <- c("var_clay_0_5","var_clay_5_15", "var_clay_15_30", "var_clay_30_60","var_clay_60_100","var_clay_100_200")

### Create empty variable
RMQS.data.clay$clay.pred.var <- NA

for(Hrow in 1:nrow(RMQS.data.clay)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.clay[Hrow,]$profondeur_hz_inf <- ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_inf > 200, 
                                                      yes= 200,
                                                      no= RMQS.data.clay[Hrow,]$profondeur_hz_inf )
        
        t <- RMQS.data.clay[Hrow,]$profondeur_hz_inf - RMQS.data.clay[Hrow,]$profondeur_hz_sup 
        
        GSM.1 <- RMQS.data.clay[Hrow, GSM[[1]]] 
        GSM.2 <- RMQS.data.clay[Hrow, GSM[[2]]]
        GSM.3 <- RMQS.data.clay[Hrow, GSM[[3]]]
        GSM.4 <- RMQS.data.clay[Hrow, GSM[[4]]]
        GSM.5 <- RMQS.data.clay[Hrow, GSM[[5]]]
        GSM.6 <- RMQS.data.clay[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[1] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            RMQS.data.clay[Hrow,]$clay.pred.var <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[2] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[2]){
            
            RMQS.data.clay[Hrow,]$clay.pred.var <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[3] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[3]){
            
            RMQS.data.clay[Hrow,]$clay.pred.var <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[4] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[4]){
            
            RMQS.data.clay[Hrow,]$clay.pred.var <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[5] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[5]){
            
            RMQS.data.clay[Hrow,]$clay.pred.var <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[6] & RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[6]){
            
            RMQS.data.clay[Hrow,]$clay.pred.var <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            alpha <- c(rep(0,6))
            
            for(i in 1:6){
                
                alpha[i] <- ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                         RMQS.data.clay[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                         RMQS.data.clay[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] &
                                         RMQS.data.clay[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i],
                                     yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                     no= ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                                    RMQS.data.clay[Hrow,]$profondeur_hz_sup <= lim.GSM.inf[i] &
                                                    RMQS.data.clay[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] & 
                                                    RMQS.data.clay[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[i],
                                                yes= (RMQS.data.clay[Hrow,]$profondeur_hz_inf - lim.GSM.sup[i])/t,
                                                no= ifelse(test=RMQS.data.clay[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[i] &
                                                               RMQS.data.clay[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                                               RMQS.data.clay[Hrow,]$profondeur_hz_inf >= lim.GSM.sup[i] & 
                                                               RMQS.data.clay[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i], 
                                                           yes= (lim.GSM.inf[i] - RMQS.data.clay[Hrow,]$profondeur_hz_sup)/t,
                                                           no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(alpha),1,tolerance=0.02)){
                
                RMQS.data.clay[Hrow,]$clay.pred.var <- sum((alpha[1]*alpha[1]*GSM.1),
                                                           (alpha[2]*alpha[2]*GSM.2),
                                                           (alpha[3]*alpha[3]*GSM.3),
                                                           (alpha[4]*alpha[4]*GSM.4),
                                                           (alpha[5]*alpha[5]*GSM.5),
                                                           (alpha[6]*alpha[6]*GSM.6),
                                                           
                                                           (2*cor.clay1.2*sqrt(GSM.1)*sqrt(GSM.2)*alpha[1]*alpha[2]),
                                                           (2*cor.clay1.3*sqrt(GSM.1)*sqrt(GSM.3)*alpha[1]*alpha[3]),
                                                           (2*cor.clay1.4*sqrt(GSM.1)*sqrt(GSM.4)*alpha[1]*alpha[4]),
                                                           (2*cor.clay1.5*sqrt(GSM.1)*sqrt(GSM.5)*alpha[1]*alpha[5]),
                                                           (2*cor.clay1.6*sqrt(GSM.1)*sqrt(GSM.6)*alpha[1]*alpha[6]),
                                                           
                                                           (2*cor.clay2.3*sqrt(GSM.2)*sqrt(GSM.3)*alpha[2]*alpha[3]),
                                                           (2*cor.clay2.4*sqrt(GSM.2)*sqrt(GSM.4)*alpha[2]*alpha[4]),
                                                           (2*cor.clay2.5*sqrt(GSM.2)*sqrt(GSM.5)*alpha[2]*alpha[5]),
                                                           (2*cor.clay2.6*sqrt(GSM.2)*sqrt(GSM.6)*alpha[2]*alpha[6]),
                                                           
                                                           (2*cor.clay3.4*sqrt(GSM.3)*sqrt(GSM.4)*alpha[3]*alpha[4]),
                                                           (2*cor.clay3.5*sqrt(GSM.3)*sqrt(GSM.5)*alpha[3]*alpha[5]),
                                                           (2*cor.clay3.6*sqrt(GSM.3)*sqrt(GSM.6)*alpha[3]*alpha[6]),
                                                           
                                                           (2*cor.clay4.5*sqrt(GSM.4)*sqrt(GSM.5)*alpha[4]*alpha[5]),
                                                           (2*cor.clay4.6*sqrt(GSM.4)*sqrt(GSM.6)*alpha[4]*alpha[6]),
                                                           
                                                           (2*cor.clay5.6*sqrt(GSM.5)*sqrt(GSM.6)*alpha[5]*alpha[6]), na.rm=TRUE)
                
            } else {
                
                print(paste (Hrow, "error with weights for", RMQS.data.clay[Hrow,]$id_tot1 ))
            }
        }
}


rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t, alpha)


# ### SILT VARIANCE (%^2) -------------------------------------------------

### Define the arguments for this loop

### Working with var.silt
GSM <- c("var_silt_0_5","var_silt_5_15", "var_silt_15_30", "var_silt_30_60","var_silt_60_100","var_silt_100_200")

### Create empty variable
RMQS.data.silt$silt.pred.var <- NA

for(Hrow in 1:nrow(RMQS.data.silt)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.silt[Hrow,]$profondeur_hz_inf <- ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_inf > 200, 
                                                      yes= 200,
                                                      no= RMQS.data.silt[Hrow,]$profondeur_hz_inf )
        
        t <- RMQS.data.silt[Hrow,]$profondeur_hz_inf - RMQS.data.silt[Hrow,]$profondeur_hz_sup 
        
        GSM.1 <- RMQS.data.silt[Hrow, GSM[[1]]] 
        GSM.2 <- RMQS.data.silt[Hrow, GSM[[2]]]
        GSM.3 <- RMQS.data.silt[Hrow, GSM[[3]]]
        GSM.4 <- RMQS.data.silt[Hrow, GSM[[4]]]
        GSM.5 <- RMQS.data.silt[Hrow, GSM[[5]]]
        GSM.6 <- RMQS.data.silt[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[1] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            RMQS.data.silt[Hrow,]$silt.pred.var <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[2] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[2]){
            
            RMQS.data.silt[Hrow,]$silt.pred.var <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[3] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[3]){
            
            RMQS.data.silt[Hrow,]$silt.pred.var <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[4] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[4]){
            
            RMQS.data.silt[Hrow,]$silt.pred.var <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[5] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[5]){
            
            RMQS.data.silt[Hrow,]$silt.pred.var <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[6] & RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[6]){
            
            RMQS.data.silt[Hrow,]$silt.pred.var <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            alpha <- c(rep(0,6))
            
            for(i in 1:6){
                
                alpha[i] <- ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                       RMQS.data.silt[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                       RMQS.data.silt[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] &
                                       RMQS.data.silt[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i],
                                   yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                   no= ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                                  RMQS.data.silt[Hrow,]$profondeur_hz_sup <= lim.GSM.inf[i] &
                                                  RMQS.data.silt[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] & 
                                                  RMQS.data.silt[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[i],
                                              yes= (RMQS.data.silt[Hrow,]$profondeur_hz_inf - lim.GSM.sup[i])/t,
                                              no= ifelse(test=RMQS.data.silt[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[i] &
                                                             RMQS.data.silt[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                                             RMQS.data.silt[Hrow,]$profondeur_hz_inf >= lim.GSM.sup[i] & 
                                                             RMQS.data.silt[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i], 
                                                         yes= (lim.GSM.inf[i] - RMQS.data.silt[Hrow,]$profondeur_hz_sup)/t,
                                                         no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(alpha),1,tolerance=0.02)){
                
                RMQS.data.silt[Hrow,]$silt.pred.var <- sum((alpha[1]*alpha[1]*GSM.1),
                                                           (alpha[2]*alpha[2]*GSM.2),
                                                           (alpha[3]*alpha[3]*GSM.3),
                                                           (alpha[4]*alpha[4]*GSM.4),
                                                           (alpha[5]*alpha[5]*GSM.5),
                                                           (alpha[6]*alpha[6]*GSM.6),
                                                           (2*cor.silt1.2*sqrt(GSM.1)*sqrt(GSM.2)*alpha[1]*alpha[2]),
                                                           (2*cor.silt1.3*sqrt(GSM.1)*sqrt(GSM.3)*alpha[1]*alpha[3]),
                                                           (2*cor.silt2.3*sqrt(GSM.2)*sqrt(GSM.3)*alpha[2]*alpha[3]),
                                                           (2*cor.silt1.4*sqrt(GSM.1)*sqrt(GSM.4)*alpha[1]*alpha[4]),
                                                           (2*cor.silt2.4*sqrt(GSM.2)*sqrt(GSM.4)*alpha[2]*alpha[4]),
                                                           (2*cor.silt3.4*sqrt(GSM.3)*sqrt(GSM.4)*alpha[3]*alpha[4]),
                                                           (2*cor.silt1.5*sqrt(GSM.1)*sqrt(GSM.5)*alpha[1]*alpha[5]),
                                                           (2*cor.silt2.5*sqrt(GSM.2)*sqrt(GSM.5)*alpha[2]*alpha[5]),
                                                           (2*cor.silt3.5*sqrt(GSM.3)*sqrt(GSM.5)*alpha[3]*alpha[5]),
                                                           (2*cor.silt4.5*sqrt(GSM.4)*sqrt(GSM.5)*alpha[4]*alpha[5]),
                                                           (2*cor.silt1.6*sqrt(GSM.1)*sqrt(GSM.6)*alpha[1]*alpha[6]),
                                                           (2*cor.silt2.6*sqrt(GSM.2)*sqrt(GSM.6)*alpha[2]*alpha[6]),
                                                           (2*cor.silt3.6*sqrt(GSM.3)*sqrt(GSM.6)*alpha[3]*alpha[6]),
                                                           (2*cor.silt4.6*sqrt(GSM.4)*sqrt(GSM.6)*alpha[4]*alpha[6]),
                                                           (2*cor.silt5.6*sqrt(GSM.5)*sqrt(GSM.6)*alpha[5]*alpha[6]), na.rm=TRUE)
                
            } else {
                
                print(paste (Hrow, "error with weights for", RMQS.data.silt[Hrow,]$id_tot1 ))
            }
        }
}

rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t,var.pred, alpha)


# ### SAND VARIANCE (%^2) -------------------------------------------------

### Define the arguments for this loop

### Working with var.sand
GSM <- c("var_sand_0_5","var_sand_5_15", "var_sand_15_30", "var_sand_30_60","var_sand_60_100","var_sand_100_200")

### Create empty variable
RMQS.data.sand$sand.pred.var <- NA

for(Hrow in 1:nrow(RMQS.data.sand)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.sand[Hrow,]$profondeur_hz_inf <- ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_inf > 200, 
                                                      yes= 200,
                                                      no= RMQS.data.sand[Hrow,]$profondeur_hz_inf )
        
        t <- RMQS.data.sand[Hrow,]$profondeur_hz_inf - RMQS.data.sand[Hrow,]$profondeur_hz_sup 
        
        GSM.1 <- RMQS.data.sand[Hrow, GSM[[1]]] 
        GSM.2 <- RMQS.data.sand[Hrow, GSM[[2]]]
        GSM.3 <- RMQS.data.sand[Hrow, GSM[[3]]]
        GSM.4 <- RMQS.data.sand[Hrow, GSM[[4]]]
        GSM.5 <- RMQS.data.sand[Hrow, GSM[[5]]]
        GSM.6 <- RMQS.data.sand[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[1] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            RMQS.data.sand[Hrow,]$sand.pred.var <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[2] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[2]){
            
            RMQS.data.sand[Hrow,]$sand.pred.var <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[3] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[3]){
            
            RMQS.data.sand[Hrow,]$sand.pred.var <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[4] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[4]){
            
            RMQS.data.sand[Hrow,]$sand.pred.var <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[5] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[5]){
            
            RMQS.data.sand[Hrow,]$sand.pred.var <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[6] & RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[6]){
            
            RMQS.data.sand[Hrow,]$sand.pred.var <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            alpha <- c(rep(0,6))
            
            for(i in 1:6){
                
                alpha[i] <- ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                       RMQS.data.sand[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                       RMQS.data.sand[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] &
                                       RMQS.data.sand[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i],
                                   yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                   no= ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_sup <= lim.GSM.sup[i] &
                                                  RMQS.data.sand[Hrow,]$profondeur_hz_sup <= lim.GSM.inf[i] &
                                                  RMQS.data.sand[Hrow,]$profondeur_hz_inf > lim.GSM.sup[i] & 
                                                  RMQS.data.sand[Hrow,]$profondeur_hz_inf <= lim.GSM.inf[i],
                                              yes= (RMQS.data.sand[Hrow,]$profondeur_hz_inf - lim.GSM.sup[i])/t,
                                              no= ifelse(test=RMQS.data.sand[Hrow,]$profondeur_hz_sup >= lim.GSM.sup[i] &
                                                             RMQS.data.sand[Hrow,]$profondeur_hz_sup < lim.GSM.inf[i] &
                                                             RMQS.data.sand[Hrow,]$profondeur_hz_inf >= lim.GSM.sup[i] & 
                                                             RMQS.data.sand[Hrow,]$profondeur_hz_inf >= lim.GSM.inf[i], 
                                                         yes= (lim.GSM.inf[i] - RMQS.data.sand[Hrow,]$profondeur_hz_sup)/t,
                                                         no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(alpha),1,tolerance=0.02)){
                
                RMQS.data.sand[Hrow,]$sand.pred.var <- sum((alpha[1]*alpha[1]*GSM.1),
                                                           (alpha[2]*alpha[2]*GSM.2),
                                                           (alpha[3]*alpha[3]*GSM.3),
                                                           (alpha[4]*alpha[4]*GSM.4),
                                                           (alpha[5]*alpha[5]*GSM.5),
                                                           (alpha[6]*alpha[6]*GSM.6),
                                                           (2*cor.sand1.2*sqrt(GSM.1)*sqrt(GSM.2)*alpha[1]*alpha[2]),
                                                           (2*cor.sand1.3*sqrt(GSM.1)*sqrt(GSM.3)*alpha[1]*alpha[3]),
                                                           (2*cor.sand2.3*sqrt(GSM.2)*sqrt(GSM.3)*alpha[2]*alpha[3]),
                                                           (2*cor.sand1.4*sqrt(GSM.1)*sqrt(GSM.4)*alpha[1]*alpha[4]),
                                                           (2*cor.sand2.4*sqrt(GSM.2)*sqrt(GSM.4)*alpha[2]*alpha[4]),
                                                           (2*cor.sand3.4*sqrt(GSM.3)*sqrt(GSM.4)*alpha[3]*alpha[4]),
                                                           (2*cor.sand1.5*sqrt(GSM.1)*sqrt(GSM.5)*alpha[1]*alpha[5]),
                                                           (2*cor.sand2.5*sqrt(GSM.2)*sqrt(GSM.5)*alpha[2]*alpha[5]),
                                                           (2*cor.sand3.5*sqrt(GSM.3)*sqrt(GSM.5)*alpha[3]*alpha[5]),
                                                           (2*cor.sand4.5*sqrt(GSM.4)*sqrt(GSM.5)*alpha[4]*alpha[5]),
                                                           (2*cor.sand1.6*sqrt(GSM.1)*sqrt(GSM.6)*alpha[1]*alpha[6]),
                                                           (2*cor.sand2.6*sqrt(GSM.2)*sqrt(GSM.6)*alpha[2]*alpha[6]),
                                                           (2*cor.sand3.6*sqrt(GSM.3)*sqrt(GSM.6)*alpha[3]*alpha[6]),
                                                           (2*cor.sand4.6*sqrt(GSM.4)*sqrt(GSM.6)*alpha[4]*alpha[6]),
                                                           (2*cor.sand5.6*sqrt(GSM.5)*sqrt(GSM.6)*alpha[5]*alpha[6]), na.rm=TRUE)
                
            } else {
                
                print(paste (Hrow, "error with weights for", RMQS.data.sand[Hrow,]$id_tot1 ))
            }
        }
}



rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t,var.pred, alpha)


# ### COARSE ELEMENTS VARIANCE (%^2) -------------------------------------------------

### Define the arguments for this loop

### Working with var.coarse
GSM <- c("coarse_0_5.preds.sd","coarse_5_15.preds.sd","coarse_15_30.preds.sd",
         "coarse_30_60.preds.sd","coarse_60_100.preds.sd","coarse_100_200.preds.sd")

### Create empty variable
RMQS.data.coarse$abondance_eg.pred.var <- NA

for(Hrow in 1:nrow(RMQS.data.coarse)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    RMQS.data.coarse[Hrow,]$prof_inf_moy <- ifelse(test=RMQS.data.coarse[Hrow,]$prof_inf_moy > 200, 
                                                   yes= 200,
                                                   no= RMQS.data.coarse[Hrow,]$prof_inf_moy )
    
    t <- RMQS.data.coarse[Hrow,]$prof_inf_moy - RMQS.data.coarse[Hrow,]$prof_sup_moy 
        
        GSM.1 <- RMQS.data.coarse[Hrow, GSM[[1]]] 
        GSM.2 <- RMQS.data.coarse[Hrow, GSM[[2]]]
        GSM.3 <- RMQS.data.coarse[Hrow, GSM[[3]]]
        GSM.4 <- RMQS.data.coarse[Hrow, GSM[[4]]]
        GSM.5 <- RMQS.data.coarse[Hrow, GSM[[5]]]
        GSM.6 <- RMQS.data.coarse[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[1] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            RMQS.data.coarse[Hrow,]$abondance_eg.pred.var <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[2] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[2]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred.var <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[3] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[3]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred.var <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[4] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[4]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred.var <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[5] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[5]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred.var <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[6] & RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[6]){
            
            RMQS.data.coarse[Hrow,]$abondance_eg.pred.var <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            alpha <- c(rep(0,6))
            
            for(i in 1:6){
                
                alpha[i] <- ifelse(test=RMQS.data.coarse[Hrow,]$prof_sup_moy <= lim.GSM.sup[i] &
                                       RMQS.data.coarse[Hrow,]$prof_sup_moy < lim.GSM.inf[i] &
                                       RMQS.data.coarse[Hrow,]$prof_inf_moy > lim.GSM.sup[i] &
                                       RMQS.data.coarse[Hrow,]$prof_inf_moy >= lim.GSM.inf[i],
                                   yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                   no= ifelse(test=RMQS.data.coarse[Hrow,]$prof_sup_moy <= lim.GSM.sup[i] &
                                                  RMQS.data.coarse[Hrow,]$prof_sup_moy <= lim.GSM.inf[i] &
                                                  RMQS.data.coarse[Hrow,]$prof_inf_moy > lim.GSM.sup[i] & 
                                                  RMQS.data.coarse[Hrow,]$prof_inf_moy <= lim.GSM.inf[i],
                                              yes= (RMQS.data.coarse[Hrow,]$prof_inf_moy - lim.GSM.sup[i])/t,
                                              no= ifelse(test=RMQS.data.coarse[Hrow,]$prof_sup_moy >= lim.GSM.sup[i] &
                                                             RMQS.data.coarse[Hrow,]$prof_sup_moy < lim.GSM.inf[i] &
                                                             RMQS.data.coarse[Hrow,]$prof_inf_moy >= lim.GSM.sup[i] & 
                                                             RMQS.data.coarse[Hrow,]$prof_inf_moy >= lim.GSM.inf[i], 
                                                         yes= (lim.GSM.inf[i] - RMQS.data.coarse[Hrow,]$prof_sup_moy)/t,
                                                         no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(alpha),1,tolerance=0.02)){
                
                RMQS.data.coarse[Hrow,]$abondance_eg.pred.var <- sum((alpha[1]*alpha[1]*GSM.1*GSM.1),
                                                           (alpha[2]*alpha[2]*GSM.2*GSM.2),
                                                           (alpha[3]*alpha[3]*GSM.3*GSM.3),
                                                           (alpha[4]*alpha[4]*GSM.4*GSM.4),
                                                           (alpha[5]*alpha[5]*GSM.5*GSM.5),
                                                           (alpha[6]*alpha[6]*GSM.6*GSM.6),
                                                           (2*cor.coarse1.2*GSM.1*GSM.2*alpha[1]*alpha[2]),
                                                           (2*cor.coarse1.3*GSM.1*GSM.3*alpha[1]*alpha[3]),
                                                           (2*cor.coarse2.3*GSM.2*GSM.3*alpha[2]*alpha[3]),
                                                           (2*cor.coarse1.4*GSM.1*GSM.4*alpha[1]*alpha[4]),
                                                           (2*cor.coarse2.4*GSM.2*GSM.4*alpha[2]*alpha[4]),
                                                           (2*cor.coarse3.4*GSM.3*GSM.4*alpha[3]*alpha[4]),
                                                           (2*cor.coarse1.5*GSM.1*GSM.5*alpha[1]*alpha[5]),
                                                           (2*cor.coarse2.5*GSM.2*GSM.5*alpha[2]*alpha[5]),
                                                           (2*cor.coarse3.5*GSM.3*GSM.5*alpha[3]*alpha[5]),
                                                           (2*cor.coarse4.5*GSM.4*GSM.5*alpha[4]*alpha[5]),
                                                           (2*cor.coarse1.6*GSM.1*GSM.6*alpha[1]*alpha[6]),
                                                           (2*cor.coarse2.6*GSM.2*GSM.6*alpha[2]*alpha[6]),
                                                           (2*cor.coarse3.6*GSM.3*GSM.6*alpha[3]*alpha[6]),
                                                           (2*cor.coarse4.6*GSM.4*GSM.6*alpha[4]*alpha[6]),
                                                           (2*cor.coarse5.6*GSM.5*GSM.6*alpha[5]*alpha[6]), na.rm=TRUE)
                
            } else {
                
                print(paste (Hrow, "error with weights for", RMQS.data.coarse[Hrow,]$id_tot1 ))
            }
        }
}


rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t,var.pred, alpha)



# ### Calculate upper and lower prediction interval, and PICP -----------------------

RMQS.data.clay.full <- RMQS.data.clay[!is.na(RMQS.data.clay$clay.pred), ]
RMQS.data.silt.full <- RMQS.data.silt[!is.na(RMQS.data.silt$silt.pred), ] 
RMQS.data.sand.full <- RMQS.data.sand[!is.na(RMQS.data.sand$sand.pred), ]
RMQS.data.coarse.full <- RMQS.data.coarse[!is.na(RMQS.data.coarse$abondance_eg.pred), ]

### calculate upper and lower prediction interval limits
RMQS.data.clay.full$clay.pred_UPL <- (RMQS.data.clay.full$clay.pred*10) + (1.64 * sqrt(RMQS.data.clay.full$clay.pred.var)*10)
RMQS.data.clay.full$clay.pred_LPL <- (RMQS.data.clay.full$clay.pred*10) - (1.64 * sqrt(RMQS.data.clay.full$clay.pred.var)*10)

RMQS.data.silt.full$silt.pred_UPL <- (RMQS.data.silt.full$silt.pred*10) + (1.64 * sqrt(RMQS.data.silt.full$silt.pred.var)*10)
RMQS.data.silt.full$silt.pred_LPL <- (RMQS.data.silt.full$silt.pred*10) - (1.64 * sqrt(RMQS.data.silt.full$silt.pred.var)*10)

RMQS.data.sand.full$sand.pred_UPL <- (RMQS.data.sand.full$sand.pred*10) + (1.64 * sqrt(RMQS.data.sand.full$sand.pred.var)*10)
RMQS.data.sand.full$sand.pred_LPL <- (RMQS.data.sand.full$sand.pred*10) - (1.64 * sqrt(RMQS.data.sand.full$sand.pred.var)*10)

RMQS.data.coarse.full$abondance_eg.pred_UPL <- RMQS.data.coarse.full$abondance_eg.pred + (1.64 * sqrt(RMQS.data.coarse.full$abondance_eg.pred.var))
RMQS.data.coarse.full$abondance_eg.pred_LPL <- RMQS.data.coarse.full$abondance_eg.pred - (1.64 * sqrt(RMQS.data.coarse.full$abondance_eg.pred.var))


### calculate PICP
PICP <- function(table,var,UPL,LPL) {
    bMat <- c(rep(NA, nrow(table)))
    bMat <- as.numeric(table[,var]<= table[,UPL] & table[,var]>= table[,LPL])
    picp <- sum(bMat)/length(bMat) 
}

picp_clay <- PICP(table = RMQS.data.clay.full,var = "clay", UPL = "clay.pred_UPL", LPL = "clay.pred_LPL" )*100
picp_silt <- PICP(table = RMQS.data.silt.full,var = "silt", UPL = "silt.pred_UPL", LPL = "silt.pred_LPL" )*100
picp_sand <- PICP(table = RMQS.data.sand.full,var = "sand", UPL = "sand.pred_UPL", LPL = "sand.pred_LPL" )*100
picp_coarse <- PICP(table = RMQS.data.coarse.full, var = "abondance_eg", UPL = "abondance_eg.pred_UPL", LPL = "abondance_eg.pred_LPL" )*100

save.image("9-Validation.RMQS.RData")

#############################################################################################################################

# ### figure for article. Figure 3 ----------------------------------------

### Clay
RMQS.data.clay.full$mean.prof <- (RMQS.data.clay.full$profondeur_hz_sup + RMQS.data.clay.full$profondeur_hz_inf)/2
RMQS.data.clay.full$clay.pred.hor <- RMQS.data.clay.full$clay.pred*10

jpeg("Figure3a.jpeg", width=1772 , height= 1570 , units="px")
par(oma = c(4,4,4,15), mar = c(8,8,8,15))
clay.plot <- ggplot(RMQS.data.clay.full, aes(clay, clay.pred.hor, color=mean.prof))
p1 <- clay.plot + geom_point(alpha = 7/10, size = 6) + theme_bw() +
    labs(y = expression(Predicted~clay~(g~kg^-1)), 
         x = expression(Observed~clay~(g~kg^-1)))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
  xlim(0, 1000) + ylim(0,1000)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))

dev.off()


### Silt
RMQS.data.silt.full$mean.prof <- (RMQS.data.silt.full$profondeur_hz_sup + RMQS.data.silt.full$profondeur_hz_inf)/2
RMQS.data.silt.full$silt.pred.hor <- RMQS.data.silt.full$silt.pred*10

jpeg("Figure3b.jpeg", width=1772 , height= 1570 , units="px")
par(oma = c(4,4,4,15), mar = c(8,8,8,15))
silt.plot <- ggplot(RMQS.data.silt.full, aes(silt, silt.pred.hor, color=mean.prof))
p2 <-silt.plot + geom_point(alpha = 7/10, size = 6) + theme_bw() +
    labs(y = expression(Predicted~silt~(g~kg^-1)), 
         x = expression(Observed~silt~(g~kg^-1)))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
  xlim(0, 1000) + ylim(0,1000)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p2
dev.off()


### sand
RMQS.data.sand.full$mean.prof <- (RMQS.data.sand.full$profondeur_hz_sup + RMQS.data.sand.full$profondeur_hz_inf)/2
RMQS.data.sand.full$sand.pred.hor <- RMQS.data.sand.full$sand.pred*10

jpeg("Figure3c.jpeg", width=1772 , height= 1570 , units="px")
par(oma = c(4,4,4,15), mar = c(8,8,8,15))
sand.plot <- ggplot(RMQS.data.sand.full, aes(sand, sand.pred.hor, color=mean.prof))
p3<-sand.plot + geom_point(alpha = 7/10, size = 6) + theme_bw() +
    labs(y = expression(Predicted~sand~(g~kg^-1)), 
         x = expression(Observed~sand~(g~kg^-1)))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
  xlim(0, 1000) + ylim(0,1000)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p3
dev.off()

### Coarse elements

RMQS.data.coarse.full$mean.prof <- (RMQS.data.coarse.full$prof_sup_moy + RMQS.data.coarse.full$prof_inf_moy)/2

jpeg("Figure3d.jpeg", width=1772 , height= 1570 , units="px")
par(oma = c(4,4,4,15), mar = c(8,8,8,15))
coarse.plot <- ggplot(RMQS.data.coarse.full, aes(abondance_eg,abondance_eg.pred, color=mean.prof))
p4 <-coarse.plot + geom_point(alpha = 7/10, size = 6) + theme_bw() +
    labs(y = "Predicted coarse elements (%)", x ="Observed coarse elements (%)")+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
  xlim(0, 100) + ylim(0,100)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p4
dev.off()

library(grid)
# Create a text
grob1 <- grobTree(textGrob("a)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob2 <- grobTree(textGrob("b)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob3 <- grobTree(textGrob("c)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob4 <- grobTree(textGrob("d)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))

p1 <- p1 + annotation_custom(grob1)
p2 <- p2 + annotation_custom(grob2)
p3 <- p3 + annotation_custom(grob3)
p4 <- p4 + annotation_custom(grob4)

library(gridExtra)

# Plot
jpeg("Figure5.jpeg", width = 3740 , height = 3600,units="px" )
grid.arrange(p1+ annotation_custom(grob1),p2+ annotation_custom(grob2),
             p3+ annotation_custom(grob3),p4+ annotation_custom(grob4),nrow=2)
dev.off()



p4+ xlim(0,20)




###############################################################################################################################

### Plot prediction error in funciton of average horizon depth

RMQS.data.clay.full$clay.pred.error <- RMQS.data.clay.full$clay.pred.hor - RMQS.data.clay.full$clay

e1 <- ggplot(RMQS.data.clay.full, aes( clay, clay.pred.error, color=mean.prof))
e1 + geom_point(alpha = 0.2, size = 6) + theme_bw() +
    labs(y = expression(Prediction~error~clay~(g~kg^-1)), 
         x = expression(Observed~clay~(g~kg^-1)))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_abline(slope=-0.629451, intercept = 151.170282, size=1, colour="black")


e1 <- ggplot(RMQS.data.clay.full, aes( mean.prof, clay.pred.error, color=clay))
e1 + geom_point(alpha = 0.2, size = 6) + theme_bw() +
    labs(y = expression(Prediction~error~clay~(g~kg^-1)), 
         x = "Mean horizon depth (cm)" )+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, 
                                 name=expression(Observed~clay~(g~kg^-1))) +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
# +    geom_abline(slope=0.22112, intercept = -25.35823, size=1, colour="black")


### Silt

RMQS.data.silt.full$silt.pred.error <- RMQS.data.silt.full$silt.pred.hor - RMQS.data.silt.full$silt

e1 <- ggplot(RMQS.data.silt.full, aes( silt, silt.pred.error, color=mean.prof))
e1 + geom_point(alpha = 0.2, size = 6) + theme_bw() +
    labs(y = expression(Prediction~error~silt~(g~kg^-1)), 
         x = expression(Observed~silt~(g~kg^-1)))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) +
    geom_abline(slope= -0.4776, intercept = 211.1320, size=1, colour="black")

summary(lm(RMQS.data.silt.full$silt.pred.error ~ RMQS.data.silt.full$silt))

e1 <- ggplot(RMQS.data.silt.full, aes( mean.prof, silt.pred.error, color=silt))
e1 + geom_point(alpha = 0.2, size = 6) + theme_bw() +
    labs(y = expression(Prediction~error~silt~(g~kg^-1)), 
         x = "Mean horizon depth (cm)" )+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, 
                                 name=expression(Observed~silt~(g~kg^-1))) +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_abline(slope=-0.08603, intercept = 23.15777, size=1, colour="black")

summary(lm(RMQS.data.silt.full$silt.pred.error ~ RMQS.data.silt.full$mean.prof))



### sand

RMQS.data.sand.full$sand.pred.error <- RMQS.data.sand.full$sand.pred.hor - RMQS.data.sand.full$sand

e1 <- ggplot(RMQS.data.sand.full, aes( sand, sand.pred.error, color=mean.prof))
e1 + geom_point(alpha = 0.2, size = 6) + theme_bw() +
    labs(y = expression(Prediction~error~sand~(g~kg^-1)), 
         x = expression(Observed~sand~(g~kg^-1)))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) +
    geom_abline(slope= -0.4776, intercept = 151.533095, size=1, colour="black")

summary(lm(RMQS.data.sand.full$sand.pred.error ~ RMQS.data.sand.full$sand))

e1 <- ggplot(RMQS.data.sand.full, aes( mean.prof, sand.pred.error, color=sand))
e1 + geom_point(alpha = 0.2, size = 6) + theme_bw() +
    labs(y = expression(Prediction~error~sand~(g~kg^-1)), 
         x = "Mean horizon depth (cm)" )+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, 
                                 name=expression(Observed~sand~(g~kg^-1))) +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_abline(slope=-0.13838, intercept = 3.57319, size=1, colour="black")

summary(lm(RMQS.data.sand.full$sand.pred.error ~ RMQS.data.sand.full$mean.prof))


### coarse

RMQS.data.coarse.full$coarse.pred.error <- RMQS.data.coarse.full$abondance_eg.pred - RMQS.data.coarse.full$abondance_eg

e1 <- ggplot(RMQS.data.coarse.full[RMQS.data.coarse.full$abondance_eg <20,], aes( abondance_eg, coarse.pred.error, color=mean.prof))
E1 <- e1 + geom_point(alpha = 0.5, size = 8) + theme_bw() +
    labs(y = "Prediction error coarse (%)", 
         x = "Observed coarse (%)")+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm")) +
    geom_abline(slope= -0.90824, intercept = 16.00567 , size=1, colour="black")

summary(lm(RMQS.data.coarse.full[RMQS.data.coarse.full$abondance_eg <20,]$coarse.pred.error ~ RMQS.data.coarse.full[RMQS.data.coarse.full$abondance_eg <20,]$abondance_eg))


plot(lm(RMQS.data.coarse.full$coarse.pred.error ~ RMQS.data.coarse.full$abondance_eg+RMQS.data.coarse.full$mean.prof))

e1 <- ggplot(RMQS.data.coarse.full, aes(prof_sup_moy, coarse.pred.error, color=abondance_eg))
E2 <- e1 + geom_point(alpha = 0.8, size = 8) + theme_bw() +
    labs(y = "Prediction error coarse (%)", 
         x = "Mean horizon depth (cm)" ) +
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, 
                                 name=expression(Observed~coarse~(g~kg^-1))) +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_abline(slope=0.0631, intercept = 1.5554 , size=1, colour="black")

summary(lm(RMQS.data.coarse.full$coarse.pred.error ~ RMQS.data.coarse.full$prof_sup_moy))


# Plot
jpeg("Figure4b.jpeg", width = 3740 , height = 1570,units="px" )
grid.arrange(E1+ annotation_custom(grob1),E2+ annotation_custom(grob2),ncol=2)
dev.off()

jpeg("Figure4c.jpeg", width = 1570 , height = 1570,units="px" )
E2
dev.off()


#### End of the script