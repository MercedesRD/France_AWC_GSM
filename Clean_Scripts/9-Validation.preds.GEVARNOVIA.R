##########################################################################################################################################

#### independent validation GEVARNOVIA soil moisture content at pF = 2.0, and pF = 4.2,
#### from SPATIAL PREDICYIONS OF SOIL MOISTURE CONTENTS

############# RESULTS of Cubist + residual cokriging (ISATIS software)

### Date 12/06/2018
### Author: Mercedes Roman Dobarco

####### Packages
library(soiltexture)
library(plyr)
library(knitr)
library(Hmisc)
library(MASS)
library(randomForest)
library(gbm)
library(Cubist)
library(caret)
library(wesanderson)
library(ggplot2)
library(randomForest)
library(quantregForest)
library(boot)
library(car)
library(sp)
library(gstat)
library(raster)
library(rgeos)
library(rgdal)
library(fBasics)
library(compositions)
library(rgr)
library(gdata)
library(ithir)
library(gridExtra)

# ### 1. Subset the data used for independent validation ------------------

### Load the data
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
lyatou_bdd <- read.csv(paste0(HomeDir,"Input/BDD_LYATOU_FPT_MVA.toute_methode.csv"), header = T, sep = ";")
lyatou_bdd <- lyatou_bdd[,1:212]

#### Eliminate INRA observations
lyatou_bdd <- lyatou_bdd[lyatou_bdd$Source != "INRA",]
names(lyatou_bdd)

#### when TEXT_controle == difference D
summary(lyatou_bdd$TEXT_controle)
summary(lyatou_bdd$CarbOrg)
summary(lyatou_bdd$DaTerreFine_mesuré)

#### subset and change name variables as they are for the PTF
bdd1 <- lyatou_bdd[lyatou_bdd$TEXT_controle == "difference D",]
bdd1$clay <- bdd1$TEXT_D_Argile * 100
bdd1$sand <- bdd1$TEXT_D_Sable * 100
bdd1$silt <- bdd1$TEXT_D_Limon * 100
bdd1$carb <- bdd1$CarbOrg
bdd1$Da_cyl <- bdd1$DaTerreFine_mesuré
bdd1$SMpF2.0 <- (bdd1$HpF2.0_labo_mesure/100) * bdd1$DaTerreFine_mesuré
bdd1$SMpF4.2 <- (bdd1$HpF4.2_labo_mesure/100) * bdd1$DaTerreFine_mesuré

#### when TEXT_controle == difference ND
#### subset and change name variables as they are for the PTF
bdd2 <- lyatou_bdd[lyatou_bdd$TEXT_controle == "difference ND",]
bdd2$clay <- bdd2$TEXT_ND_Argile * 100
bdd2$sand <- bdd2$TEXT_ND_Sable * 100
bdd2$silt <- bdd2$TEXT_ND_Limon * 100
bdd2$carb <- bdd2$CarbOrg
bdd2$Da_cyl <- bdd2$DaTerreFine_mesuré
bdd2$SMpF2.0 <- (bdd2$HpF2.0_labo_mesure/100) * bdd2$DaTerreFine_mesuré
bdd2$SMpF4.2 <- (bdd2$HpF4.2_labo_mesure/100) * bdd2$DaTerreFine_mesuré

#### when TEXT_controle == OK (D)
#### subset and change name variables as they are for the PTF
bdd3 <- lyatou_bdd[lyatou_bdd$TEXT_controle == "OK (D)",]
bdd3$clay <- bdd3$TEXT_D_Argile * 100
bdd3$sand <- bdd3$TEXT_D_Sable * 100
bdd3$silt <- bdd3$TEXT_D_Limon * 100
bdd3$carb <- bdd3$CarbOrg
bdd3$Da_cyl <- bdd3$DaTerreFine_mesuré
bdd3$SMpF2.0 <- (bdd3$HpF2.0_labo_mesure/100) * bdd3$DaTerreFine_mesuré
bdd3$SMpF4.2 <- (bdd3$HpF4.2_labo_mesure/100) * bdd3$DaTerreFine_mesuré

#### when TEXT_controle == OK (ND)
#### subset and change name variables as they are for the PTF
bdd4 <- lyatou_bdd[lyatou_bdd$TEXT_controle == "OK (ND)",]
bdd4$clay <- bdd4$TEXT_ND_Argile * 100
bdd4$sand <- bdd4$TEXT_ND_Sable * 100
bdd4$silt <- bdd4$TEXT_ND_Limon * 100
bdd4$carb <- bdd4$CarbOrg
bdd4$Da_cyl <- bdd4$DaTerreFine_mesuré
bdd4$SMpF2.0 <- (bdd4$HpF2.0_labo_mesure/100) * bdd4$DaTerreFine_mesuré
bdd4$SMpF4.2 <- (bdd4$HpF4.2_labo_mesure/100) * bdd4$DaTerreFine_mesuré

##### Now join the 4 subsets
bdd.validation <- rbind(bdd1,bdd2,bdd3,bdd4)
rm(bdd1,bdd2,bdd3,bdd4,lyatou_bdd)

bdd.validation.w <- bdd.validation[!is.na(bdd.validation$HpF2.0_labo_mesure) | !is.na(bdd.validation$HpF4.2_labo_mesure),]

############################################################################################################################

# #### Load raster predictions --------------------------------------------


### Load the spatial predictions of soil moisture content for both potentials

### LOAD the mean values of SMFC (pf  = 2.0) for the 6 depth layers
### Elementary volumetric soil moisture content (cm3/cm3)
setwd(paste0(HomeDir,"Clean_Output/6-Calculate_AWC"))
SMFC_0_5 <- raster("SMFC_0_5.tif")
SMFC_5_15 <- raster("SMFC_5_15.tif")
SMFC_15_30 <- raster("SMFC_15_30.tif")
SMFC_30_60 <- raster("SMFC_30_60.tif")
SMFC_60_100 <- raster("SMFC_60_100.tif")
SMFC_100_200 <- raster("SMFC_100_200.tif")

### Load the mean values of SMPWP (pF = 4.2) for the 6 depth layers
### Elementary volumetric soil moisture content (cm3/cm3)
SMPWP_0_5 <- raster("SMPWP_0_5.tif")
SMPWP_5_15 <- raster("SMPWP_5_15.tif")
SMPWP_15_30 <- raster("SMPWP_15_30.tif")
SMPWP_30_60 <- raster("SMPWP_30_60.tif")
SMPWP_60_100 <- raster("SMPWP_60_100.tif")
SMPWP_100_200 <- raster("SMPWP_100_200.tif")

### Load the uncertainty of the predictions
setwd(paste0(HomeDir,"Clean_Output/7.1-AWC_Taylor_0_5"))
Total_var_w20.0_5 <- raster("Total_var_w20.0_5.tif")
Total_var_w42.0_5 <- raster("Total_var_w42.0_5.tif")
setwd(paste0(HomeDir,"Clean_Output/7.2-AWC_Taylor_5_15"))
Total_var_w20.5_15 <- raster("Total_var_w20.5_15.tif")
Total_var_w42.5_15 <- raster("Total_var_w42.5_15.tif")
setwd(paste0(HomeDir,"Clean_Output/7.3-AWC_Taylor_15_30"))
Total_var_w20.15_30 <- raster("Total_var_w20.15_30.tif")
Total_var_w42.15_30 <- raster("Total_var_w42.15_30.tif")
setwd(paste0(HomeDir,"Clean_Output/7.4-AWC_Taylor_30_60"))
Total_var_w20.30_60 <- raster("Total_var_w20.30_60.tif")
Total_var_w42.30_60 <- raster("Total_var_w42.30_60.tif")
setwd(paste0(HomeDir,"Clean_Output/7.5-AWC_Taylor_60_100"))
Total_var_w20.60_100 <- raster("Total_var_w20.60_100.tif")
Total_var_w42.60_100 <- raster("Total_var_w42.60_100.tif")
setwd(paste0(HomeDir,"Clean_Output/7.6-AWC_Taylor_100_200"))
Total_var_w20.100_200 <- raster("Total_var_w20.100_200.tif")
Total_var_w42.100_200 <- raster("Total_var_w42.100_200.tif")

### change owrking directory
dir.create(paste0(HomeDir,"Clean_Output/9-Validation.preds.GEVARNOVIA"))
setwd(paste0(HomeDir,"Clean_Output/9-Validation.preds.GEVARNOVIA"))

### Create raster stack
predictions <- stack(SMFC_0_5,SMFC_5_15,SMFC_15_30,SMFC_30_60,SMFC_60_100, SMFC_100_200,
                   SMPWP_0_5,SMPWP_5_15,SMPWP_15_30,SMPWP_30_60,SMPWP_60_100, SMPWP_100_200,
                   Total_var_w20.0_5 ,Total_var_w20.5_15, Total_var_w20.15_30,
                   Total_var_w20.30_60,Total_var_w20.60_100,Total_var_w20.100_200,
                   Total_var_w42.0_5,Total_var_w42.5_15,Total_var_w42.15_30,
                   Total_var_w42.30_60,Total_var_w42.60_100,Total_var_w42.100_200)
proj4string(predictions)

#### Extract predictions at independent validation points

### transform into spatial
str(bdd.validation.w)
bdd.validation.s <- bdd.validation.w

### select only those observations where the upper limit of the horizon is above 200 cm
bdd.validation.s <- bdd.validation.s[bdd.validation.s$Horizon_cote_superieure  <= 200,]

### Transform into spatial
coordinates(bdd.validation.s) <- ~ Lambert.93.X  + Lambert.93.Y
proj4string(bdd.validation.s) <- CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

### Plot the locations of the horizons
departments<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "DEPARTEMENT")
departments <- spTransform(departments, CRS("+init=epsg:2154"))
Fr.ggmap <- fortify(departments, region="CODE_DEPT") # pour le mettre au format data.frame
INDEP.bdd <- ggplot(Fr.ggmap, aes(long, lat)) + geom_polygon(aes(group=group), fill=NA, colour="black", size=.3)+
    geom_point(data=bdd.validation.w, aes(x=Lambert.93.X, y=Lambert.93.Y),colour="mediumblue" ,size=3, pch=19) + theme_bw() +coord_equal()

hist(bdd.validation.w$Horizon_profondeur_moyenne)
hist(bdd.validation.w$Horizon_cote_superieure)
hist(bdd.validation.w$Horizon_cote_inferieure)

### Extract predictions
pred_at_arvalis <- extract(y=bdd.validation.s,x=predictions, method="simple",df=TRUE )

### Attach to dataframe
### not spatial
bdd.validation.s <- as.data.frame(bdd.validation.s)

###Cbind the predictions
bdd.validation.s<- cbind(bdd.validation.s, pred_at_arvalis)


# ### Calculate wheighed average of SMFC for GERANOVIA horizons and VALIDATE  ----------

### Define the arguments for this loop
### Working with silt
GSM <- c("SMFC_0_5","SMFC_5_15","SMFC_15_30","SMFC_30_60","SMFC_60_100", "SMFC_100_200")

### Create empty variable
bdd.validation.s$SMFC.pred <- NA

bdd.validation.s[1:10, c("Horizon_cote_superieure","Horizon_cote_inferieure","Horizon_profondeur_moyenne")]


for(Hrow in 1:nrow(bdd.validation.s)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    bdd.validation.s[Hrow,]$Horizon_cote_inferieure <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_inferieure > 200, 
                                                      yes= 200,
                                                      no= bdd.validation.s[Hrow,]$Horizon_cote_inferieure )
    
    t <- bdd.validation.s[Hrow,]$Horizon_cote_inferieure - bdd.validation.s[Hrow,]$Horizon_cote_superieure 
        
        GSM.1 <- bdd.validation.s[Hrow, GSM[[1]]] 
        GSM.2 <- bdd.validation.s[Hrow, GSM[[2]]]
        GSM.3 <- bdd.validation.s[Hrow, GSM[[3]]]
        GSM.4 <- bdd.validation.s[Hrow, GSM[[4]]]
        GSM.5 <- bdd.validation.s[Hrow, GSM[[5]]]
        GSM.6 <- bdd.validation.s[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[1] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            bdd.validation.s[Hrow,]$SMFC.pred <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[2] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[2]){
            
            bdd.validation.s[Hrow,]$SMFC.pred <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[3] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[3]){
            
            bdd.validation.s[Hrow,]$SMFC.pred <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[4] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[4]){
            
            bdd.validation.s[Hrow,]$SMFC.pred <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[5] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[5]){
            
            bdd.validation.s[Hrow,]$SMFC.pred <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[6] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[6]){
            
            bdd.validation.s[Hrow,]$SMFC.pred <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            weights <- c(rep(0,6))
            
            for(i in 1:6){
                
                weights[i] <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                         bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                         bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] &
                                         bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i],
                                     yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                     no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                                    bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.inf[i] &
                                                    bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] & 
                                                    bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[i],
                                                yes= (bdd.validation.s[Hrow,]$Horizon_cote_inferieure - lim.GSM.sup[i])/t,
                                                no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[i] &
                                                               bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                                               bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.sup[i] & 
                                                               bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i], 
                                                           yes= (lim.GSM.inf[i] - bdd.validation.s[Hrow,]$Horizon_cote_superieure)/t,
                                                           no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(weights),1,tolerance=0.02)){
                
                bdd.validation.s[Hrow,]$SMFC.pred <- sum(weights*GSM.preds)
                
            } else {
                
                print(paste (Hrow, "error with weights for", bdd.validation.s[Hrow,]$id_tot1 ))
            }
        }
}



rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t,GSM.preds, weights)

summary(bdd.validation.s$SMFC.pred) ### there is one missing value, that bellongs to the last horizon

bdd.validation.s <- bdd.validation.s[!is.na(bdd.validation.s$SMFC.pred),]

#### Validation statistics for weighed averages
par(mfrow=c(1,1))
SMFC_valid <-goof(observed = bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMpF2.0,
                  predicted = bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMFC.pred, plot.it = TRUE, type="DSM")


plot(bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMpF2.0,
     bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMFC.pred, xlim=c(0,0.6), ylim=c(0,0.5),
     pch=19,cex=0.7,col="cyan3", xlab=" ",   ylab=" ", main="SMFC");abline(0,1)

length(bdd.validation.s[!is.na(bdd.validation.s$SMFC.pred),]$SMpF2.0) ### 309

# ### Calculate wheighed average of SMPWP for GERANOVIA horizons and VALIDATE  ----------

### Define the arguments for this loop
### Working with silt
GSM <- c("SMPWP_0_5","SMPWP_5_15","SMPWP_15_30","SMPWP_30_60","SMPWP_60_100", "SMPWP_100_200")

### Create empty variable
bdd.validation.s$SMPWP.pred <- NA

bdd.validation.s[1:10, c("Horizon_cote_superieure","Horizon_cote_inferieure","Horizon_profondeur_moyenne")]


for(Hrow in 1:nrow(bdd.validation.s)){
    
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    bdd.validation.s[Hrow,]$Horizon_cote_inferieure <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_inferieure > 200, 
                                                              yes= 200,
                                                              no= bdd.validation.s[Hrow,]$Horizon_cote_inferieure )
        
        t <- bdd.validation.s[Hrow,]$Horizon_cote_inferieure - bdd.validation.s[Hrow,]$Horizon_cote_superieure 
        
        GSM.1 <- bdd.validation.s[Hrow, GSM[[1]]] 
        GSM.2 <- bdd.validation.s[Hrow, GSM[[2]]]
        GSM.3 <- bdd.validation.s[Hrow, GSM[[3]]]
        GSM.4 <- bdd.validation.s[Hrow, GSM[[4]]]
        GSM.5 <- bdd.validation.s[Hrow, GSM[[5]]]
        GSM.6 <- bdd.validation.s[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[1] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            bdd.validation.s[Hrow,]$SMPWP.pred <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[2] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[2]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[3] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[3]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[4] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[4]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[5] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[5]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[6] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[6]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            weights <- c(rep(0,6))
            
            for(i in 1:6){
                
                weights[i] <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                         bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                         bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] &
                                         bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i],
                                     yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                     no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                                    bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.inf[i] &
                                                    bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] & 
                                                    bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[i],
                                                yes= (bdd.validation.s[Hrow,]$Horizon_cote_inferieure - lim.GSM.sup[i])/t,
                                                no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[i] &
                                                               bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                                               bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.sup[i] & 
                                                               bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i], 
                                                           yes= (lim.GSM.inf[i] - bdd.validation.s[Hrow,]$Horizon_cote_superieure)/t,
                                                           no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(weights),1,tolerance=0.02)){
                
                bdd.validation.s[Hrow,]$SMPWP.pred <- sum(weights*GSM.preds)
                
            } else {
                
                print(paste (Hrow, "error with weights for", bdd.validation.s[Hrow,]$id_tot1 ))
            }
        }
    }



rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t,GSM.preds, weights)

summary(bdd.validation.s$SMPWP.pred) 

bdd.validation.s <- bdd.validation.s[!is.na(bdd.validation.s$SMPWP.pred),]

#### Validation statistics for weighed averages
par(mfrow=c(1,1))
SMPWP_valid <-goof(observed = bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMpF4.2,
                  predicted = bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMPWP.pred, plot.it = TRUE, type="DSM")


plot(bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMpF4.2,
     bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMPWP.pred, xlim=c(0,0.42), ylim=c(0,0.4),
     pch=19,cex=0.7,col="dodgerblue1", xlab=" ",   ylab=" ", main="SMPWP");abline(0,1)
length(bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMpF4.2) ## 308


#################################################################################################################################################


# ### Correlations --------------------------------------------------------

cor.SMFC1.2 <- cor.test(x=getValues(predictions[[1]]), y=getValues(predictions[[2]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC1.3 <- cor.test(x=getValues(predictions[[1]]), y=getValues(predictions[[3]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC1.4 <- cor.test(x=getValues(predictions[[1]]), y=getValues(predictions[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC1.5 <- cor.test(x=getValues(predictions[[1]]), y=getValues(predictions[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC1.6 <- cor.test(x=getValues(predictions[[1]]), y=getValues(predictions[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC2.3 <- cor.test(x=getValues(predictions[[2]]), y=getValues(predictions[[3]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC2.4 <- cor.test(x=getValues(predictions[[2]]), y=getValues(predictions[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC2.5 <- cor.test(x=getValues(predictions[[2]]), y=getValues(predictions[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC2.6 <- cor.test(x=getValues(predictions[[2]]), y=getValues(predictions[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC3.4 <- cor.test(x=getValues(predictions[[3]]), y=getValues(predictions[[4]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC3.5 <- cor.test(x=getValues(predictions[[3]]), y=getValues(predictions[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC3.6 <- cor.test(x=getValues(predictions[[3]]), y=getValues(predictions[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC4.5 <- cor.test(x=getValues(predictions[[4]]), y=getValues(predictions[[5]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC4.6 <- cor.test(x=getValues(predictions[[4]]), y=getValues(predictions[[6]]), na.rm=TRUE, method="pearson")$estimate
cor.SMFC5.6 <- cor.test(x=getValues(predictions[[5]]), y=getValues(predictions[[6]]), na.rm=TRUE, method="pearson")$estimate

cor.SMPWP1.2 <- cor.test(x=getValues(predictions[[7]]), y=getValues(predictions[[8]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP1.3 <- cor.test(x=getValues(predictions[[7]]), y=getValues(predictions[[9]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP1.4 <- cor.test(x=getValues(predictions[[7]]), y=getValues(predictions[[10]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP1.5 <- cor.test(x=getValues(predictions[[7]]), y=getValues(predictions[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP1.6 <- cor.test(x=getValues(predictions[[7]]), y=getValues(predictions[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP2.3 <- cor.test(x=getValues(predictions[[8]]), y=getValues(predictions[[9]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP2.4 <- cor.test(x=getValues(predictions[[8]]), y=getValues(predictions[[10]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP2.5 <- cor.test(x=getValues(predictions[[8]]), y=getValues(predictions[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP2.6 <- cor.test(x=getValues(predictions[[8]]), y=getValues(predictions[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP3.4 <- cor.test(x=getValues(predictions[[9]]), y=getValues(predictions[[10]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP3.5 <- cor.test(x=getValues(predictions[[9]]), y=getValues(predictions[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP3.6 <- cor.test(x=getValues(predictions[[9]]), y=getValues(predictions[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP4.5 <- cor.test(x=getValues(predictions[[10]]), y=getValues(predictions[[11]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP4.6 <- cor.test(x=getValues(predictions[[10]]), y=getValues(predictions[[12]]), na.rm=TRUE, method="pearson")$estimate
cor.SMPWP5.6 <- cor.test(x=getValues(predictions[[11]]), y=getValues(predictions[[12]]), na.rm=TRUE, method="pearson")$estimate



### Now, Calculate the variance of the average, taking into account the correlation between layers......

### Start here Tuesday morning


# ### Calculate variance of SMFC ------------------------------------------

### Define the arguments for this loop

### Working with var.SMFC
GSM <- c( "Total_var_w20.0_5" ,"Total_var_w20.5_15", "Total_var_w20.15_30",
          "Total_var_w20.30_60","Total_var_w20.60_100","Total_var_w20.100_200")

### Create empty variable
bdd.validation.s$SMFC.pred.var <- NA

for(Hrow in 1:nrow(bdd.validation.s)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    bdd.validation.s[Hrow,]$Horizon_cote_inferieure <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_inferieure > 200, 
                                                              yes= 200,
                                                              no= bdd.validation.s[Hrow,]$Horizon_cote_inferieure )
        
        t <- bdd.validation.s[Hrow,]$Horizon_cote_inferieure - bdd.validation.s[Hrow,]$Horizon_cote_superieure 
        
        GSM.1 <- bdd.validation.s[Hrow, GSM[[1]]] 
        GSM.2 <- bdd.validation.s[Hrow, GSM[[2]]]
        GSM.3 <- bdd.validation.s[Hrow, GSM[[3]]]
        GSM.4 <- bdd.validation.s[Hrow, GSM[[4]]]
        GSM.5 <- bdd.validation.s[Hrow, GSM[[5]]]
        GSM.6 <- bdd.validation.s[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[1] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            bdd.validation.s[Hrow,]$SMFC.pred.var <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[2] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[2]){
            
            bdd.validation.s[Hrow,]$SMFC.pred.var <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[3] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[3]){
            
            bdd.validation.s[Hrow,]$SMFC.pred.var <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[4] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[4]){
            
            bdd.validation.s[Hrow,]$SMFC.pred.var <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[5] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[5]){
            
            bdd.validation.s[Hrow,]$SMFC.pred.var <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[6] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[6]){
            
            bdd.validation.s[Hrow,]$SMFC.pred.var <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            alpha <- c(rep(0,6))
            
            for(i in 1:6){
                
                alpha[i] <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                       bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                       bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] &
                                       bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i],
                                   yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                   no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                                  bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.inf[i] &
                                                  bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] & 
                                                  bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[i],
                                              yes= (bdd.validation.s[Hrow,]$Horizon_cote_inferieure - lim.GSM.sup[i])/t,
                                              no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[i] &
                                                             bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                                             bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.sup[i] & 
                                                             bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i], 
                                                         yes= (lim.GSM.inf[i] - bdd.validation.s[Hrow,]$Horizon_cote_superieure)/t,
                                                         no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(alpha),1,tolerance=0.02)){
                
                bdd.validation.s[Hrow,]$SMFC.pred.var <- sum((alpha[1]*alpha[1]*GSM.1),
                                                           (alpha[2]*alpha[2]*GSM.2),
                                                           (alpha[3]*alpha[3]*GSM.3),
                                                           (alpha[4]*alpha[4]*GSM.4),
                                                           (alpha[5]*alpha[5]*GSM.5),
                                                           (alpha[6]*alpha[6]*GSM.6),
                                                           
                                                           (2*cor.SMFC1.2*sqrt(GSM.1)*sqrt(GSM.2)*alpha[1]*alpha[2]),
                                                           (2*cor.SMFC1.3*sqrt(GSM.1)*sqrt(GSM.3)*alpha[1]*alpha[3]),
                                                           (2*cor.SMFC1.4*sqrt(GSM.1)*sqrt(GSM.4)*alpha[1]*alpha[4]),
                                                           (2*cor.SMFC1.5*sqrt(GSM.1)*sqrt(GSM.5)*alpha[1]*alpha[5]),
                                                           (2*cor.SMFC1.6*sqrt(GSM.1)*sqrt(GSM.6)*alpha[1]*alpha[6]),
                                                           
                                                           (2*cor.SMFC2.3*sqrt(GSM.2)*sqrt(GSM.3)*alpha[2]*alpha[3]),                                                                                                                    (2*cor.SMFC2.4*sqrt(GSM.2)*sqrt(GSM.4)*alpha[2]*alpha[4]),
                                                           (2*cor.SMFC2.4*sqrt(GSM.2)*sqrt(GSM.4)*alpha[2]*alpha[4]),
                                                           (2*cor.SMFC2.5*sqrt(GSM.2)*sqrt(GSM.5)*alpha[2]*alpha[5]),
                                                           (2*cor.SMFC2.6*sqrt(GSM.2)*sqrt(GSM.6)*alpha[2]*alpha[6]),
                                                           
                                                           (2*cor.SMFC3.4*sqrt(GSM.3)*sqrt(GSM.4)*alpha[3]*alpha[4]),
                                                           (2*cor.SMFC3.5*sqrt(GSM.3)*sqrt(GSM.5)*alpha[3]*alpha[5]),
                                                           (2*cor.SMFC3.6*sqrt(GSM.3)*sqrt(GSM.6)*alpha[3]*alpha[6]),
                                                           
                                                           (2*cor.SMFC4.5*sqrt(GSM.4)*sqrt(GSM.5)*alpha[4]*alpha[5]),
                                                           (2*cor.SMFC4.6*sqrt(GSM.4)*sqrt(GSM.6)*alpha[4]*alpha[6]),
                                                           
                                                           (2*cor.SMFC5.6*sqrt(GSM.5)*sqrt(GSM.6)*alpha[5]*alpha[6]), na.rm=TRUE)
                
            } else {
                
                print(paste (Hrow, "error with weights for", bdd.validation.s[Hrow,]$id_tot1 ))
            }
        }
   }

rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t, alpha, GSM.preds)


# ### SMPWP VARIANCE (%^2) -------------------------------------------------

### Define the arguments for this loop

### Working with var.SMPWP
GSM <- c("Total_var_w42.0_5","Total_var_w42.5_15","Total_var_w42.15_30",
         "Total_var_w42.30_60","Total_var_w42.60_100","Total_var_w42.100_200")

### Create empty variable
bdd.validation.s$SMPWP.pred.var <- NA

for(Hrow in 1:nrow(bdd.validation.s)){
    
    print(Hrow)
    ### If the RMQS horizon is deeper than 2 m?
    ### I constrain to 200 cm depth
    bdd.validation.s[Hrow,]$Horizon_cote_inferieure <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_inferieure > 200, 
                                                              yes= 200,
                                                              no= bdd.validation.s[Hrow,]$Horizon_cote_inferieure )
        
        t <- bdd.validation.s[Hrow,]$Horizon_cote_inferieure - bdd.validation.s[Hrow,]$Horizon_cote_superieure 
        
        GSM.1 <- bdd.validation.s[Hrow, GSM[[1]]] 
        GSM.2 <- bdd.validation.s[Hrow, GSM[[2]]]
        GSM.3 <- bdd.validation.s[Hrow, GSM[[3]]]
        GSM.4 <- bdd.validation.s[Hrow, GSM[[4]]]
        GSM.5 <- bdd.validation.s[Hrow, GSM[[5]]]
        GSM.6 <- bdd.validation.s[Hrow, GSM[[6]]]
        
        lim.GSM.sup <- c(0,5,15,30,60,100)
        lim.GSM.inf <- c(5,15,30,60,100,200)
        GSM.preds <- c(GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6)
        
        #### 1. the RMQS layer falls within a GSM layer
        
        ## 1.1 In the 0-5 cm
        
        if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[1] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[1]) {  
            ## we are lucky, and our RMQS horizon falls in the first layer
            bdd.validation.s[Hrow,]$SMPWP.pred.var <- GSM.1
            
            ## 1.2 In the 5-15 cm    
        } else if(bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[2] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[2]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred.var <- GSM.2
            
            ## 1.3 In the 15-30 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[3] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[3]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred.var <- GSM.3
            
            ## 1.4 In the 30-60 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[4] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[4]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred.var <- GSM.4
            
            ## 1.5 In the 60-100 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[5] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[5]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred.var <- GSM.5
            
            ## 1.6 In the 100-200 cm
        } else if (bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[6] & bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[6]){
            
            bdd.validation.s[Hrow,]$SMPWP.pred.var <- GSM.6
            
        } else {
            
            ### 2. the RMQS layer covers several GSM horizons
            
            ### Define weights for each GSM layer
            alpha <- c(rep(0,6))
            
            for(i in 1:6){
                
                alpha[i] <- ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                       bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                       bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] &
                                       bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i],
                                   yes= (lim.GSM.inf[i]- lim.GSM.sup[i])/t, 
                                   no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.sup[i] &
                                                  bdd.validation.s[Hrow,]$Horizon_cote_superieure <= lim.GSM.inf[i] &
                                                  bdd.validation.s[Hrow,]$Horizon_cote_inferieure > lim.GSM.sup[i] & 
                                                  bdd.validation.s[Hrow,]$Horizon_cote_inferieure <= lim.GSM.inf[i],
                                              yes= (bdd.validation.s[Hrow,]$Horizon_cote_inferieure - lim.GSM.sup[i])/t,
                                              no= ifelse(test=bdd.validation.s[Hrow,]$Horizon_cote_superieure >= lim.GSM.sup[i] &
                                                             bdd.validation.s[Hrow,]$Horizon_cote_superieure < lim.GSM.inf[i] &
                                                             bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.sup[i] & 
                                                             bdd.validation.s[Hrow,]$Horizon_cote_inferieure >= lim.GSM.inf[i], 
                                                         yes= (lim.GSM.inf[i] - bdd.validation.s[Hrow,]$Horizon_cote_superieure)/t,
                                                         no= 0 )))
                
                ### check that the sum of weights is 1
                
            } 
            
            if (all.equal(sum(alpha),1,tolerance=0.02)){
                
                bdd.validation.s[Hrow,]$SMPWP.pred.var <- sum((alpha[1]*alpha[1]*GSM.1),
                                                           (alpha[2]*alpha[2]*GSM.2),
                                                           (alpha[3]*alpha[3]*GSM.3),
                                                           (alpha[4]*alpha[4]*GSM.4),
                                                           (alpha[5]*alpha[5]*GSM.5),
                                                           (alpha[6]*alpha[6]*GSM.6),
                                                           
                                                           (2*cor.SMPWP1.2*sqrt(GSM.1)*sqrt(GSM.2)*alpha[1]*alpha[2]),
                                                           (2*cor.SMPWP1.3*sqrt(GSM.1)*sqrt(GSM.3)*alpha[1]*alpha[3]),
                                                           (2*cor.SMPWP1.4*sqrt(GSM.1)*sqrt(GSM.4)*alpha[1]*alpha[4]),
                                                           (2*cor.SMPWP1.5*sqrt(GSM.1)*sqrt(GSM.5)*alpha[1]*alpha[5]),
                                                           (2*cor.SMPWP1.6*sqrt(GSM.1)*sqrt(GSM.6)*alpha[1]*alpha[6]),
                                                           
                                                           (2*cor.SMPWP2.3*sqrt(GSM.2)*sqrt(GSM.3)*alpha[2]*alpha[3]),
                                                           (2*cor.SMPWP2.4*sqrt(GSM.2)*sqrt(GSM.4)*alpha[2]*alpha[4]),
                                                           (2*cor.SMPWP2.5*sqrt(GSM.2)*sqrt(GSM.5)*alpha[2]*alpha[5]),
                                                           (2*cor.SMPWP2.6*sqrt(GSM.2)*sqrt(GSM.6)*alpha[2]*alpha[6]),
                                                           
                                                           (2*cor.SMPWP3.4*sqrt(GSM.3)*sqrt(GSM.4)*alpha[3]*alpha[4]),
                                                           (2*cor.SMPWP3.5*sqrt(GSM.3)*sqrt(GSM.5)*alpha[3]*alpha[5]),
                                                           (2*cor.SMPWP3.6*sqrt(GSM.3)*sqrt(GSM.6)*alpha[3]*alpha[6]),
                                                           
                                                           (2*cor.SMPWP4.5*sqrt(GSM.4)*sqrt(GSM.5)*alpha[4]*alpha[5]),
                                                           (2*cor.SMPWP4.6*sqrt(GSM.4)*sqrt(GSM.6)*alpha[4]*alpha[6]),
                                                           
                                                           (2*cor.SMPWP5.6*sqrt(GSM.5)*sqrt(GSM.6)*alpha[5]*alpha[6]), na.rm=TRUE)
                
            } else {
                
                print(paste (Hrow, "error with weights for", bdd.validation.s[Hrow,]$id_tot1 ))
            }
        }
    }

rm(GSM, GSM.1, GSM.2, GSM.3, GSM.4, GSM.5, GSM.6, Hrow, i, lim.GSM.inf, lim.GSM.sup, t,var.pred, alpha, GSM.preds)



# ### Calculate PICP ------------------------------------------------------

### calculate upper and lower prediction interval limits
bdd.validation.s$SMFC.pred_UPL <- bdd.validation.s$SMFC.pred + (1.64 * sqrt(bdd.validation.s$SMFC.pred.var))
bdd.validation.s$SMFC.pred_LPL <- bdd.validation.s$SMFC.pred - (1.64 * sqrt(bdd.validation.s$SMFC.pred.var))

bdd.validation.s$SMPWP.pred_UPL <- bdd.validation.s$SMPWP.pred + (1.64 * sqrt(bdd.validation.s$SMPWP.pred.var))
bdd.validation.s$SMPWP.pred_LPL <- bdd.validation.s$SMPWP.pred - (1.64 * sqrt(bdd.validation.s$SMPWP.pred.var))

### calculate PICP
PICP <- function(table,var,UPL,LPL) {
    bMat <- c(rep(NA, nrow(table)))
    bMat <- as.numeric(table[,var]<= table[,UPL] & table[,var]>= table[,LPL])
    picp <- sum(bMat)/length(bMat) 
}

picp_SMFC  <- PICP(table = bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),],var = "SMpF2.0", UPL = "SMFC.pred_UPL",  LPL = "SMFC.pred_LPL" )*100
picp_SMPWP <- PICP(table = bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),],var = "SMpF4.2", UPL = "SMPWP.pred_UPL", LPL = "SMPWP.pred_LPL" )*100

####################################################################################################################################################

### Plot something

### AWC
observed_w20 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMpF2.0 
predicted_w20 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMFC.pred
UPL_w20 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMFC.pred_UPL
LPL_w20 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),]$SMFC.pred_LPL

observed_w42 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMpF4.2 
predicted_w42 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMPWP.pred
UPL_w42 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMPWP.pred_UPL
LPL_w42 <- bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),]$SMPWP.pred_LPL

# ### Apply PTF on observed PSD data --------------------------------------
#names(bdd.validation.s)
### Copy dataframe
bdd.validation.PTF <- bdd.validation.s
### load PTF
#load("D:/romandobarco/AWC/SOLHYDRO/clean_output/1.PTF_SOLHYDRO/ContinuousPTF.RData")
### apply PTF 
bdd.validation.PTF$smFC.ptf <- predict(object=lm_w20_ClSd, newdata = bdd.validation.PTF,interval = "prediction", level=0.90 )
bdd.validation.PTF$smPWP.ptf <- predict(object=lm_w42_ClSd, newdata = bdd.validation.PTF,interval = "prediction", level=0.90)


#par(mar=c(5,5,1,1))
#jpeg(filename = "SMFC_Obs_Pred.jpg", width=7, height = 7, units = "in", res=400)

par(mfrow=c(2,2), oma = c(5,5,2,2), mar = c(12,12,10,10), las=1,mgp=c(0, 4, 0))
jpeg("Figure5.jpeg", width = 3740 , height = 3600,units="px")
par(mfrow=c(2,2), oma = c(5,5,2,2), mar = c(12,12,10,10), las=1,mgp=c(0, 4, 0))
#goof(observed = observed_w20, predicted = predicted_w20, type="DSM",plot.it = T)
plot(observed_w20, predicted_w20, col="black",las=1, 
     cex=6, pch=17, ylim=c(0.1,0.5),cex.lab=3,cex.axis=6,cex.main=2,
     xlab="", ylab="")
segments(x0 =observed_w20, y0 =LPL_w20, x1 = observed_w20, y1=UPL_w20, col="cyan3", lwd=3  )
points(observed_w20, predicted_w20, col="black", cex=5, pch=17); abline(0,1,col="black")
mtext(expression("Measured SMFC (cm"^3*"/cm"^3*")"), side=1, line=12, cex=5)
mtext(expression("Predicted SMFC (cm"^3*"/cm"^3*")"), side=2, line=12, cex=5, las=0)
legend("topleft", "a)", bty="n", cex=6) 

#goof(observed = observed_w42, predicted = predicted_w42, type="DSM",plot.it = T)
plot(observed_w42, predicted_w42, col="black",las=1,
     cex=6, pch=17, ylim=c(0,0.44),cex.lab=3,cex.axis=6,cex.main=2,
     xlab="", ylab="")
segments(x0 =observed_w42, y0 =LPL_w42,x1 = observed_w42, y1=UPL_w42, col="cyan3", lwd=3   )
points(observed_w42, predicted_w42, col="black", cex=5, pch=17); abline(0,1,col="black")
mtext(expression("Measured SMPWP (cm"^3*"/cm"^3*")"), side=1, line=12, cex=5)
mtext(expression("Predicted SMPWP (cm"^3*"/cm"^3*")"), side=2, line=12, cex=5, las=0)
legend("topleft", "b)", bty="n", cex=6 ) 
#rm(LPL_w20, LPL_w42, UPL_w42, UPL_w20, observed_w20, observed_w42, predicted_w20, predicted_w42)
###############################################################################################################################################################
#save.image("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/8.2-Validation.preds.GEVARNOVIA/8.2-Validation.preds.GEVARNOVIA.RData")
#################################################################################################################################
### Make the plots
#par(mfrow=c(2,2), las=1, mar=c(5,5,2,1))
plot(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$SMpF2.0,
     bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$smFC.ptf[,1],
     col="black",las=1, cex=6, pch=17, ylim=c(0,0.5),cex.lab=3,cex.axis=6,cex.main=2,
     xlab="", ylab="")
segments(x0 =bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$SMpF2.0,
         y0 =bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$smFC.ptf[,2],
         x1 = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$SMpF2.0,
         y1=bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$smFC.ptf[,3],
         col="cyan3", lwd=3   )
points(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$SMpF2.0,
       bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$smFC.ptf[,1], col="black", cex=5, pch=17); abline(0,1,col="black")
mtext(expression("Measured SMFC (cm"^3*"/cm"^3*")"), side=1, line=12, cex=5)
mtext(expression("Predicted SMFC (cm"^3*"/cm"^3*")"), side=2, line=12, cex=5, las=0)
legend("topleft", "c)", bty="n", cex=6 ) 

plot(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$SMpF4.2,
     bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$smPWP.ptf[,1],
     col="black",las=1, cex=6, pch=17, ylim=c(0,0.44),cex.lab=3,cex.axis=6,cex.main=2,
     xlab="",
     ylab="");
segments(x0 =bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$SMpF4.2,
         y0 =bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$smPWP.ptf[,2],
         x1 = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$SMpF4.2,
         y1=bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$smPWP.ptf[,3],
         col="cyan3", lwd=3  )
points(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$SMpF4.2,
       bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$smPWP.ptf[,1],
       col="black", cex=5, pch=17); abline(0,1,col="black")
mtext(expression("Measured SMPWP (cm"^3*"/cm"^3*")"), side=1, line=12, cex=5)
mtext(expression("Predicted SMPWP (cm"^3*"/cm"^3*")"), side=2, line=12, cex=5 ,las=0)
legend("topleft", "d)", bty="n", cex=6 ) 
dev.off()

### Calculate validation statistics

### calculate upper and lower prediction interval limits
bdd.validation.PTF$SMFC.ptf.pred <- bdd.validation.PTF$smFC.ptf[,1]
bdd.validation.PTF$SMFC.ptf_UPL <- bdd.validation.PTF$smFC.ptf[,3]
bdd.validation.PTF$SMFC.ptf_LPL <- bdd.validation.PTF$smFC.ptf[,2]

bdd.validation.PTF$SMPWP.ptf.pred <- bdd.validation.PTF$smPWP.ptf[,1]
bdd.validation.PTF$SMPWP.ptf_UPL <- bdd.validation.PTF$smPWP.ptf[,3]
bdd.validation.PTF$SMPWP.ptf_LPL <- bdd.validation.PTF$smPWP.ptf[,2]

### Validation statistics
SMFC_valid.ptf <-goof(observed = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$SMpF2.0,
     predicted = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),]$SMFC.ptf.pred, type="DSM",plot.it = T)

SMPWP_valid.ptf <-goof(observed = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$SMpF4.2,
     predicted = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),]$SMPWP.ptf.pred, type="DSM",plot.it = T)

### calculate PICP
PICP <- function(table,var,UPL,LPL) {
    bMat <- c(rep(NA, nrow(table)))
    bMat <- as.numeric(table[,var]<= table[,UPL] & table[,var]>= table[,LPL])
    picp <- sum(bMat)/length(bMat) 
}

picp_SMFC.ptf  <- PICP(table = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),],var = "SMpF2.0", UPL = "SMFC.ptf_UPL",  LPL = "SMFC.ptf_LPL" )*100
picp_SMPWP.ptf <- PICP(table = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),],var = "SMpF4.2", UPL = "SMPWP.ptf_UPL", LPL = "SMPWP.ptf_LPL" )*100

#############################################################################################################################################

### Clean image and save
rm("lm_w20_ClSdBd" , "lm_w20_ClSdBdSOC", "lm_w20_ClSdSOC"  ,      "lm_w20_sub_ClSd" ,
   "lm_w20_sub_ClSdBd"  ,   "lm_w20_sub_ClSdBdSOC" , "lm_w20_sub_ClSdSOC"  ,  "lm_w20_top_ClSd"   ,   
    "lm_w20_top_ClSdBd"  ,   "lm_w20_top_ClSdBdSOC",  "lm_w20_top_ClSdSOC"  ,  "lm_w25_ClSd"     , 
   "lm_w25_ClSdBd" ,        "lm_w25_ClSdBdSOC" ,"lm_w25_ClSdSOC",        "lm_w25_sub_ClSd" ,   
   "lm_w25_sub_ClSdBd"  ,   "lm_w25_sub_ClSdBdSOC" , "lm_w25_sub_ClSdSOC" ,   "lm_w25_top_ClSd" ,     
    "lm_w25_top_ClSdBd"   ,  "lm_w25_top_ClSdBdSOC" , "lm_w25_top_ClSdSOC" ,     "lm_w42_ClSdBd" ,
   "lm_w42_ClSdBdSOC" ,  "lm_w42_ClSdSOC"    ,    "lm_w42_sub_ClSd"  ,  "lm_w42_sub_ClSdBd" , 
   "lm_w42_sub_ClSdBdSOC" , "lm_w42_sub_ClSdSOC"  ,  "lm_w42_top_ClSd" ,     
    "lm_w42_top_ClSdBd"  ,   "lm_w42_top_ClSdBdSOC" , "lm_w42_top_ClSdSOC", INDEP.bdd , GSM.preds,pred_at_arvalis, bdd.validation.s.pwp,
   "SMFC_0_5" ,  "SMFC_100_200" ,    "SMFC_15_30" , "SMFC_30_60" ,  "SMFC_5_15" , "SMFC_60_100" ,
   "SMPWP_0_5", "SMPWP_100_200" , "SMPWP_15_30" , "SMPWP_30_60" ,  "SMPWP_5_15"  , "SMPWP_60_100" ,         
    "Total_var_w20.0_5" ,  "Total_var_w20.100_200", "Total_var_w20.15_30"  , "Total_var_w20.30_60",  "Total_var_w20.5_15" ,  
   "Total_var_w20.60_100",  "Total_var_w42.0_5" ,    "Total_var_w42.100_200" ,"Total_var_w42.15_30"  , "Total_var_w42.30_60",   "Total_var_w42.5_15",   
    "Total_var_w42.60_100")

save.image("9-Validation.preds.GEVARNOVIA.RData")

### End of script