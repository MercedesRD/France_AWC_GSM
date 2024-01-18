

### Apply PTFs to RMQS horizon data and compare with preds from GSM maps

### Load packages

library(UsingR)
library(xtable)
library(plotrix)
library(ggplot2)
library(Hmisc)
library(rgdal)
library(FactoMineR)
library(maps)
library(maptools)
library(rgeos)
library(RODBC)
library(RPostgreSQL)
library(calibrate)
library(gdata)
library(sp)
library(doBy)
library(gdata)
library(debug)
library(plyr)
library(viridis)
library(soiltexture)
library(sp)
library(raster)
library(rgdal)
library(dplyr)

### Load data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/RMQS_sites/RMQS.validation.RData")

# #### Attach land use (clc_EXPERT) ---------------------------------------
occupation <- read.csv("D:/romandobarco/AWC/AWC_GSM_Dec2017/RMQS_sites/occupations.csv", sep=";", header=TRUE,na.strings = "")
id_profil_all_sites.rmqs <- read.csv("D:/romandobarco/AWC/AWC_GSM_Dec2017/RMQS_sites/id_profil_tous_sites_rmqs.csv", sep=";", header=TRUE,na.strings = "")
id_profil_sites_foret.rmqs <- read.csv("D:/romandobarco/AWC/AWC_GSM_Dec2017/RMQS_sites/id_profil_sites_forestiers.csv", sep=";", header=TRUE,na.strings = "")

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

### Apply PTFs
str(RMQS.data)
### eliminate rows without clay and sand
RMQS.data <- RMQS.data[(!is.na(RMQS.data$clay)) & (!is.na(RMQS.data$abondance_eg)), ]

### Constrain lower horizon limit to 2 m
RMQS.data[RMQS.data$profondeur_hz_inf > 200,]$profondeur_hz_inf <- 200
RMQS.data$hor_thickness <- RMQS.data$profondeur_hz_inf - RMQS.data$profondeur_hz_sup
str(RMQS.data)
hist(RMQS.data$hor_thickness)
RMQS.data[RMQS.data$hor_thickness > 150,]

ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=0.3) +
    geom_point(data=RMQS.data, aes(x_reel, y_reel), size=3, pch=20) + theme_bw() + coord_equal()

### Calculate AWC by horizon
load("D:/romandobarco/AWC/SOLHYDRO/clean_output/1.PTF_SOLHYDRO/ContinuousPTF.RData")
RMQS.data$clay <- RMQS.data$clay/10
RMQS.data$sand <- RMQS.data$sand/10
RMQS.data$silt <- RMQS.data$silt/10
RMQS.data$SMFC_ptf <- predict(object = lm_w20_ClSd, newdata = RMQS.data )
RMQS.data$SMPWP_ptf <- predict(object = lm_w42_ClSd, newdata = RMQS.data )
RMQS.data$awc_mm_ptf <- (RMQS.data$SMFC_ptf - RMQS.data$SMPWP_ptf) * (1-(RMQS.data$abondance_eg/100)) * RMQS.data$hor_thickness * 10
hist(RMQS.data$awc_mm_ptf)

RMQS.data[is.na(RMQS.data$x_reel),]$x_reel <- RMQS.data[is.na(RMQS.data$x_reel),]$x ### fill missing coordinates
RMQS.data[is.na(RMQS.data$y_reel),]$y_reel <- RMQS.data[is.na(RMQS.data$y_reel),]$y

RMQS.data$id_profil.x <- factor(RMQS.data$id_profil.x)
### Calculate total AWC by profile
### Calculate average clay, silt, and sand

RMQS.awc  <- group_by(.data =RMQS.data,id_profil.x  )
RMQS.awc <- summarise(RMQS.awc,
                      x_reel = mean(x_reel),
                      y_reel = mean(y_reel),
                      count = n(),
                      forest = first(forest),
                      awc_mm_ptf= sum(awc_mm_ptf),
                      prof_profil = sum(hor_thickness))
str(RMQS.awc)

library(gridExtra)
library(ggplot2)
library(viridis)
p1 <-ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=0.3) +
    geom_point(data=RMQS.awc, aes(x_reel, y_reel, color=awc_mm_ptf), size=10, pch=20) + theme_bw() + coord_equal()+
    viridis::scale_color_viridis(option = "D",direction = -1)+
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))
  

p2 <- ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=0.3) +
        geom_point(data=RMQS.awc, aes(x_reel, y_reel, color=prof_profil), size=10, pch=20) + theme_bw() + coord_equal()+
    viridis::scale_color_viridis(option = "A",direction = -1)+
    theme(legend.position="bottom") +
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
    

#####################################################################################################

### extract predictions at RMQS sites
### copy
RMQS.awc.sp <- RMQS.awc

### Transform to spatial and assign CRS
coordinates(RMQS.awc.sp) <- ~x_reel +y_reel
proj4string(RMQS.awc.sp) <- "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

### Load raster files
AWCr_mm_0_200 <- raster("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/9.AWC_Total_depth/AWCr_mm_0_200.tif")
AWCr_mm_0_200_sd <- raster("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/9.AWC_Total_depth/AWCr_mm_0_200.sd.tif")
soil_depth_mm <- raster("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/soil_depth_mm.tif")

#preds.s <- stack(AWCr_mm_0_200, AWCr_mm_0_200_sd, soil_depth_mm)


### extract predictions from raster stack
AWCr_mm_0_200.col <- extract(AWCr_mm_0_200, RMQS.awc.sp,method="simple", df=TRUE)
AWCr_mm_0_200_sd.col <- extract(AWCr_mm_0_200_sd, RMQS.awc.sp,method="simple", df=TRUE)
soil_depth_mm.col <- extract(soil_depth_mm, RMQS.awc.sp,method="simple", df=TRUE)

### bind
RMQS.awc <- cbind(RMQS.awc, AWCr_mm_0_200.col)
RMQS.awc <- cbind(RMQS.awc, AWCr_mm_0_200_sd.col)
RMQS.awc <- cbind(RMQS.awc, soil_depth_mm.col)
RMQS.awc$soil_depth_cm <- RMQS.awc$soil_depth_mm/10
RMQS.awc$GSM.AWC.UPL <- RMQS.awc$AWCr_mm_0_200 + (1.64* RMQS.awc$AWCr_mm_0_200.sd)
RMQS.awc$GSM.AWC.LPL <- RMQS.awc$AWCr_mm_0_200 - (1.64* RMQS.awc$AWCr_mm_0_200.sd)

awc.plot <- ggplot(RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,],
                   aes(awc_mm_ptf, AWCr_mm_0_200, color=soil_depth_cm))
p3 <- awc.plot + 
    theme_bw() +
    labs(y = expression("GSM AWC (mm)"), 
         x = expression("PTF AWC (mm)"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Predicted soil depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = awc_mm_ptf, y = GSM.AWC.UPL, xend = awc_mm_ptf, yend = GSM.AWC.LPL),
                 alpha = 4/10,size=2,
                 data = RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,])+
    geom_point(alpha = 8/10, size=9)
p3


soil.depth.plot <- ggplot(RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,],
                   aes(awc_mm_ptf, soil_depth_cm))
p4 <- soil.depth.plot + geom_point(alpha = 9/10, size=6) +
    theme_bw() +
    labs(y = expression("GSM soil depth (cm)"), 
         x = expression("RMQS soil depth (cm)"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"))+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))
p4

grid.arrange(p1,p2,p3,p4,ncol=2)

grob1 <- grobTree(textGrob("a)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob2 <- grobTree(textGrob("b)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob3 <- grobTree(textGrob("c)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob4 <- grobTree(textGrob("d)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))

##########################################################################
# Plot
jpeg("FigureExtra.jpeg", width = 3740 , height = 3600,units="px" )
grid.arrange(p1+ annotation_custom(grob1),p2+ annotation_custom(grob2),
             p3+ annotation_custom(grob3),p4+ annotation_custom(grob4),nrow=2)
dev.off()



library(ithir)
RMQS.awc <- RMQS.awc[complete.cases(RMQS.awc),]
goof(observed =RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,]$awc_mm_ptf, predicted = RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,]$AWCr_mm_0_200, plot.it = F )
goof(observed =RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,]$prof_profil, predicted = RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,]$soil_depth_cm, plot.it = TRUE )

### calculate PICP
PICP <- function(table,var,UPL,LPL) {
    bMat <- c(rep(NA, nrow(table)))
    bMat <- as.numeric(table[,var]<= table[,UPL] & table[,var]>= table[,LPL])
    picp <- sum(bMat)/length(bMat) 
}

picp_AWC <- PICP(table = RMQS.awc[RMQS.awc$AWCr_mm_0_200 > 1,],
                 var = "awc_mm_ptf", UPL = "GSM.AWC.UPL", LPL = "GSM.AWC.LPL" )*100

### end of the script