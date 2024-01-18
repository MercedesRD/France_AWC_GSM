#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  RMQS profile data on coarse elements
###
###  Preparation of validation data for validation
###
###  Author: Mercedes Roman Dobarco
###  Date: 8/12/2017

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
library(calibrate)
library(gdata)
library(sp)
library(doBy)
library(gdata)
library(debug)
library(plyr)
library(soiltexture)

### Set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
setwd(paste0(HomeDir,"Clean_Output/3.3-RMQS_EG"))

### Query in pgAdminIII

# "SELECT id_site, id_profil, no_horizon, 
# prof_sup_min, prof_sup_max, prof_sup_moy, prof_inf_min, prof_inf_max, 
# prof_inf_moy, abondance_eg, abondance_eg_prin, abondance_eg_sec,
# nom_eg1_h, nom_eg2_h, taille_eg1_h, taille_eg2_h, acidite_eg1_h,
# carbonate_eg1_h, acidite_eg2_h, carbonate_eg2_h, trans_eg1_h, 
# trans_eg2_h, forme_eg1_h, forme_eg2_h, orient_eg1_h, orient_eg2_h
# 
# FROM data.horizon inner join data_rmqs.l_profil_intervention using (id_profil)
# inner join data_rmqs.intervention using(id_intervention)
# 
# WHERE type_profil_rmqs  = 'F' AND (id_campagne = '1') AND
# id_site not between 3000 AND 10000  AND
# id_site not in (13012, 13005) AND 
# id_site not between 7000 AND 8000 
#                     
# order by id_profil, no_horizon;"

### Read csv
EG_RMQS <- read.csv(paste0(HomeDir,"Input/RMQS_sites/rmqs_coarse_fosses.07062018.csv"))

### how many observations do I have?
dim(EG_RMQS[!is.na(EG_RMQS$abondance_eg),])
dim(EG_RMQS[!is.na(EG_RMQS$abondance_eg_prin),])
dim(EG_RMQS[!is.na(EG_RMQS$abondance_eg_sec),])

### when abondance_eg is missing but not abondance_eg_prin
EG_RMQS[is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg <- 
    EG_RMQS[is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_prin

### when abondance_eg is missing but not abondance_eg_sec
EG_RMQS[is.na(EG_RMQS$abondance_eg) & is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg <- 
    EG_RMQS[is.na(EG_RMQS$abondance_eg) & is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_sec

### Are there negative values? ### No
subset(EG_RMQS, (abondance_eg<0 | abondance_eg_prin<0 | abondance_eg_sec<0))

# remplacement des NA du champ abon_eg_sec par des 0 
EG_RMQS[is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_sec <- 0

### when there is abondance_eg but not eg_prin
EG_RMQS[!is.na(EG_RMQS$abondance_eg) & is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_prin <- 
    EG_RMQS[!is.na(EG_RMQS$abondance_eg) & is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg -
    EG_RMQS[!is.na(EG_RMQS$abondance_eg) & is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_sec

# vérification de la concordance entre les champs abondance_eg, abondance_eg_prin  et abondance_eg_sec
subset(EG_RMQS, (((abondance_eg_prin + abondance_eg_sec) - abondance_eg) > 0)) ### There are two differeing observations.

all.equal(EG_RMQS[!is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg,
          EG_RMQS[!is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_prin +
              EG_RMQS[!is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_sec) 
              
plot(EG_RMQS[!is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg,
          EG_RMQS[!is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_prin +
              EG_RMQS[!is.na(EG_RMQS$abondance_eg) & !is.na(EG_RMQS$abondance_eg_prin) &  !is.na(EG_RMQS$abondance_eg_sec),]$abondance_eg_sec)          

hist(EG_RMQS$abondance_eg)

##################################################################################################################################

### Add coordinates to the data
### SQL query
#                         SELECT id_profil,st_x(the_geom_wgs84),st_y(the_geom_wgs84)
#                         FROM data.poi
#                         WHERE id_profil is not null;

# Select coordiantes data for id_profil in those selected
coord.wgs84 <- read.csv(paste0(HomeDir,"Input/RMQS_sites/rmqs_coords_fosses.07062018.csv"))
coord.wgs84 <- subset(coord.wgs84, id_profil %in% EG_RMQS$id_profil) ### this gives 2178 unique profiles and locations

### however, we have more id_profiles in the IGCS granulo data, 46738
length(unique(EG_RMQS$id_profil))

### which ones do not have coordinates?
setdiff(unique(EG_RMQS$id_profil), unique(coord.wgs84$id_profil))

# Select coordinantes data with y > 90
coord.wgs84.errors <- subset(coord.wgs84, st_y <= -90) ## No observations
dim(coord.wgs84.errors); rm(coord.wgs84.errors)
coord.wgs84 <- subset(coord.wgs84, st_y > -90)

# Manage coordinates systems
id.profil <- coord.wgs84$id_profil
# Define coordiante system (WGS84, EPSG = 4326)
coord.wgs84 <- SpatialPoints(coord.wgs84[,c("st_x","st_y")], CRS("+init=epsg:4326")) 
plot(coord.wgs84)
# Projection from WGS84 (EPSG = 4326) to Lambert 93 (EPSG = 2154)
coord.L93 <- spTransform(coord.wgs84, CRS("+init=epsg:2154"))
coord.L93 <- cbind(id.profil, coord.L93@coords)

coord.L93 <- coord.L93 [order(coord.L93 [,1],decreasing =F),] 
coord.L93 <- as.data.frame(coord.L93)
coord.L93 <- rename.vars(coord.L93, "id.profil", "id_profil")
coord.L93 <- rename.vars(coord.L93, "st_x", "x")
coord.L93 <- rename.vars(coord.L93, "st_y", "y")
setdiff(unique(EG_RMQS$id_profil), unique(coord.L93$id_profil))

### give coordinates to RMQS_EG data
EG_RMQS <- merge(EG_RMQS,coord.L93, by="id_profil", all.x=TRUE, all.y= FALSE, suffix=c("", "bis") )

### Exclude a site without coordinates
EG_RMQS <- EG_RMQS[!is.na(EG_RMQS$x),]
summary(EG_RMQS$y)
EG_RMQS <- EG_RMQS[EG_RMQS$y >4000000,]

####################################################################################################################################

### Plot the data

departments<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "DEPARTEMENT")
departments <- spTransform(departments, CRS("+init=epsg:2154"))

Fr.ggmap <- fortify(departments, region="CODE_DEPT") # pour le mettre au format data.frame
ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=0.3) +
    geom_point(data=EG_RMQS, aes(x, y), size=2, colour="darkcyan", pch=20) + theme_bw() + coord_equal()

ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=0.3) +
    geom_point(data=EG_RMQS, aes(x, y, color=abondance_eg), size=3, pch=20) + theme_bw() + coord_equal()+
    viridis::scale_color_viridis(option = "A",direction = -1)

rm(login, password, ds3, id.profil, coord.wgs84, coord.L93)
save.image("RMQS_EG.RData")
###
### End of the script