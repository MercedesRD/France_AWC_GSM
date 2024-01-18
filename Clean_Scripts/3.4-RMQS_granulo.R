#############################################################################################################################################
#############################################################################################################################################
###
###  Prediction of available water capacity in France - GSM depth intervals
###
###  RMQS profile data on particle size distribution
###
###  Preparation of validation data for validation
###
###  Author: Mercedes Roman Dobarco
###  Date: 8/10/2018
###  REAPEATED: 07/06/2018

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
library(soiltexture)


### Request of data in Donesol:

# SELECT  annee, insee_commune, nom_commune, code_dept, no_campagne,  
# id_site, site_officiel, x_reel, y_reel, x_theo, y_theo, 
# type_profil_rmqs, id_profil, no_horizon, profondeur_hz_sup, profondeur_hz_inf, 
# profondeur_profil, id_prelevement_donesol, prof_prelev_sommet, 
# prof_prelev_base, argile, limon_fin, limon_grossier, sable_fin, 
# sable_grossier, prc_masse_tot_eg
# FROM dm_donnees_ponctuelles.synthese_analyses_rmqs_fosses_prelev
# 
# WHERE no_campagne = '1' AND argile is not null AND
# id_site not between 3000 AND 10000 and id_site not in (13012, 13005) AND
# id_site not between 7000 AND 8000 
# 
# order by id_site, no_horizon;


### Set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
dir.create(paste0(HomeDir,"Clean_Output/3.4-RMQS_granulo"))
setwd(paste0(HomeDir,"Clean_Output/3.4-RMQS_granulo"))

### Read the data
granulo_RMQS <- read.csv(paste0(HomeDir,"Input/RMQS_sites/rmqs_granulo_fosses.07062018.csv"))

### Eliminate sites with strange coordinates
granulo_RMQS <- granulo_RMQS[granulo_RMQS$y_reel > 6000000,]

### Order
granulo_RMQS <- orderBy(~ id_profil + no_horizon, granulo_RMQS)

### how many sites do I have?
length(unique(granulo_RMQS$id_profil)) ### 1622
length(unique(granulo_RMQS$id_site)) ### 1622

rmqs.gr.profil <- unique(granulo_RMQS$id_profil)
rmqs.gr.profil <- sort(decreasing = FALSE, rmqs.gr.profil)

departments<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "DEPARTEMENT")
departments <- spTransform(departments, CRS("+init=epsg:2154"))

Fr.ggmap <- fortify(departments, region="CODE_DEPT") # pour le mettre au format data.frame
ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=0.3) +
    geom_point(data=granulo_RMQS, aes(x_reel, y_reel), size=3, colour="darkcyan", pch=20) + theme_bw() + coord_equal()

### there are sites on islands and Corse, that will not have predictions, but for the moment I leave it as it is.
granulo_RMQS$clay <- granulo_RMQS$argile
granulo_RMQS$silt <- granulo_RMQS$limon_fin + granulo_RMQS$limon_grossier
granulo_RMQS$sand <- granulo_RMQS$sable_fin + granulo_RMQS$sable_grossier

ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=0.3) +
    geom_point(data=granulo_RMQS, aes(x_reel, y_reel, color=silt), size=3, pch=20) + theme_bw() + coord_equal()+
    viridis::scale_color_viridis(option = "A",direction = -1)

save.image("RMQS_granulo.RData")

### End of the script