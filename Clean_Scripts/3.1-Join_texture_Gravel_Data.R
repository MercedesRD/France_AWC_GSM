#############################################################################################################################################
###  Prediction of available water capacity in France - GSM depth intervals
###
###  In this script: IGCS data on particle size distribution for France, GSM depth intervals
######## 1. Join texture and coarse elements data - eliminate RMQS profiles
######## 2.Transformation of dataframe, adding different depths as new columns, one row by profile
######## 3. Create raster stack with candidate, predictor covariates
######## 4. Extraction of covariate values at IGCS profiles coordinates

####### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(raster)
library(lattice)
library(spatstat)
library(automap)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RODBC)
library(geoR)


### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

######## Load coarse elements data
load(paste0(HomeDir,"Input/eg_cor2.RData"))

### Link the actual coordinates by id_profil

## Change coordinates in the eg IGCS data
login <- "mromandobar"
password <- "u,Blyd7r"
ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)
coord.wgs84 <- sqlQuery(ds3,"
                        SELECT id_profil,x(the_geom_wgs84),y(the_geom_wgs84)
                        FROM data.poi")

# Select coordinantes data with y > 90
coord.wgs84.errors <- subset(coord.wgs84, y <= -90) ## 29 observations
dim(coord.wgs84.errors); rm(coord.wgs84.errors)
coord.wgs84 <- subset(coord.wgs84, y > -90)

# Define coordiante system (WGS84, EPSG = 4326)
coordinates(coord.wgs84) <- ~ x + y
proj4string(coord.wgs84) <- CRS("+init=epsg:4326")

# Projection from WGS84 (EPSG = 4326) to Lambert 93 (EPSG = 2154)
coord.L93 <- spTransform(coord.wgs84, CRS("+init=epsg:2154"))
## transform into dataframe
coord.L93.df <- as.data.frame(coord.L93)
rm(coord.wgs84, coord.L93)

### substitute in eg_cor2
eg_cor2.df <- merge(eg_cor2,coord.L93.df, by.x="ID", by.y="id_profil", all.x=TRUE, all.y=FALSE)
eg_cor2.df <- eg_cor2.df[,c("ID", "x","y", "l0", "l5", "l15", "l30", "l60", "l100")]

dim(eg_cor2.df[!(is.na(eg_cor2.df$x)| is.na(eg_cor2.df$y)),])  ### Some profiles miss the coordinates now
dim(eg_cor2.df[(is.na(eg_cor2.df$x)| is.na(eg_cor2.df$y)),])

eg_cor2.df <- eg_cor2.df[!(is.na(eg_cor2.df$x)| is.na(eg_cor2.df$y)),]
eg_cor2.df <- eg_cor2.df[eg_cor2.df$y > 4442001, ]

### Transformation into spatial
eg_cor2.sp <- eg_cor2.df 
coordinates(eg_cor2.sp) <- ~ x + y
proj4string(eg_cor2.sp) <- CRS("+init=epsg:2154")
plot(eg_cor2.sp, pch=20, col="darkcyan")

departments<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "DEPARTEMENT")
departments <- spTransform(departments, CRS("+init=epsg:2154"))
Fr.ggmap <- fortify(departments, region="CODE_DEPT") # pour le mettre au format data.frame

### Plot as "Initial_EG.jpg"
ggplot(Fr.ggmap, aes(long, lat)) + geom_polygon(aes(group=group), fill=NA, colour="black", size=.3)+
    geom_point(data=eg_cor2.df, aes(x=x, y=y),colour="coral",size=1) + theme_bw()+coord_equal()+
    ggtitle("Initial EG data")


### there are RMQS observations

### Extract RMQS profiles form Donesol3

# Identification des données rmqs (type de profil)
# type de profil :
# B : composite Biosoil
# C : Composite
# D : prélèvement volumétrique
# F : Fosse
# M : 2ème composite
# N : 3ème composite
# O : 4ème composite
# S : sondage
# T : théorique
# X : 2ème sondage
# Y : 2ème fosse

login <- "mromandobar"
password <- "u,Blyd7r"
ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)

prof.rmqs <- sqlQuery(ds3,"SELECT id_profil, type_profil_rmqs
                      FROM data_rmqs.l_profil_intervention") 

length(unique(prof.rmqs$id_profil))

coord.wgs84 <- sqlQuery(ds3,"SELECT id_profil,x(the_geom_wgs84),y(the_geom_wgs84)
                        FROM data.poi")

### Attach the coordinates to the RMQS sites
prof.rmqs <- merge(prof.rmqs, coord.wgs84, by="id_profil", all.x=TRUE, all.y=FALSE)
### Exclude theoretical profiles
prof.rmqs <- prof.rmqs[prof.rmqs$type_profil_rmqs != "T", ]

### which ones do not have coordinates?
sum(is.na(prof.rmqs$x)) ### 1075 points do not have coordinates in the table 22/11/2017

# Select coordinantes data with y > 90
coord.wgs84.errors <- subset(prof.rmqs, y <= -90) ## 29 observations
dim(coord.wgs84.errors); rm(coord.wgs84.errors)
prof.rmqs <- subset(prof.rmqs, y > -90)
ggplot(prof.rmqs,aes(x,y)) + geom_point(size=0.1) +  coord_equal() ### there are sites from Guyane, etc.

# Define coordiante system (WGS84, EPSG = 4326)
coordinates(prof.rmqs) <- ~ x + y
proj4string(prof.rmqs) <- CRS("+init=epsg:4326")
plot(prof.rmqs)

# Projection from WGS84 (EPSG = 4326) to Lambert 93 (EPSG = 2154)
prof.rmqs.L93 <- spTransform(prof.rmqs, CRS("+init=epsg:2154"))

## eliminate coords in wgs84
rm(prof.rmqs, coord.wgs84)

### superimpose both
plot(departments, main= "Coarse elements (%) & RMQS locations")
points(eg_cor2.sp, pch=1, cex=0.7, col="coral")
points(prof.rmqs.L93, pch=3, cex=0.5,col="blue")
### save as "EG_data_RMQS_sites.jpg"

### To dataframe
prof.rmqs <- as.data.frame(prof.rmqs.L93)
### which of these are in eg_cor?
length(setdiff(prof.rmqs$id_profil, eg_cor2.df$ID)) ### 4252 are not
length(setdiff(eg_cor2.df$ID, prof.rmqs$id_profil)) ### 64424 are not

#### Eliminate RMQS profiles
### One last thing. there are some IGCS profiles that are actually RMQS profiles.

### subset those that are in rmqs.prof
rmqs_id_profil <- prof.rmqs$id_profil
eg_cor2.RMQS <- eg_cor2.df[eg_cor2.df$ID %in% rmqs_id_profil,]  #### 1592 observations from RMQS profiles
points(eg_cor2.RMQS$x,   eg_cor2.RMQS$y ,pch="+", col=2)

### Plot as "Initial_EG_RMQS_ID.jpg"
ggplot(Fr.ggmap, aes(long, lat)) + geom_polygon(aes(group=group), fill=NA, colour="black", size=.3)+
    geom_point(data=eg_cor2.RMQS, aes(x=x, y=y),colour="red",size=1.5) + theme_bw()+coord_equal()+
    ggtitle("EG data with RMQS profile ID")

# add a flag for the origin of the profil
eg_cor2.df$src <- "IGCS"
eg_cor2.df[eg_cor2.df$ID %in% rmqs_id_profil,]$src <- "RMQS"

# Plot as "EG_data_by_source"
ggplot(eg_cor2.df,aes(x,y)) + geom_point() + facet_grid(~src) + coord_equal()

### Yet, there are some points that are still RMQS withing the dataset (they are on a regular GRID).
### It is good to eliminate spatial duplicates

### eliminate spatial duplicates of RMQs locations with a buffer
tmp <- gBuffer(prof.rmqs.L93,width = 25)
join <- over(eg_cor2.sp,tmp)
eg.igcs <- eg_cor2.sp[is.na(join),] ### Exclude points in the buffer area of rmqs profiles

# consider then internal double

ids <- dup.coords(
    as.geodata(obj = as.data.frame(eg.igcs),
               coords = which(colnames(as.data.frame(eg.igcs)) %in% c("x","y")) ,
               data.col  = which(colnames(as.data.frame(eg.igcs)) == c("ID","l0","l5","l15","l30","l60","l100"))  
    )
)

ggplot(Fr.ggmap, aes(long, lat)) + geom_polygon(aes(group=group), fill=NA, colour="black", size=.3)+
    geom_point(data=ids, aes(x=x, y=y),colour="red",size=1.5) + theme_bw()+coord_equal()+
    ggtitle("Spatial duplicate points")

# extract the id of the duplicates
ids$id <- as.numeric(row.names(ids))
# select on duplicate
listeProf <- aggregate(id~dup,ids, FUN = max)
mask <- as.numeric( row.names(ids[!ids$id %in% listeProf$id,]) )
eg.igcs.2 <- eg.igcs[-mask,]
eg.igcs.2.df <- as.data.frame(eg.igcs.2)

ggplot(Fr.ggmap, aes(long, lat)) + geom_polygon(aes(group=group), fill=NA, colour="black", size=.3)+
    geom_point(data=as.data.frame(eg.igcs.2), aes(x=x, y=y),colour="coral",size=1) + theme_bw()+coord_equal()+
    ggtitle("E.G. data clean of RMQS spatial duplicates")

rm(eg_cor2.RMQS,  listeProf, ids,join, mask,  tmp,  eg_cor2.df, eg_cor2.sp, eg.igcs)
rm(coord.L93.df)

# Save backup of this data
setwd(paste0(HomeDir,"Clean_Output/3.1-Join_Texture_Gravel"))
save.image("eg_cor2.igcs.RData")


# ##############################################################################################################################

# 1. Texture data -------------------------------------------------

### IGCS data on texture for the GSM depth intervals

### load the data
load(paste0(HomeDir,"Clean_Output/2-spline_granulo_GSM/granulo_GSM.RData"))

### change name of table
texture_gsm <- data.spl.gsm; rm(data.spl.gsm)


## Change coordinates in the eg IGCS data
login <- "mromandobar"
password <- "u,Blyd7r"
ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)
coord.wgs84 <- sqlQuery(ds3,"
                        SELECT id_profil,x(the_geom_wgs84),y(the_geom_wgs84)
                        FROM data.poi")

# Select coordinantes data with y > 90
coord.wgs84.errors <- subset(coord.wgs84, y <= -90) ## 29 observations
dim(coord.wgs84.errors); rm(coord.wgs84.errors)
coord.wgs84 <- subset(coord.wgs84, y > -90)

# Define coordiante system (WGS84, EPSG = 4326)
coordinates(coord.wgs84) <- ~ x + y
proj4string(coord.wgs84) <- CRS("+init=epsg:4326")

# Projection from WGS84 (EPSG = 4326) to Lambert 93 (EPSG = 2154)
igcs.L93 <- spTransform(coord.wgs84, CRS("+init=epsg:2154"))
## transform into dataframe
igcs.L93.df <- as.data.frame(igcs.L93)
rm(coord.wgs84, igcs.L93)

### substitute in the texture_gsm
texture_gsm <- merge(texture_gsm, igcs.L93.df, by="id_profil", all.x=TRUE, all.y=FALSE)

texture_gsm <- texture_gsm[, c("id_profil","x.y" ,"y.y","prof_sup_moy", "prof_inf_moy",
                                         "argile", "limon","sable","prof_sup_GSM", "prof_inf_GSM")]

colnames(texture_gsm) <- c("id_profil","x" ,"y","prof_sup_moy", "prof_inf_moy",
                                "argile", "limon","sable","prof_sup_GSM", "prof_inf_GSM")

### which of these are in rmqs?
length(setdiff(prof.rmqs$id_profil, texture_gsm$id_profil)) ### [1] 5844 ###NONE of RMQS profile ID is IN texture_GSM
length(setdiff(texture_gsm$id_profil, prof.rmqs$id_profil))

### which of these are in igcs.eg?
length(setdiff(eg.igcs.2.df$ID, texture_gsm$id_profil)) ### [1] 46556  of EG profiles are not included in the texture data
length(setdiff(texture_gsm$id_profil, eg.igcs.2.df$ID)) ### [1] 21164 of texture profiles are not in EG data
length(unique(texture_gsm$id_profil))

### What is the intersection?
length(intersect(unique(texture_gsm$id_profil),unique(eg.igcs.2.df$ID))) ### [1] 16629 are common in two dataframes

### check if they overlap spatially with RMQS data

### Transform table
str(texture_gsm)
###Split by depths
texture_gsm_0_5 <- texture_gsm[texture_gsm$prof_sup_GSM==0,]
texture_gsm_5_15 <- texture_gsm[texture_gsm$prof_sup_GSM==5,]
texture_gsm_15_30 <- texture_gsm[texture_gsm$prof_sup_GSM==15,]
texture_gsm_30_60 <- texture_gsm[texture_gsm$prof_sup_GSM==30,]
texture_gsm_60_100 <- texture_gsm[texture_gsm$prof_sup_GSM==60,]
texture_gsm_100_200 <- texture_gsm[texture_gsm$prof_sup_GSM==100,]

### change column names
texture_gsm_0_5 <- texture_gsm_0_5[,c("id_profil", "x", "y", "argile", "sable", "limon")]
colnames(texture_gsm_0_5) <- c("id_profil", "x" ,"y", "argile_0_5","sable_0_5", "limon_0_5" )

texture_gsm_5_15 <- texture_gsm_5_15[,c("id_profil", "argile", "sable", "limon", "x", "y")]
colnames(texture_gsm_5_15) <- c("id_profil", "argile_5_15","sable_5_15", "limon_5_15", "x" ,"y" )

texture_gsm_15_30 <- texture_gsm_15_30[,c("id_profil", "argile", "sable", "limon", "x", "y")]
colnames(texture_gsm_15_30) <- c("id_profil", "argile_15_30","sable_15_30", "limon_15_30", "x" ,"y" )

texture_gsm_30_60 <- texture_gsm_30_60[,c("id_profil", "argile", "sable", "limon", "x", "y")]
colnames(texture_gsm_30_60) <- c("id_profil", "argile_30_60","sable_30_60", "limon_30_60", "x" ,"y" )

texture_gsm_60_100 <- texture_gsm_60_100[,c("id_profil", "argile", "sable", "limon", "x", "y")]
colnames(texture_gsm_60_100) <- c("id_profil", "argile_60_100","sable_60_100", "limon_60_100", "x" ,"y" )

texture_gsm_100_200 <- texture_gsm_100_200[,c("id_profil", "argile", "sable", "limon", "x", "y")]
colnames(texture_gsm_100_200) <- c("id_profil", "argile_100_200","sable_100_200", "limon_100_200", "x" ,"y" )

### Merge all dataframes
texture_gsm.m <- merge(texture_gsm_0_5, texture_gsm_5_15[,1:4], by="id_profil", all=TRUE, suffixes=c("", ".x") )
texture_gsm.m <- merge(texture_gsm.m, texture_gsm_15_30[,1:4], by="id_profil", all=TRUE, suffixes=c("", ".x") )
texture_gsm.m <- merge(texture_gsm.m, texture_gsm_30_60[,1:4], by="id_profil", all=TRUE, suffixes=c("", ".x") )
texture_gsm.m <- merge(texture_gsm.m, texture_gsm_60_100[,1:4], by="id_profil", all=TRUE, suffixes=c("", ".x") )
texture_gsm.m <- merge(texture_gsm.m, texture_gsm_100_200[,1:4], by="id_profil", all=TRUE, suffixes=c("", ".x") )

rm(texture_gsm_0_5,texture_gsm_5_15, texture_gsm_15_30, texture_gsm_30_60, texture_gsm_60_100, texture_gsm_100_200 )
str(texture_gsm.m)

### Check again for RMQS sites in the texture data, just in case
plot(texture_gsm.m$x , texture_gsm.m$y ,pch="+", col="grey", main="Texture observations")

rmqs_id_profil <- prof.rmqs$id_profil
### subset those that are in the RMQS
texture_gsm.m.RMQS <- texture_gsm.m[texture_gsm.m$id_profil %in% rmqs_id_profil,]  #### 0 observations from RMQS

# add a flag for the origin of the profil
texture_gsm.m$src <- "IGCS"
texture_gsm.m[texture_gsm.m$ID %in% rmqs_id_profil]$src <- "RMQS"
ggplot(texture_gsm.m,aes(x,y)) + geom_point() + facet_grid(~src) + coord_equal()

rm(texture_gsm.m.RMQS)
rm(prof.rmqs)
rm(ds3, login, password, texture_gsm, rmqs_id_profil)

######################################################################################################################################

### Transform to spatial
texture_gsm.sp <- texture_gsm.m
coordinates(texture_gsm.sp) <- ~ x + y
proj4string(texture_gsm.sp) <- CRS("+init=epsg:2154")
plot(texture_gsm.sp, cex=0.5)
points(prof.rmqs.L93, col="deeppink3", cex=0.5, pch="o")

### eliminate spatial duplicates of RMQs locations with a buffer
tmp <- gBuffer(prof.rmqs.L93,width = 25)
join <- over(texture_gsm.sp,tmp)
texture_gsm.sp2 <- texture_gsm.sp[is.na(join),] ### Exclude points in the buffer area of rmqs profiles
### Onlyy 51 profiles were within the buffer, I decide to take them all off

# consider then internal double
ids <- dup.coords(
    as.geodata(obj = as.data.frame(texture_gsm.sp2),
               coords = which(colnames(as.data.frame(texture_gsm.sp2)) %in% c("x","y")) ,
               data.col  = which(colnames(as.data.frame(texture_gsm.sp2)) == c("id_profil","argile_0_5","sable_0_5","limon_0_5",
                                                                       "argile_5_15","sable_5_15","limon_5_15",
                                                                       "argile_15_30","sable_15_30","limon_15_30",
                                                                       "argile_30_60","sable_30_60","limon_30_60",
                                                                       "argile_60_100","sable_60_100","limon_60_100",
                                                                       "argile_100_200","sable_100_200","limon_100_200"))  
    )
)
### There were 533 replicated data

### Plot them
departments<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "DEPARTEMENT")
departments <- spTransform(departments, CRS("+init=epsg:2154"))
Fr.ggmap <- fortify(departments, region="CODE_DEPT") # pour le mettre au format data.frame

ggplot(Fr.ggmap, aes(long, lat)) + geom_polygon(aes(group=group), fill=NA, colour="black", size=.3)+
    geom_point(data=ids, aes(x=x, y=y, color = "red"), size=2) + theme_bw()+coord_equal()+
    ggtitle("Granulo data spatial duplicates")

# extract the id of the duplicates
ids$id <- as.numeric(row.names(ids))
# select on duplicate
listeProf <- aggregate(id~dup,ids, FUN = max)
mask <- as.numeric( row.names(ids[!ids$id %in% listeProf$id,]) )
texture_gsm.sp2 <- texture_gsm.sp2[-mask,]

plot(texture_gsm.sp2, pch=20, cex=0.7, main="Texture data clean of spatial duplicates")
#points(eg_cor2.sp[eg_cor2.sp@data$ID %in% rmqs_id_profil,],pch="+", col=2 ,cex=0.7)

rm(ids, listeProf, ds3, join, login, mask, password, tmp)

texture_gsm.df <- as.data.frame(texture_gsm.sp2)
eg.igcs.df <- eg.igcs.2.df

rm(eg.igcs.2.df, eg_cor2.df, eg_cor2.RMQS, coord.L93.df, igcs.L93.df, prof.rmqs.L93,texture_gsm.2.df, eg_cor2.sp, eg.igcs, texture_gsm.sp)

### Change names
eg.igcs.sp <- eg.igcs.2; rm(eg.igcs.2)
texture_gsm.sp <- texture_gsm.sp2 ; rm(texture_gsm.sp2)

### save image
setwd(paste0(HomeDir,"Clean_Output/3.1-Join_Texture_Gravel"))
save.image("igcs.eg_texture.RData")

##############################################################################################################################

#### Here, 22/11/2017

### Continue, Monday
load("igcs.eg_texture.RData")


#########################################################################################################
texture_gsm.df$id_profil <- factor(texture_gsm.df$id_profil)
eg.igcs.df$ID <- factor(eg.igcs.df$ID)

### Common sites
### What is the intersection?
common_sites_ID <- intersect(texture_gsm.df$id_profil,eg.igcs.df$ID) ### [1] 16337 are common in two dataframes

### sites in EG but not in granulo
sites.EG.not.granulo <- setdiff(eg.igcs.df$ID, texture_gsm.df$id_profil) ### [1] 46848  of EG profiles are not included in the texture data
### sites in granulo but not in EG
sites.granulo.not.EG <- setdiff(texture_gsm.df$id_profil, eg.igcs.df$ID) ### 20885 sites are inclluded in granulo but not in EG

### Merge common sites
common_sites <- merge(texture_gsm.df[texture_gsm.df$id_profil %in% common_sites_ID,],
                      eg.igcs.df[eg.igcs.df$ID %in% common_sites_ID, ], by.x = "id_profil", by.y = "ID", all=TRUE, suffixes=c("", ".x"))

length(unique(common_sites$id_profil))
length(common_sites$id_profil)

### eliminate duplicate rows
common_sites[duplicated(common_sites),]
common_sites <- common_sites[!duplicated(common_sites),]

sites.EG <- merge(texture_gsm.df[texture_gsm.df$id_profil %in% sites.EG.not.granulo,],
                      eg.igcs.df[eg.igcs.df$ID %in% sites.EG.not.granulo, ],
                  by.x = "id_profil", by.y = "ID", all=TRUE, suffixes=c("", ".x"))

sites.EG[is.na(sites.EG$x),]$x <- sites.EG[is.na(sites.EG$x),]$x.x
sites.EG[is.na(sites.EG$y),]$y <- sites.EG[is.na(sites.EG$y),]$y.x


sites.granulo <- merge(texture_gsm.df[texture_gsm.df$id_profil %in% sites.granulo.not.EG,],
                  eg.igcs.df[eg.igcs.df$ID %in% sites.granulo.not.EG, ],
                  by.x = "id_profil", by.y = "ID", all=TRUE, suffixes=c("", ".x"))

sites.granulo <- sites.granulo[!duplicated(sites.granulo),]


### Order by id_profil
common_sites <- common_sites[order(common_sites$id_profil),]
sites.granulo <- sites.granulo[order(sites.granulo$id_profil),]
sites.EG <- sites.EG[order(sites.EG$id_profil),]

### bind all of them
soil_data <- rbind(common_sites,sites.granulo, sites.EG )

###Select columns of interest
soil_data<- soil_data[, c("id_profil","x","y","argile_0_5","sable_0_5","limon_0_5","argile_5_15","sable_5_15","limon_5_15",
                          "argile_15_30","sable_15_30","limon_15_30","argile_30_60","sable_30_60","limon_30_60","argile_60_100",
                          "sable_60_100","limon_60_100","argile_100_200","sable_100_200","limon_100_200","l0","l5","l15","l30",
                          "l60","l100")]     
colnames(soil_data) <- c("id_profil","x","y","argile_0_5","sable_0_5","limon_0_5","argile_5_15","sable_5_15","limon_5_15",
                         "argile_15_30","sable_15_30","limon_15_30","argile_30_60","sable_30_60","limon_30_60",
                         "argile_60_100","sable_60_100","limon_60_100","argile_100_200","sable_100_200","limon_100_200",
                         "coarse_0_5","coarse_5_15","coarse_15_30","coarse_30_60","coarse_60_100","coarse_100_200")

rm(common_sites, common_sites_ID, sites.EG, sites.granulo,  sites.granulo.not.EG, sites.EG.not.granulo)

#### Eliminating too close points

########## check for spatial duplicates in texture data
soil_data.sp <- soil_data 
coordinates(soil_data.sp) <- ~ x + y
proj4string(soil_data.sp) <- CRS("+init=epsg:2154")
plot(soil_data.sp)

dim(zerodist(soil_data.sp, zero = 50.0)) ### 3514 points are closer than 50 m
soil_data <- soil_data[-(zerodist(soil_data.sp, zero = 50.0)[,2]),]
soil_data.sp <- soil_data
coordinates(soil_data.sp) =  ~ x+y
proj4string(soil_data.sp) <- CRS("+init=epsg:2154")
zerodist(soil_data.sp, zero = 50.0)

par(mfrow=c(1,1))
ggplot(Fr.ggmap, aes(long, lat)) + geom_polygon(aes(group=group), fill=NA, colour="black", size=.3)+
    geom_point(data=soil_data, aes(x=x, y=y),colour="coral", size=1) + theme_bw()+coord_equal()+
    ggtitle("Texture and coarse fragments data")
save.image("soil_data.RData")

#########################################################################################################################################
length(unique(soil_data[(!is.na(soil_data$argile_0_5)| !is.na(soil_data$argile_5_15)|
           !is.na(soil_data$argile_15_30)| !is.na(soil_data$argile_30_60) | !is.na(soil_data$argile_60_100)|
           !is.na(soil_data$argile_100_200)),]$id_profil))

### end of the script