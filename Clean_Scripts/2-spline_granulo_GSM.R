#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
#   Donesol3 data cleaning and pretreatment for the prediction of soil texture fractions and coarse elements, 
#   for the estimation of available water capacity in Metropolitan France
#   Author: Mercedes Roman Dobarco - 2017
#   Reusing a script written by Marine Lacoste - 2013
#   US INFOSOL Orleans
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


# chargement des packages nécessaires
# Load packages
# memory.limit(10000000000000)
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
library(infosolR)
library(doParallel)
library(parallel)


### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

# chargements des fonctions nécessaires
# Load functions
setwd(paste0(HomeDir,"Clean_scripts"))
source("functions_texture.R")
source("spline.lib.r")
rm(test,depth,splin)

## set the output directory in a different folder to store the output
setwd(paste0(HomeDir,"Clean_Output/2-spline_granulo_AWC"))

# Connect to pgadmin database
login <- "mromandobar"
password <- "u,Blyd7r"

# Spline function to determined granulo and coarse elements for the GSM depth intervals
############################################################################################

### Prepare data before doing splines
############################################################################################

# chargement des données sols traitées à partir du script extraction_data_AWC.R
# Load awc data 
load(paste0(HomeDir,"Clean_Output/1-extraction_data_AWC/preparation_data_granulo_4.RData"))

# Check class of the column with id_profil
class(data.awc$id_profil)
data.awc$id_profil <- as.numeric(as.character(data.awc$id_profil))

# Check the consistence between the max depth of granulo data and the soil depth 
max.depth.granulo <- summaryBy(prof_inf_moy ~ id_profil, data.awc, FUN = max) # max depth for each soil profil, according to their horizons
colnames(max.depth.granulo) <- c("id_profil", "max.depth.granulo")

### Rememeber that in the last script I assigned prof_sol_p (when it existed) to prof.sol.est.....

### From the script extraction_data_AWC.R
# # profondeur du sol, raison de l'arrêt du sondage, abondance et distribution des racines dans le profil
# data <- merge(data, prof.sol,  by = "id_profil", all.x = T, all.y = F, suffixes = c("",".x"))
# max.depth.hz <- summaryBy(prof_inf_moy~id_profil, data,FUN=max)  ###For each profile, I extract the maximum lower depth form all its horizons
# data <- merge(data, max.depth.hz, by="id_profil", all.x=T, all.y=F)
# data <- rename.vars(data, "prof_inf_moy.max", "prof.sol.est")
# data$prof.sol.est[!is.na(data$prof_sol_p)] <- data$prof_sol_p[!is.na(data$prof_sol_p)] ## When there is soil depth in the profile description, we assign that one instead. So.... for just 8777 profiles

data.awc <- merge(data.awc, max.depth.granulo, by="id_profil", all.x=T)
liste1 <- ((data.awc$max.depth.granulo <= data.awc$prof_sol_p) & !is.na(data.awc$prof_sol_p)) # soil profil with granulo shallower than soil depth
liste2 <- ((data.awc$max.depth.granulo > data.awc$prof_sol_p) & !is.na(data.awc$prof_sol_p)) # soil profil with granulo deeper than soil depth
liste3 <- is.na(data.awc$prof_sol_p) #  soil profil without soil depth

length(liste1[liste1==TRUE]) ### 14036 observations
length(liste2[liste2==TRUE]) ### 9888 observations
length(liste3[liste3==TRUE]) ### 97533 observations

### Order and check data
data.awc <- orderBy(~ id_profil + prof_sup_moy, data.awc)

### Correction of soil depth
### I keep the shallowest soil depth
data.awc$prof.sol.cor[liste1] <- data.awc$max.depth.granulo[liste1]  
data.awc$prof.sol.cor[liste2] <- data.awc$prof_sol_p[liste2]
data.awc$prof.sol.cor[liste3] <- data.awc$max.depth.granulo[liste3]

dim(data.awc[data.awc$prof.sol.est == 1,]) ### There are 11 observations where prof.sol.est is 1, and 
data.awc[data.awc$id_profil == 153427,] ### this profile has only one horizon, and the depth is of 1 cm.... I prefer to eliminate it

### Save into table
# setwd("D:/romandobarco/AWC/splines/splines_awc/Output")
# write.table(data.awc[data.awc$prof.sol.est == 1,], file ="prof_sol_1.csv", sep=";", col.names = T,row.names = F )
### Eliminate that profil with one horizon and 0-1 depth
data.awc <- data.awc[data.awc$id_profil != 153427,]
### correct the other observations
data.awc[data.awc$prof.sol.est == 1,]$prof.sol.cor <- data.awc[data.awc$prof.sol.est == 1,]$max.depth.granulo

# Select observations where prof_inf_moy <= soil depth
data.awc.2 <- subset(data.awc, prof_inf_moy <= prof.sol.cor)

### Clean workspace
rm(liste1,liste2,liste3,granulo, granulo.hor,granulo.hor.bck,max.depth.granulo)
rm(granulo.hor.till, granulo.hor.till.sp)

data.awc.2 <- orderBy(~ id_profil + prof_sup_moy, data.awc.2) ###
# save(data.awc.2, file="data.awc.2.RData")


###################################################################################################################################################

### Clean some profiles

### Load data
#setwd("D:/romandobarco/AWC/splines/splines_awc/Output")
# load("data.awc.2.RData")

# détection des profils à spliner
# Select soil profil to spline
data.awc.2$tillage[is.na(data.awc.2$tillage)] <- 0 # if no tillage
maxG <- summaryBy(prof_inf_moy~id_profil, data.awc.2, FUN=max) # max depth for soil horizons
minG <- summaryBy(prof_sup_moy~id_profil, data.awc.2, FUN=min) # min depth for soil horizons

str(maxG)
hist(maxG$prof_inf_moy.max, breaks=40)
summary(maxG$prof_inf_moy.max)
str(minG)
hist(minG$prof_sup_moy.min, breaks=40)
summary(minG$prof_sup_moy.min)

### eliminate that profile where prof_sup_moy.min is greater than 0
elim.profiles <- minG[minG$prof_sup_moy.min > 0,]$id_profil  ### 2301 profiles don't start at 0
data.awc.2 <- data.awc.2[data.awc.2$id_profil %nin% elim.profiles,]

## Redo
maxG <- summaryBy(prof_inf_moy~id_profil, data.awc.2, FUN=max) # max depth for soil horizons
minG <- summaryBy(prof_sup_moy~id_profil, data.awc.2, FUN=min) # min depth for soil horizons


####################################################################################################################################################

data.awc.2 <- orderBy(~ id_profil + prof_sup_moy, data.awc.2)

# selection of the variable to spline
nb.p <- summaryBy(argile~id_profil, data.awc.2, FUN=length) # number of horizon per soil profiles
liste.multip <- subset(nb.p, argile.length > 1)$id_profil # soil profiles with several soil horizons
liste.1p <- subset(nb.p, argile.length == 1)$id_profil # soil profiles with only one soil horizon

data.1p.gsm <- subset(data.awc.2, id_profil %in% liste.1p) # soil profiles with only one soil horizon
data.sp.gsm <- subset(data.awc.2, id_profil %in% liste.multip) # soil profiles with several soil horizons

## check
isTRUE(dim(data.1p.gsm)[1] + dim(data.sp.gsm)[1] == dim(data.awc.2)[1])

########################################################################################################################################################

##############################        Soil profiles with several horizons

##############################       Prepare first horizon on tilled soils

data.sp.gsm <- orderBy(~id_profil+prof_sup_moy, data.sp.gsm)
length(unique(data.sp.gsm$id_profil)) ### 32972 profiles

### comment
### The function prep.hz.multi does not keep the information of many columns, and does not deal well with all profiles, even when these are tilled.

# Prepare soil profils with multiple soil horizons --- Adjust tilled horizons
var <- c("argile", "limon", "sable")

### Separate tilled from non-tilled profiles
data.sp.gsm.tilled <- data.sp.gsm[data.sp.gsm$tillage==1,]
data.sp.gsm.non_tilled <- data.sp.gsm[data.sp.gsm$tillage==0,]

### Correcting some profiles
########## Profile id == 59479
data.sp.gsm.tilled[data.sp.gsm.tilled$id_profil==59479,]$prof_inf_moy <- c(7,45,72,102)

########## Profile id == 29201
data.sp.gsm.non_tilled[data.sp.gsm.non_tilled$id_profil ==29201,]  ### this one I want to eliminate
data.sp.gsm.non_tilled <- data.sp.gsm.non_tilled[data.sp.gsm.non_tilled$id_profil != 29201,]

length(unique(data.sp.gsm.tilled$id_profil)) ### 21687
length(unique(data.sp.gsm.non_tilled$id_profil)) ### 11284


### Apply the function that adds two thin layers before and after tilled horizons
data.phsp.gsm.tilled <- prep.hz.multi(table=data.sp.gsm.tilled, champ.id="id_profil", champ.prof.sup="prof_sup_moy", champ.prof.inf="prof_inf_moy", champ.variable=var,
                             champ.n.horizon="no_horizon", champ.till.depth="prof.till.max", champ.occup="tillage", liste.occup = 1)
data.phsp.gsm.tilled <- orderBy(~ id_profil+prof_sup_moy, data.phsp.gsm.tilled)

length(unique(data.phsp.gsm.tilled$id_profil)) ### 21646

setdiff(unique(data.sp.gsm.tilled$id_profil),unique(data.phsp.gsm.tilled$id_profil))
data.sp.gsm.tilled[data.sp.gsm.tilled$id_profil==12192,]

### old organization
# setwd("D:/romandobarco/AWC/splines/splines_awc")
# dir.create("D:/romandobarco/AWC/splines/splines_awc/GSM")
# setwd("D:/romandobarco/AWC/splines/splines_awc/GSM")

rm(splin,depth,elim.profiles,test, data.awc)
save.image(paste0(HomeDir,"Clean_Output/2-spline_granulo_AWC/preparation_splines_GSM.RData"))
load("preparation_splines_GSM.RData")

####################################################################################################################

### check that the function was applied well
tilled_prof <- unique(data.phsp.gsm.tilled$id_profil)

###
data.sp.gsm[data.sp.gsm$id_profil==tilled_prof[[15000]],]
data.phsp.gsm.tilled[data.phsp.gsm.tilled$id_profil==tilled_prof[[15000]],]

rm(tilled_prof)
### Most of the profiles are alright. There might be some (e.g., id_profil == 3698) where the horizon below the proughed horizon has been erased (????)

### Recalculate the id_tot1
data.phsp.gsm.tilled <- subset(data.phsp.gsm.tilled, select=-c(id_tot1))
data.phsp.gsm.tilled <- add.id(data.phsp.gsm.tilled, c("id_profil","no_horizon"), "id_tot1")

### Eliminate many variables and merge with data.sp.gsm by "id_tot1"
data.phsp.gsm.tilled <- data.phsp.gsm.tilled[ ,c("id_profil","argile","limon","sable","no_horizon","prof_sup_moy","prof_inf_moy","tillage","prof.till.max","id_tot1")]

### Merge with the previous dataframe, data.sp.gsm.tilled
data.phsp.gsm.tilled <- merge(data.phsp.gsm.tilled,data.sp.gsm.tilled, by="id_tot1", all.x=TRUE, all.y=FALSE, suffixes = c("",".x"))

### Eliminate double variables
data.phsp.gsm.tilled <- data.phsp.gsm.tilled[, c("id_profil","id_tot","id_tot1","id_prelevement","argile","limon","sable","no_prelevement","prof_sommet",
                                   "prof_base","no_horizon","cpcs_nom","rp_95_nom","rp_2008_nom","prof_sup_moy","prof_inf_moy","type_prof",
                                   "id_tot2","id_tot3","tillage","prof.till.max","x","y","CODE_DEPT","CODE_REG","NOM_DEPT","NOM_REGION",
                                   "INSEE_COM","abondance_eg","abondance_eg_prin","abondance_eg_sec","depth","prof_sol_p","abond_rac_p",
                                   "distrib_rac_p","arret","prof.sol.est","max.depth.granulo", "prof.sol.cor")]

data.phsp.gsm.tilled <- orderBy(~id_profil+prof_sup_moy, data.phsp.gsm.tilled)

### Bind dataframes with non-tilled and tilled profiles
isTRUE(length(unique(data.phsp.gsm.tilled$id_profil)) + length(unique(data.sp.gsm.non_tilled$id_profil)) == length(unique(data.sp.gsm$id_profil)))
data.phsp.gsm <- rbind(data.phsp.gsm.tilled, data.sp.gsm.non_tilled)
data.phsp.gsm <- orderBy(~id_profil+prof_sup_moy, data.phsp.gsm)


##########################################################################################################################################################

###################################### Splines

### correction of problematic profiles
data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)
rm(data.sp.gsm.non_tilled, data.sp.gsm.tilled)
rm(nb.p,minG,maxG)

#### Splined variables
var <- c("argile", "limon", "sable")

################# Problematic profiles
### i = 2985
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[2985],]
data.sp.gsm[data.sp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[2985],]

### eliminate profile
data.phsp.gsm <- data.phsp.gsm[data.phsp.gsm$id_profil!=unique(data.phsp.gsm$id_profil)[2985],]
data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)

### i = 8457
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[8457],]
data.sp.gsm[data.sp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[8457],]  ## I don't see any problem ....

### eliminate profile
data.phsp.gsm <- data.phsp.gsm[data.phsp.gsm$id_profil!=unique(data.phsp.gsm$id_profil)[8457],]
data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)

### i = 8521
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[8521],]
data.sp.gsm[data.sp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[8521],]  ## I don't see any problem ....

### eliminate profile
data.phsp.gsm <- data.phsp.gsm[data.phsp.gsm$id_profil!=unique(data.phsp.gsm$id_profil)[8521],]
data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)

### i = 11009
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[11009],]
data.sp.gsm[data.sp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[11009],]  ## Problem with horizon limits

###Try to correct horizon limits
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[11009],]$prof_sup_moy
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[11009],]$prof_inf_moy
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[11009] & data.phsp.gsm$prof_inf_moy ==9.50e+01,]
data.phsp.gsm <- data.phsp.gsm[!(data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[11009] & data.phsp.gsm$prof_inf_moy == 9.50e+01),]

data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)

### i = 14169
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[14169],]
data.sp.gsm[data.sp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[14169],]  ## Problem with horizon limits

### eliminate profile
data.phsp.gsm <- data.phsp.gsm[data.phsp.gsm$id_profil!=unique(data.phsp.gsm$id_profil)[14169],]
data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)
 
# ### i = 16641
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[16641],]
data.sp.gsm[data.sp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[16641],]  ## Problem with horizon limits
 
### eliminate profile
data.phsp.gsm <- data.phsp.gsm[data.phsp.gsm$id_profil!=unique(data.phsp.gsm$id_profil)[16641],]
data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)
 
### i = 30399
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[30399],]
data.sp.gsm[data.sp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[30399],]  ## Problem with horizon limits
 
### Correct profile horizon limits
data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[30399],]$prof_inf_moy <- data.phsp.gsm[data.phsp.gsm$id_profil==unique(data.phsp.gsm$id_profil)[30399],]$prof_base
data.phsp.gsm.profiles <- unique(data.phsp.gsm$id_profil)
# 

# Spline of the data for the soil profils with more than one soil horizon
#length(unique(data.phsp.gsm$id_profil))
for(i in 1:length(unique(data.phsp.gsm$id_profil))){
    print(i)
    data.phsp.gsm.i  <- subset(data.phsp.gsm, id_profil == unique(data.phsp.gsm$id_profil)[i])
    data.phsp.gsm.i  <- unique(data.phsp.gsm.i)
    data.phsp.gsm.i <- orderBy(~ id_profil+prof_sup_moy, data.phsp.gsm.i)
    
    if(max(data.phsp.gsm.i$prof_inf_moy) > 0 & max(data.phsp.gsm.i$prof_inf_moy) <= 5) std.depth.i <- c(0, max(data.phsp.gsm.i$prof_inf_moy))
    if(max(data.phsp.gsm.i$prof_inf_moy) > 5 & max(data.phsp.gsm.i$prof_inf_moy) <= 15) std.depth.i <- c(0, 5, max(data.phsp.gsm.i$prof_inf_moy))
    if(max(data.phsp.gsm.i$prof_inf_moy) > 15 & max(data.phsp.gsm.i$prof_inf_moy) <= 30) std.depth.i <- c(0, 5, 15, max(data.phsp.gsm.i$prof_inf_moy))
    if(max(data.phsp.gsm.i$prof_inf_moy) > 30 & max(data.phsp.gsm.i$prof_inf_moy) <= 60) std.depth.i <- c(0, 5, 15, 30, max(data.phsp.gsm.i$prof_inf_moy))
    if(max(data.phsp.gsm.i$prof_inf_moy) > 60 & max(data.phsp.gsm.i$prof_inf_moy) <= 100) std.depth.i <- c(0, 5, 15, 30, 60, max(data.phsp.gsm.i$prof_inf_moy))
    if(max(data.phsp.gsm.i$prof_inf_moy) > 100 & max(data.phsp.gsm.i$prof_inf_moy) <= 200) std.depth.i <- c(0, 5, 15, 30, 60, 100,max(data.phsp.gsm.i$prof_inf_moy))
    if(max(data.phsp.gsm.i$prof_inf_moy) >  200) std.depth.i <- c(0,5,15,30,60,100,200)
    
    data.splsp.i <- pred.spl.multiPV(data.phsp.gsm.i, std.depth=std.depth.i, champ.id="id_profil", champ.prof.sup="prof_sup_moy", champ.prof.inf="prof_inf_moy", champ.variables=var)
    data.splsp.i <- orderBy(~prof_sup_moy, data.splsp.i)
    
    # Link with GSM standard depths
    gsm.d <- c(0,5,15,30,60,100,200)
    data.splsp.i$prof_sup_GSM <-  gsm.d[1:nrow( data.splsp.i)]
    data.splsp.i$prof_inf_GSM <-  gsm.d[2:(nrow( data.splsp.i)+1)]
    
    if (i == 1){
        data.splsp.gsm <-  data.splsp.i
    }else{
        data.splsp.gsm <- rbind( data.splsp.gsm, data.splsp.i)
    }
}

data.splsp.gsm <- orderBy(~id_profil+prof_sup_moy,  data.splsp.gsm)

rm(i,gsm.d,data.phsp.gsm.profiles, std.depth.i,data.splsp.i,data.phsp.gsm.i,var, liste.1p, liste.multip, data.phsp.gsm.tilled) 
save.image(paste0(HomeDir,"Clean_Output/2-spline_granulo_AWC/splines_GSM.RData"))

### I left it here on 31/03/2017

### I start again the 18/04/2017
#setwd("D:/romandobarco/AWC/Splines/2-splines_awc/Output/GSM")
load(paste0(HomeDir,"Clean_Output/2-spline_granulo_AWC/splines_GSM.RData"))


############################################################################################################

### reorder
data.splsp.gsm <- orderBy(~id_profil+prof_sup_moy, data.splsp.gsm)
### check NA
data.splsp.gsm <- subset(data.splsp.gsm, !(is.na(argile)|is.na(limon)|is.na(sable))) # check NA

### Negative data?
data.splsp.gsm$argile[data.splsp.gsm$argile < 0 ] <- 0
data.splsp.gsm$limon[data.splsp.gsm$limon < 0 ] <- 0
data.splsp.gsm$sable[data.splsp.gsm$sable < 0 ] <- 0


#########################################################################################################################################################

######################################  soil profiles with only one horizon

# Management of soil profil with only one soil horizon
data.ph1p.gsm <- subset(data.1p.gsm, select = c("id_profil", "no_horizon", "prof_sup_moy", "prof_inf_moy","argile", "limon","sable"))

for(i in 1:length(unique(data.ph1p.gsm$id_profil))){
    print(i)
    data.ph1p.gsm.i  <- subset(data.ph1p.gsm, id_profil == unique(data.ph1p.gsm$id_profil)[i])
    data.ph1p.gsm.0 <- data.ph1p.gsm.i 
    
    if(max(data.ph1p.gsm.i$prof_inf_moy) > 0 & max(data.ph1p.gsm.i$prof_inf_moy) <= 5) std.depth.i <- c(0, max(data.ph1p.gsm.i$prof_inf_moy))
    if(max(data.ph1p.gsm.i$prof_inf_moy) > 5 & max(data.ph1p.gsm.i$prof_inf_moy) <= 15) std.depth.i <- c(0, 5, max(data.ph1p.gsm.i$prof_inf_moy))
    if(max(data.ph1p.gsm.i$prof_inf_moy) > 15 & max(data.ph1p.gsm.i$prof_inf_moy) <= 30) std.depth.i <- c(0, 5, 15, max(data.ph1p.gsm.i$prof_inf_moy))
    if(max(data.ph1p.gsm.i$prof_inf_moy) > 30 & max(data.ph1p.gsm.i$prof_inf_moy) <= 60) std.depth.i <- c(0, 5, 15, 30, max(data.ph1p.gsm.i$prof_inf_moy))
    if(max(data.ph1p.gsm.i$prof_inf_moy) > 60 & max(data.ph1p.gsm.i$prof_inf_moy) <= 100) std.depth.i <- c(0, 5, 15, 30, 60, max(data.ph1p.gsm.i$prof_inf_moy))
    if(max(data.ph1p.gsm.i$prof_inf_moy) > 100 & max(data.ph1p.gsm.i$prof_inf_moy) <= 200) std.depth.i <- c(0, 5, 15, 30, 60, 100,max(data.ph1p.gsm.i$prof_inf_moy))
    
    if(length(std.depth.i)>1){
        for(j in 2:(length(std.depth.i)-1)){
            data.ph1p.gsm.i <- rbind(data.ph1p.gsm.i,data.ph1p.gsm.0)
        }
    }
    
    data.ph1p.gsm.i$prof_sup_moy <- std.depth.i[1:(length(std.depth.i)-1)]
    data.ph1p.gsm.i$prof_inf_moy <- std.depth.i[2:length(std.depth.i)]
    
    gsm.d <- c(0,5,15,30,60,100,200)
    data.ph1p.gsm.i$prof_sup_GSM <-  gsm.d[1:(length(std.depth.i)-1)]
    data.ph1p.gsm.i$prof_inf_GSM <-  gsm.d[2:length(std.depth.i)]
    
    if (i == 1){
        data.spl1p.gsm <- data.ph1p.gsm.i
    }else{
        data.spl1p.gsm <- rbind(data.spl1p.gsm, data.ph1p.gsm.i)
    }
}

rm(i,j,gsm.d,data.ph1p.gsm.i, data.ph1p.gsm.0)

### eliminate column "no_horizon"
data.spl1p.gsm <- data.spl1p.gsm[ ,-2]

### Order
data.spl1p.gsm <- orderBy(~id_profil+prof_sup_moy, data.spl1p.gsm)

# Merge data (soil profil with several soil horizons and soil profiles with only one soil horizon)
data.spl.gsm <- rbind(data.splsp.gsm, data.spl1p.gsm)
data.spl.gsm <- orderBy(~id_profil+prof_sup_moy, data.spl.gsm)

length(unique(data.spl.gsm$id_profil)) ## 37793

### list of profile_id
profils_id.list <- unique(data.spl.gsm$id_profil)

### Add coordinates in L93
coordinates_L93 <- data.awc.2[,c("id_profil", "x", "y")]
length(unique(coordinates_L93$id_profil)) ### 37840

### Delete duplicated rows
champs <- c("id_profil","x","y")
coordinates_L93 <- del.duplicated1(coordinates_L93, champs, "profiles_coord")
### Merge
data.spl.gsm <- merge(data.spl.gsm, coordinates_L93, by="id_profil", all.x=T, all.y=F )

# Save data
save(data.spl.gsm, file = "granulo_GSM.RData")

rm(champs, profils_id.list, std.depth.i)
#### summary of dataframes

### data.awc.2 : Dataframe with all profiles from the script extraction_data_AWC, corrected soil depth (prof.sol.cor),
###              and top of first horizon at 0 cm

### data.1p.gsm : soil profiles with only one soil horizon, subset from data.awc.2

### data.sp.gsm : soil profiles with >1 soil horizons, subset from data.awc.2

### data.phsp.gsm : Correction of ploughed profiles with the function prep.hz.multi, of data.sp.gsm

### data.splsp.gsm : Predictions from the splines for the GSM depth intervals, for the dataframe data.phsp.gsm

### data.spl1p.gsm : Observations for the GSM depth intervals from those profiles with only one horizon described.

### data.spl.gsm : Texture fraction estimates (measured or splined) for the GSM depth intervals, with coordinates in L93 (EPSG : 2154)

### save image
save.image(paste0(HomeDir,"Clean_Output/2-spline_granulo_AWC/preparation_data_granulo_GSM.RData"))

#### end of script