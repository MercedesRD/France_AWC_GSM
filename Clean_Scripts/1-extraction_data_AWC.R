#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
#   Processing of the Donesol3 Data for available water capacity mapping for France
#   Extraction of predictor data for the soil pedotransfer functions (sand, silt, and clay) and coarse elements
#   Mercedes Roman Dobarco, February 2017
#   US INFOSOL Orleans
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


# Load package
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
library(plyr)
library(soiltexture)


### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

### The folder with scripts
setwd(paste0(HomeDir,"Clean_scripts"))

# loading bneeded functions
# Load functions
source("functions_texture.R")
# source("spline.lib.r")

## set the output directory in a different folder to store the output
setwd(paste0(HomeDir,"Clean_Output/1-extraction_data_AWC"))

########################################################################################################################
## Processing sand, silt, clay, and coarse elements data
########################################################################################################################

###########################################################################################################################
# Identifiants pour la connection à DONESOL
# ATTENTION : A COMPLETER
login <- "mromandobar"
password <- "u,Blyd7r"

### Compare all horizons and texture data
### Extract texture data
granulo <- extract.data("granulo",login,password)
prel <- extract.data("prelevements",login,password)

##########################################################################################################################

#### For those prelevements with the same ID, calculate the average of all observations

## how many prelevement_id do we find in both tables?
length(unique(granulo$id_prelevement))
dim(granulo)[1] ### we have more observations than prelevements (several prelevement_id)

### order by id_prelevement
granulo <- orderBy(~ id_prelevement, granulo)

### check those prelevements that have more than one observation
prel.dupl <- data.matrix(granulo[,"id_prelevement"],rownames.force=FALSE)
liste.champ.dupl <- duplicated(prel.dupl)
length(liste.champ.dupl[liste.champ.dupl == TRUE])### we have 446 duplicated id_prelevement

### What is their ID?
prel.dupl.id <- granulo[liste.champ.dupl,"id_prelevement"]
###subset
dupl.prel.df <- granulo[granulo$id_prelevement %in% prel.dupl.id, ]
granulo <- granulo[granulo$id_prelevement %nin% prel.dupl.id, ]
write.table(dupl.prel.df, "dupl_prel_id.csv", sep=";", col.names = T,row.names = F)

### they were all prelevements that had either a value for clay or a value for sand
### A possibility is to calculate silt as the difference from 1000 of the sum of the two fractions, 
### after assigning clay to the observation with sand, etc.

argile <- summaryBy(argile ~ id_prelevement, dupl.prel.df, FUN = mean, na.rm=TRUE)
limon <- summaryBy(limon ~ id_prelevement, dupl.prel.df, FUN = mean,na.rm=TRUE)
sable <- summaryBy(sable ~ id_prelevement, dupl.prel.df,FUN = mean,na.rm=TRUE)

dupl.prel.df2 <- merge(argile, limon, by="id_prelevement")
dupl.prel.df2 <- merge(dupl.prel.df2, sable, by="id_prelevement")
dupl.prel.df2$limon.mean <- 1000 - (dupl.prel.df2$argile.mean + dupl.prel.df2$sable.mean)
rm(argile, sable, limon)
colnames(dupl.prel.df2) <- c("id_prelevement", "argile", "limon", "sable")

### Bind to the granulo dataframe
granulo <- rbind(granulo,dupl.prel.df2)
rm(prel.dupl,liste.champ.dupl, prel.dupl.id,dupl.prel.df, dupl.prel.df2 )      
  
##################################################################################################################################################

#### clean granulo data

## Eliminate observations that have missing values in two fractions, or in three fractions
granulo.2 <- subset(granulo, !((is.na(argile) & is.na(limon)) | (is.na(argile) & is.na(sable)) | (is.na(limon) & is.na(sable)) | (is.na(argile) & is.na(limon) & is.na(sable))))

listea <- is.na(granulo.2$argile) & !is.na(granulo.2$limon) & !is.na(granulo.2$sable)
listeb <- !is.na(granulo.2$argile) & is.na(granulo.2$limon) & !is.na(granulo.2$sable)
listec <- !is.na(granulo.2$argile) & !is.na(granulo.2$limon) & is.na(granulo.2$sable)

# granulo.2 <- granulo.2[!(listea),]
# granulo.2 <- granulo.2[!(listeb),]
# granulo.2 <- granulo.2[!(listec),]

length(listea[listea == TRUE]) ### 16 observations
length(listeb[listeb == TRUE]) ### 0 observations
length(listec[listec == TRUE]) ### 18 observations

### calculate the fraction missing from the other two

granulo.2$argile[listea] <- 1000 - (granulo.2$limon[listea] + granulo.2$sable[listea])
granulo.2$limon[listeb] <- 1000 - (granulo.2$argile[listeb] + granulo.2$sable[listeb])
granulo.2$sable[listec] <- 1000 - (granulo.2$argile[listec] + granulo.2$limon[listec])


### Do the three fractions sum 1000?

### Use the soil.texture package  to normalize to 100
granulo.3 <- granulo.2
granulo.3$CLAY <- granulo.3$argile/10
granulo.3$SILT <- granulo.3$limon/10
granulo.3$SAND <- granulo.3$sable/10

granulo.4 <- TT.normalise.sum(tri.data=granulo.3)

### Assign to main dataframe
granulo.2$argile <- granulo.4$CLAY * 10
granulo.2$limon <- granulo.4$SILT * 10
granulo.2$sable <- granulo.4$SAND * 10

rm(listea, listeb, listec);gc()
granulo.2$argile <- round(granulo.2$argile)
granulo.2$limon <- round(granulo.2$limon)
granulo.2$sable <- round(granulo.2$sable)
rm(granulo.3,granulo.4)

### eliminate observations with negative values
granulo.2 <- subset(granulo.2, argile >= 0 &  limon >= 0 & sable >= 0)

#####################################################################################################################

### Merge granulo and prelevement
length(unique(prel$id_prelevement))
length(unique(granulo.2$id_prelevement))
granulo.2 <- merge(granulo.2, prel, by = "id_prelevement", all.x = T, all.y = F)

### join field id_profil and id_horizon
granulo.2 <- add.id(granulo.2, c("id_profil","no_horizon"), "id_tot1")

### Load all horizons, independently if they have or not measurements of texture
hor <- extract.data("horizons",login,password)
### join field id_profil and id_horizon into new column
hor <- add.id(hor, c("id_profil","no_horizon"), "id_tot1")

### compare levels for id_tot1 between table hor and granulo
length(setdiff(granulo.2$id_tot1, hor$id_tot1)) ### 0 observations do have texture data, but are not in horizons table
length(setdiff(hor$id_tot1, granulo.2$id_tot1)) ### 389060 horizon-profile combination do not have texture data 
### March 22, 388942 horizon-profile combination do not have texture data (>???) why different?

# ### check if those horizon-profile combinations belong to RMQS or are IGCS fiction profile
# granulo.tmp <- merge(granulo, hor, by="id_tot1", all.x=T)
# hor.tmp <- hor
# ## which ones belong to RMQS?
# nrow(granulo.tmp[!is.na(granulo.tmp$type_profil_rmqs),]) ### 8975 belong to RMQS program
# nrow(granulo.tmp[granulo.tmp$type_prof=="1",]) ### 1783 are fictional IGCS profile
# nrow(granulo.tmp[(granulo.tmp$type_prof=="1" & !is.na(granulo.tmp$type_profil_rmqs)),]) ## 52 are both
# rm(granulo.tmp)
# 
# ### eliminate RMQS profiles
# nrow(hor[!is.na(hor$type_profil_rmqs),]) ### 18524 observations of the RMQS program
# hor.tmp <- hor.tmp[is.na(hor.tmp$type_profil_rmqs),]
# 
# ### Eliminate non-existing profiles in the IGCS, and all profile types from RMQS campaign
# hor.tmp <- subset(hor.tmp, type %nin% c("1","T", "B", "C", "D", "F", "M", "N", "O", "S","X", "Y")) # "1" = profil fictif pour l'IGCS
# hor.tmp <- orderBy(~ id_profil + prof_sommet, hor)
# hor.tmp <- droplevels(hor)
# hor.tmp <- hor.tmp[,-10:11] ### Now I don't need the variables type profile rmqs and type
# 
# ### compare levels for id_tot1 between table hor and granulo
# length(setdiff(granulo$id_tot1, hor.tmp$id_tot1)) ### 9299 observations do have texture data, but are not in horizons table (maybe they belong to the RMQS or are fictif IGCS profiles)
# length(setdiff(hor.tmp$id_tot1, granulo$id_tot1)) ### 373312 horizon-profile combination do not have texture data (these we have to eliminate)
# rm(hor.tmp)

### Join horizons and granulo tables
granulo.hor <- merge(granulo.2, hor, by="id_tot1", all.x=T, all.y=F)

### eliminate RMQS profiles
nrow(granulo.hor[!is.na(granulo.hor$type_profil_rmqs),]) ### 8975 observations of the RMQS program
granulo.hor <- granulo.hor[is.na(granulo.hor$type_profil_rmqs),]

### Eliminate non-existing profiles in the IGCS, and all profile types from RMQS campaign
granulo.hor <- subset(granulo.hor, type %nin% c("1","T", "B", "C", "D", "F", "M", "N", "O", "S","X", "Y")) # "1" = profil fictif pour l'IGCS

### change names and elimnate some columns that are repeated
granulo.hor <- granulo.hor[,c("id_tot1","id_prelevement","argile","limon","sable","no_prelevement","prof_sommet","prof_base",
                              "id_profil.x","no_horizon.x","cpcs_nom", "rp_95_nom","rp_2008_nom","prof_sup_moy","prof_inf_moy",
                              "type_prof","type_profil_rmqs","type")]

colnames(granulo.hor) <- c("id_tot1","id_prelevement","argile","limon","sable","no_prelevement","prof_sommet","prof_base",
                           "id_profil","no_horizon","cpcs_nom", "rp_95_nom","rp_2008_nom","prof_sup_moy","prof_inf_moy",
                           "type_prof","type_profil_rmqs","type")

granulo.hor <- orderBy(~ id_profil + prof_sommet, granulo.hor) ## we order by id_profil and the upper depth of the prelevements
granulo.hor <- droplevels(granulo.hor)
granulo.hor <- granulo.hor[,-18] ### Now I don't need the variables type profile rmqs and type
granulo.hor <- granulo.hor[,-17]

##############################################################################################################################
### we start with 147483 observations
granulo.hor.bck <- granulo.hor ## copy
# granulo.hor <- granulo.hor.bck

# ### Count number of observations by profile-horizon combination
# granulo.hor2 <- granulo.hor
# granulo.hor2$id_tot1 <- as.factor(granulo.hor2$id_tot1) ## transformed to factor
# 
# # ### are there more than one measurement by horizon-profile?
# library(plyr)
# no_hor_prof <- ddply(.data = granulo.hor2,.variables ="id_tot1",summarise, count=length(argile))
# hist(no_hor_prof$count)
# summary(no_hor_prof$count)
# # ### Most have one, but there are some cases with 2, or even 11 observations
# 
# ### For several horizons, there are several prelevements by horizons, each sampled at different depths (!= prof_sommet and != prof_base)
# ### For what I've seen, these are mostly in very deep horizons (>2 m)
# length(unique(no_hor_prof[no_hor_prof$count>1,]$id_tot1)) ## 1097 profile-horizon have more than 1 observation
# length(unique(no_hor_prof[no_hor_prof$count>2,]$id_tot1)) ## 119 profile-horizon have more than 2 observation
# length(unique(no_hor_prof[no_hor_prof$count>3,]$id_tot1)) ## 38 profile-horizon have more than 3 observation
# length(unique(no_hor_prof[no_hor_prof$count>4,]$id_tot1)) ## 24 profile-horizon have more than 3 observation
# length(unique(no_hor_prof[no_hor_prof$count>5,]$id_tot1)) ## 19 profile-horizon have more than 3 observation
# length(unique(no_hor_prof[no_hor_prof$count>8,]$id_tot1)) ## 5 profile-horizon have more than 8 observation
# length(unique(no_hor_prof[no_hor_prof$count>10,]$id_tot1)) ## 1 profile-horizon have more than 10 observation
# length(unique(no_hor_prof[no_hor_prof$count>7,]$id_tot1)) ## 8 profile-horizon have more than 8 observation
# 
# which.noes <- unique(no_hor_prof[no_hor_prof$count==6,]$id_tot1)
# which.noes <- droplevels(which.noes)
# granulo.hor2[granulo.hor2$id_tot1 %in% which.noes,]
# granulo.hor2[granulo.hor2$id_profil == 66031,]
#rm(no_hor_prof, granulo.hor2)
rm(hor, prel)

###########################################################################################################################################

### Calculate average texture values by horizon-profile
# length(unique(granulo.hor$id_tot1)) 
# argile <- summaryBy(argile ~ id_tot1, granulo.hor, FUN = mean)
# limon <- summaryBy(limon ~ id_tot1, granulo.hor, FUN = mean)
# sable <- summaryBy(sable ~ id_tot1, granulo.hor, FUN = mean)
# granulo <- merge(argile, limon, by = "id_tot1")
# granulo <- merge(granulo, sable, by = "id_tot1")
# granulo <- rename.vars(granulo, "argile.mean","argile")
# granulo <- rename.vars(granulo, "limon.mean","limon")
# granulo <- rename.vars(granulo, "sable.mean","sable")

### Before joining the other table... we need to fix the depth of the horizons! 
### Repeat calculating the average later

###############################################################################################################################################

### FIX horizon depths

# Correction des valeurs 999
# If the fields "prof_base","prof_inf_moy", and "prof_sommet" == 999 --> change to NA

### repeat
# granulo.hor <- granulo.hor.bck

##########################################################################################################################################
champs <- c("prof_base","prof_inf_moy","prof_sup_moy","prof_sommet")
granulo.hor <- cor.999(granulo.hor, champs, "granulo") ### THERE WAS ONLY ONE!

# Horizons avec champs de profondeur sans valeur (prof_sommet ou prof_base) 
# Checking the values of the horizons' limits and correction if necessary

##########################################################################################################################################

### Fix missing values
granulo.hor <- lim.hz.NA.1(granulo.hor,"granulo")

### Despite using cor.999, there is one observation with prof_sommet == 999. I eliminate it
granulo.hor <- subset(granulo.hor, granulo.hor$prof_sommet != 999)

### Is there only one prelevement_id by observation?
length(unique(granulo.hor$id_prelevement)) ### YES, 147408
# no_prel_id <- ddply(.data = granulo.hor,.variables ="id_prelevement",summarise, count=length(argile))
# summary(no_prel_id$count) ### 1
# rm(no_prel_id) 

## Correct when there is prof_base and prof_inf_moy are NA
granulo.hor <- lim.hz.NA.2.2(granulo.hor,"granulo")
## Correct when prof_inf_moy = NA but prof_base != NA
granulo.hor <- lim.hz.NA.3(granulo.hor,"granulo")

#########################################################################################################################
### Eliminate O horizons, where limit < 0

##########################################################################################################################################
granulo.hor <- lim.hz.neg(granulo.hor,"granulo")

# problème dans les numéros d'horizons : présence de valeurs négatives
# serait dû à une mauvaise saisie
a <-subset(granulo.hor, no_horizon < 0)
dim(a)
length(unique(a$id_profil)) #4
unique(a$id_profil) # 80885 80954 80984 81001

granulo.hor  <- subset(granulo.hor , id_profil %nin% c(80885, 80954, 80984, 81001));rm(a)

# corection des identifiants (au cas où)
granulo.hor <- subset(granulo.hor, select=-c(id_tot1))
granulo.hor <- add.id(granulo.hor, c("id_profil","no_horizon"), "id_tot1")

# check that argile, sand, or silt is not > 1000

##########################################################################################################################################
nrow(granulo.hor[granulo.hor$argile > 1000,])
nrow(granulo.hor[granulo.hor$limon > 1000,])
nrow(granulo.hor[granulo.hor$sable > 1000,])

# Horizons avec inversion des limites d'horizons
# Check horizon limits inversion

##########################################################################################################################################
granulo.hor <- coherence.lim(granulo.hor,"granulo")

# supressions des doublons : champs "id_profil","prof_sommet","prof_base","argile" "sable", "limon"
# check if they are duplicated data

##########################################################################################################################################
champs <- c("id_profil","prof_sommet","prof_base","argile" ,"sable", "limon")
granulo.hor <- del.duplicated1(granulo.hor, champs, "granulo")

# Supression des horizons quand les limites des analyses ne sont pas comprises dans les limites des horizons
# Delete horizons if limits of analysis are not comprise in horizons limits

##########################################################################################################################################
granulo.hor <- del.decal(granulo.hor, "granulo")

### Order by id_profil and no_horizon
granulo.hor <- orderBy(~ id_profil + no_horizon, granulo.hor)
granulo.hor <- droplevels(granulo.hor)

##########################################################################################################################################

#### Calculate average texture by horizon, when there are horizons with several prelevements

# library(plyr)
# no_hor_prof <- ddply(.data = granulo.hor,.variables ="id_tot1",summarise, count=length(argile))
# hist(no_hor_prof$count)
# summary(no_hor_prof$count)
# 
# ### Most have one, but there are some cases with 2, or even 11 observations

### Make temporal copy
granulo.hor.tmp <- granulo.hor

### check those prelevements that have more than one observation
prof_hor.dupl.liste <- duplicated(granulo.hor.tmp$id_tot1)
length(prof_hor.dupl.liste[prof_hor.dupl.liste == TRUE]) ### we have 1112 duplicated profile_horizon

### What is their ID?
prof_hor.dupl.id <- granulo.hor.tmp$id_tot1[prof_hor.dupl.liste]
### Identify the profile-horizon unique id
horizon_prof_id <- unique(prof_hor.dupl.id) ### 915 unique IDs

###subset those with several prelevements form horizon, from those with a single prelevement by horizon
dupl.prof_hor.df <- granulo.hor.tmp[granulo.hor.tmp$id_tot1 %in% horizon_prof_id, ]
granulo.hor.tmp  <- granulo.hor.tmp[granulo.hor.tmp$id_tot1 %nin% horizon_prof_id, ]

### Calculate average value by horizon
argile <- summaryBy(argile ~ id_tot1, dupl.prof_hor.df, FUN = mean, na.rm=TRUE)
limon  <- summaryBy(limon ~ id_tot1,  dupl.prof_hor.df, FUN = mean, na.rm=TRUE)
sable  <- summaryBy(sable ~ id_tot1,  dupl.prof_hor.df, FUN = mean, na.rm=TRUE)
granulo.prof_hor_single <- merge(argile, limon, by = "id_tot1")
granulo.prof_hor_single <- merge(granulo.prof_hor_single, sable, by = "id_tot1")

### Merge with the duplicated dataframe
dupl.prof_hor.df <- merge(dupl.prof_hor.df,granulo.prof_hor_single, by="id_tot1" ) ### this step changes the order of columns
rm(argile,limon,sable, prof_hor.dupl.id,prof_hor.dupl.liste,granulo.prof_hor_single, champs)

### susbtitute the prelevement values by the average ones
dupl.prof_hor.df$argile <- dupl.prof_hor.df$argile.mean
dupl.prof_hor.df$limon  <- dupl.prof_hor.df$limon.mean
dupl.prof_hor.df$sable  <- dupl.prof_hor.df$sable.mean
dupl.prof_hor.df <- dupl.prof_hor.df[,1:16]

### Assign same values to lower and upper horizon limits for the same id_tot1
for (i in 1:length(horizon_prof_id)){
    ### subset each profile-horizon combination
    dupl.prof_hor.i <- dupl.prof_hor.df[dupl.prof_hor.df$id_tot1 == horizon_prof_id[[i]],]
    
    ### Ask is the horizon lower depth is the same 
    if(isTRUE(length(unique(dupl.prof_hor.i$prof_inf_moy))==1)) {
        print("Same horizon lower depth")
    } else {
        ##If it is not the same
        print("Different horizon lower depth")
        ### Identify the deepest lower depth
        lower_depth.i <- max(dupl.prof_hor.i$prof_inf_moy)
        ### Assign to the field prof_inf_moy
        dupl.prof_hor.df[dupl.prof_hor.df$id_tot1 == horizon_prof_id[[i]],]$prof_inf_moy <- lower_depth.i
        
    }
    
    ### Ask is the horizon upper depth is the same 
    if(isTRUE(length(unique(dupl.prof_hor.i$prof_sup_moy))==1)) {
        print("Same horizon upper depth")
    } else {
        ##If it is not the same
        print("Different horizon upper depth")
        ### Identify the minimum upper depth
        upper_depth.i <- min(dupl.prof_hor.i$prof_sup_moy)
        ### Assign to the field prof_sup_moy
        dupl.prof_hor.df[dupl.prof_hor.df$id_tot1 == horizon_prof_id[[i]],]$prof_sup_moy <- upper_depth.i
        
    }
    
}
rm(i,dupl.prof_hor.i, upper_depth.i,lower_depth.i);gc()

dupl.prof_hor.df$argile <- round(dupl.prof_hor.df$argile)
dupl.prof_hor.df$limon  <- round(dupl.prof_hor.df$limon)
dupl.prof_hor.df$sable  <- round(dupl.prof_hor.df$sable)

### Eliminate duplicates, based on id_tot1, and texture values
champs <- c("id_profil", "no_horizon","argile", "limon", "sable", "prof_sup_moy", "prof_inf_moy")
dupl.prof_hor.df <- del.duplicated1(dupl.prof_hor.df, champs, "prof_hor")

### bind to the main dataframe
dupl.prof_hor.df <- dupl.prof_hor.df[,c("id_prelevement","argile","limon","sable","no_prelevement","prof_sommet","prof_base",
                                        "id_profil","no_horizon","cpcs_nom", "rp_95_nom","rp_2008_nom","prof_sup_moy","prof_inf_moy",
                                        "type_prof", "id_tot1")] 

granulo.hor <- rbind(granulo.hor.tmp,dupl.prof_hor.df)

### Order by id_profil and no_horizon
granulo.hor <- orderBy(~ id_profil + no_horizon, granulo.hor)
granulo.hor <- droplevels(granulo.hor)
rm(dupl.prof_hor.df, horizon_prof_id, granulo.hor.tmp, champs)

########################################################################################################################

### Are there different measurements that by mistake have assigned different no_horizon for a same profile,
### but share horizon limits? (i.e., they were different prelevements from a same horizon, 
### but were recorded as different horizons instead).

### create new variable
granulo.hor<- add.id(granulo.hor, c("id_profil","prof_sup_moy","prof_inf_moy"),"id_tot")

### how many unique values do we find for it?
length(unique(granulo.hor$id_tot)) ### 140149 ### March 23, 140625
dim(granulo.hor)[1] ### 140292 observations ### March 22, 140768

## Locate duplicated cases
liste.dupl <- duplicated(granulo.hor$id_tot)
length(liste.dupl[liste.dupl == TRUE]) ### 143 duplicated cases 

### subset these observations with duplicated id_tot
prof_hor.dupl.id <- granulo.hor$id_tot[liste.dupl]
prof_hor_id <- unique(prof_hor.dupl.id)

### subset those with several prelevements form horizon, from those with a single prelevement by horizon
dupl.prof_hor.df <- granulo.hor[granulo.hor$id_tot %in% prof_hor_id, ]
granulo.hor.tmp <- granulo.hor[granulo.hor$id_tot %nin% prof_hor_id, ]

write.table(dupl.prof_hor.df, "prof_eq_depth_diff_hor.csv", sep=";", col.names = T,row.names = F)

### It seems that there are different prelevements from a same horizon that have been recorded with a different horizon number

### Calculate average value by horizon
argile <- summaryBy(argile ~ id_tot, dupl.prof_hor.df, FUN = mean, na.rm=TRUE)
limon  <- summaryBy(limon ~ id_tot,  dupl.prof_hor.df, FUN = mean, na.rm=TRUE)
sable  <- summaryBy(sable ~ id_tot,  dupl.prof_hor.df, FUN = mean, na.rm=TRUE)
granulo.prof_depth_single <- merge(argile, limon, by = "id_tot")
granulo.prof_depth_single <- merge(granulo.prof_depth_single, sable, by = "id_tot")

### Merge with the duplicated dataframe
dupl.prof_hor.df <- merge(dupl.prof_hor.df,granulo.prof_depth_single, by="id_tot" ) ### this step changes the order of columns
rm(argile,limon,sable, prof_hor.dupl.id,granulo.prof_depth_single)

### susbtitute the prelevement values by the averaged ones
dupl.prof_hor.df$argile <- dupl.prof_hor.df$argile.mean
dupl.prof_hor.df$limon  <- dupl.prof_hor.df$limon.mean
dupl.prof_hor.df$sable  <- dupl.prof_hor.df$sable.mean
dupl.prof_hor.df        <- dupl.prof_hor.df[,1:17] 
dupl.prof_hor.df$argile <- round(dupl.prof_hor.df$argile)
dupl.prof_hor.df$limon  <- round(dupl.prof_hor.df$limon)
dupl.prof_hor.df$sable  <- round(dupl.prof_hor.df$sable)

### Eliminate duplicates, based on id_tot, and texture values
champs <- c("id_profil", "argile", "limon", "sable", "prof_sup_moy", "prof_inf_moy")
dupl.prof_hor.df <- del.duplicated1(dupl.prof_hor.df, champs, "prof_diff_hor_same_lim")

### bind to the main dataframe

names(granulo.hor.tmp)
# "id_prelevement","argile","limon","sable","no_prelevement","prof_sommet","prof_base",   
# "id_profil","no_horizon","cpcs_nom","rp_95_nom","rp_2008_nom","prof_sup_moy","prof_inf_moy", 
# "type_prof","id_tot1","id_tot"

### change order of variables
dupl.prof_hor.df <- dupl.prof_hor.df[,c("id_prelevement","argile","limon","sable","no_prelevement","prof_sommet","prof_base",   
                                        "id_profil","no_horizon","cpcs_nom","rp_95_nom","rp_2008_nom","prof_sup_moy","prof_inf_moy", 
                                        "type_prof","id_tot1","id_tot")] 

granulo.hor <- rbind(granulo.hor.tmp,dupl.prof_hor.df)

### Order by id_profil and no_horizon
granulo.hor <- orderBy(~ id_profil + no_horizon, granulo.hor)
granulo.hor <- droplevels(granulo.hor)
rm(dupl.prof_hor.df, champs,liste.dupl,prof_hor_id,granulo.hor.tmp)

### more IDs
### create more fields de ID
granulo.hor <- add.id(granulo.hor, c("id_profil","prof_sup_moy"), "id_tot2")
granulo.hor <- add.id(granulo.hor, c("id_profil","prof_inf_moy"), "id_tot3")

# We could ..... supression des horizons partageant une même limite d'horizon : pas possible de trier les horizons une à un pour définir quel est l'horizon à conserver. 
# Plusieurs possibilités : erreur de saisi, présence de glosse ou de zones différenciées au sein d'un horizon : plusieurs analyses pour un mçme horizon.
# Delete horizons if they share one limit

### how many have repeated prof_sup_moy?
liste2 <- granulo.hor$id_tot2[duplicated(granulo.hor$id_tot2)]
dim(granulo.hor[granulo.hor$id_tot2 %in% liste2,]) 
head(granulo.hor[granulo.hor$id_tot2 %in% liste2,])

liste3 <- granulo.hor$id_tot3[duplicated(granulo.hor$id_tot3)]
dim(granulo.hor[granulo.hor$id_tot3 %in% liste3,])
head(granulo.hor[granulo.hor$id_tot3 %in% liste3,]) 

### For texture I prefered to calculate the mean for the same horizon... although... 
### if different horizons (id_tot2 and id_tot3) have the same upper or lower limit, we need to  eliminate some measurements

## I separate those from the main dataset
granulo.hor_del <- subset(granulo.hor, id_tot2 %in% liste2 | id_tot3 %in% liste3) ## 423 observations
### Eliminate those horizons that share prof_sup_moy and/or prof_inf_moy
granulo.hor <- subset(granulo.hor, id_tot2 %nin% liste2 & id_tot3 %nin% liste3)
rm(granulo.hor_del)
# problème d'horizons supperposés mais ne partageant pas de limite de prelevement?
# check horizon superposition
hz.sup <- hz.superpose(granulo.hor, "id_profil", "prof_sup_moy", "prof_inf_moy",c("id_profil", "prof_sup_moy", "prof_inf_moy"), "id_tot")

liste.hz.sup <- hz.sup[[2]]
write.table(liste.hz.sup, "overlapping_hor.csv", sep=";", col.names = T,row.names = F)
granulo.hor <- hz.sup[[1]]

rm(liste1, liste2, liste3, granulo.2,granulo.hor_del)
rm(liste.hz.sup, hz.sup)
save(list =ls(), file = "preparations_donnees_granulo.RData")
write.table(granulo.hor, "granulo_data.csv", sep=";", col.names = T,row.names = F)

########################################################################################################################

###### Continue March 27, 2017.

###### We check land use, and type of crop to assume that tilled soils have homogeneous texture on topsoils

###### Load session image

load("preparations_donnees_granulo.RData")

### Check that there is an id_prelevement by observation

length(unique(granulo.hor$id_prelevement)) ### 139140
dim(granulo.hor) ### [1] 139140     19
## ok


# Correction des horizons de surface

################################################################################################################################
granulo.hor.till <- cor.topsoil3(granulo.hor,login,password)
granulo.hor.till <- orderBy(~ id_profil + no_horizon, granulo.hor.till)    

save(list =ls(), file = "preparations_donnees_granulo_1.RData")
write.table(granulo.hor.till, "granulo_data_till.csv", sep=";", col.names = T,row.names = F)
granulo.hor.till[1:50, c(3,4,5,12,11,10,13,14,16,20,21)]

########################################################################################################################
## Ajout des cordonnées des profils et sélection des profils en France métropolitaine
## Add coordinates and selection of soil profils in France (without island)

########################################################################################################################

# Load coordinate data

#######################################################################################################################
# login <- login
# password <- password

ds3 <- odbcConnect(dsn="Donesol3_ns64",uid=login,pwd=password)
# perso <- odbcConnect(dsn="baracoaPerso",uid=login,pwd=password)
sqlQuery(ds3,"SET client_encoding TO \'UTF8\'") 

coord.wgs84 <- sqlQuery(ds3,"
                        SELECT id_profil,x(the_geom_wgs84),y(the_geom_wgs84)
                        FROM data.poi")

# Select coordiantes data for id_profil in those selected
coord.wgs84 <- subset(coord.wgs84, id_profil %in% granulo.hor.till$id_profil) ### this gives 46255 unique profiles and locations

### however, we have more id_profiles in the IGCS granulo data, 46738
length(unique(granulo.hor.till$id_profil))

### which ones do not have coordinates?
setdiff(unique(granulo.hor$id_profil), unique(coord.wgs84$id_profil))

# Select coordinantes data with y > 90
coord.wgs84.errors <- subset(coord.wgs84, y <= -90) ## No observations
dim(coord.wgs84.errors); rm(coord.wgs84.errors)
coord.wgs84 <- subset(coord.wgs84, y > -90)

# Manage coordinates systems
id.profil <- coord.wgs84$id_profil
# Define coordiante system (WGS84, EPSG = 4326)
coord.wgs84 <- SpatialPoints(coord.wgs84[,c("x","y")], CRS("+init=epsg:4326")) 
plot(coord.wgs84)
# Projection from WGS84 (EPSG = 4326) to Lambert 93 (EPSG = 2154)
coord.L93 <- spTransform(coord.wgs84, CRS("+init=epsg:2154"))
coord.L93 <- cbind(id.profil, coord.L93@coords)

coord.L93 <- coord.L93 [order(coord.L93 [,1],decreasing =F),] 
coord.L93 <- as.data.frame(coord.L93)
coord.L93 <- rename.vars(coord.L93, "id.profil", "id_profil")

### Transform to spatial
coord.L93.sp <- coord.L93
coordinates(coord.L93.sp) <- ~ x + y
proj4string(coord.L93.sp) <- CRS("+init=epsg:2154")
rm(coord.wgs84)

### Assign the department to which it belongs
### Load shapefile with departments
departments<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "DEPARTEMENT")
departments <- spTransform(departments, CRS("+init=epsg:2154"))

plot(departments, col="darkcyan")  
points(coord.L93.sp, pch=1, cex=0.5)

###Extract the department
departments2 <- departments[,c("CODE_DEPT", "CODE_REG", "NOM_DEPT", "NOM_REGION")]
dept_profiles <- over(coord.L93.sp,departments2)
coord.L93.sp@data <- cbind(coord.L93.sp@data, dept_profiles)

# delete soil profils in corse
dept <- c("2A", "2B")
coord.L93.sp <- coord.L93.sp[coord.L93.sp@data$CODE_DEPT %nin% dept, ]
unique(coord.L93.sp@data$CODE_DEPT)
### those that are in the middle of the sea
coord.L93.sp <- coord.L93.sp[!(is.na(coord.L93.sp@data$CODE_DEPT)), ]

### Attach the communes
communes<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "COMMUNE")
communes <- spTransform(communes, CRS("+init=epsg:2154"))
plot(communes, col="pink")  
points(coord.L93.sp, pch=1, cex=0.5)

communes2 <- communes[, c("INSEE_COM")]
comm_profiles <- over(coord.L93.sp,communes2)
coord.L93.sp@data <- cbind(coord.L93.sp@data, comm_profiles)
rm(comm_profiles,dept_profiles,dept, id.profil, ds3)

# Delete soil profils from other islands
iles <- c(56009, 56114, 56152, 56241, 29155, 56069, 29083, 85113, 85163, 17019, 17121, 17207, 17286, 17318, 17051, 17161, 17297, 17360, 17369, 17093, 17140, 17485, 17411, 17486, 17323, 17337, 17385, 83069)
iles <- as.factor(iles)
coord.L93.sp <- coord.L93.sp[coord.L93.sp@data$INSEE_COM %nin% iles,]

# supression des profils se trouvant sur les îles de la commune de Hyères
prof.iles <- c(40052, 40054)
prof.iles <- as.factor(prof.iles)
coord.L93.sp <- coord.L93.sp[coord.L93.sp@data$INSEE_COM %nin% prof.iles,]

coord.L93.sp@data <- drop.levels(coord.L93.sp@data)

# supression des profils hors France (no INSEE code)
coord.L93.sp <- coord.L93.sp[!(is.na(coord.L93.sp@data$INSEE_COM)), ]

# supression des profils sans coordonnées
# ???

### Merge to soil data
coord.L93 <- as.data.frame(coord.L93.sp)

granulo.hor.till.sp <- merge(granulo.hor.till, coord.L93, by.x="id_profil", by.y="id_profil", all.x=F, all.y=F)

# derniers problèmes de  profondeur d'horizon
a <- subset(granulo.hor.till.sp, prof_sup_moy >= prof_inf_moy)
dim(a) # 
rm(a)
#granulo.hor <- subset(granulo.hor, prof_sup_moy < prof_inf_moy)
rm(iles, prof.iles, departments, communes)
save(list = ls(), file = "preparations_donnees_granulo_2.RData")
odbcCloseAll()

##########################################################################################################################################
## Récupération des autres données sol
##########################################################################################################################################

# Récupération des données sur les éléments grossiers estimés (% massiques) pour les horizons
### ??? aren't these volumique?

### for my study on AWC, I will use abondance_eg as (estimated) coarse elements (% volume).

eg.est <- extract.data("eg.est",login,password)
eg.est <- orderBy(~id_profil+no_horizon, eg.est)
eg.est <- add.id(eg.est, c("id_profil","no_horizon"), "id_tot1")
#eg.est <- add.id(eg.est, c("id_profil","prof_sommet","prof_base"), "id_tot")
eg.est <- subset(eg.est, id_tot1 %in% granulo.hor.till.sp$id_tot1)

# attention : quand il n'y a pas d'eg secondaire, le pourcentage d'eg peut-être renseigné dans le champ abon_eg uniquement, ou dans le champ abond_eg_prin en plus
# modification de la table pour que le champ  abond_eg_prin soit toujours renseigné.
eg.est$abondance_eg_prin[!is.na(eg.est$abondance_eg) & is.na(eg.est$abondance_eg_prin) &  is.na(eg.est$abondance_eg_sec)] <- eg.est$abondance_eg[!is.na(eg.est$abondance_eg) & is.na(eg.est$abondance_eg_prin) &  is.na(eg.est$abondance_eg_sec)]
eg.est$abondance_eg_prin[!is.na(eg.est$abondance_eg) & is.na(eg.est$abondance_eg_prin) &  !is.na(eg.est$abondance_eg_sec)] <- eg.est$abondance_eg[!is.na(eg.est$abondance_eg) & is.na(eg.est$abondance_eg_prin) &  !is.na(eg.est$abondance_eg_sec)] - eg.est$abondance_eg_sec[!is.na(eg.est$abondance_eg) & is.na(eg.est$abondance_eg_prin) &  !is.na(eg.est$abondance_eg_sec)]
# remplacement des NA du champ abon_eg_sec par des 0
eg.est$abondance_eg_sec[is.na(eg.est$abondance_eg_sec)] <- 0

eg.est <- subset(eg.est, !(abondance_eg<0 | abondance_eg_prin<0 | abondance_eg_sec<0))

# vérification de la concordance entre les champs abondance_eg, abondance_eg_prin  et abondance_eg_sec
# eg.est <- subset(eg.est, !(((abondance_eg_prin + abondance_eg_sec) - abondance_eg) > 0))


# eg estimées par horizon
data <- merge(granulo.hor.till.sp, eg.est[, c("id_profil", "no_horizon", "abondance_eg", "abondance_eg_prin", "abondance_eg_sec", "id_tot1")], by = "id_tot1", all.x = T, all.y = F, suffixes = c("",".x"))
data  <- subset(data, select = -c(id_profil.x, no_horizon.x))
# rm(eg.est); gc()


## Save data
save(list = ls(), file = "preparation_data_granulo_3.RData")


# Carte des profils sélectionnés
##########################################################################################################################################

Fr.ggmap <- fortify(departments2, region="CODE_DEPT") # pour le mettre au format data.frame
maptot <- ggplot(Fr.ggmap, aes(long, lat))+ geom_polygon(aes(group=group), fill=NA, colour="black", size=.3) +geom_point(data=data, aes(x, y), size=0.5) + theme_bw()
save(maptot, file = "map_tot.RData")
pdf("maptot.pdf")
maptot
dev.off()
rm(coord.L93, eg.est,Fr.ggmap,communes2,departments2, coord.L93.sp,maptot)

########################################################################################################################################

### Add data on max soil depth

# Récupération des données sur la profondeur du sol
prof.sol <- extract.data("prof.sol",login,password)
write.table(subset(prof.sol, prof_sol_p <= 0), "probleme_prof_sol.csv", sep=";")
prof.sol <- subset(prof.sol, id_profil %in% data$id_profil & prof_sol_p > 0 & !is.na(prof_sol_p))
length(unique(prof.sol$id_profil)) ### Only 8674 profiles have soil depth??!!!!
length(unique(data$id_profil)) ### 40225

# Calcul de la profondeur moyenne des horizons
depth <- rep(NA, nrow(data))
depth <- data$prof_sup_moy + ((data$prof_inf_moy - data$prof_sup_moy)/2)
hist(depth);summary(depth)
depth <- cbind(as.character(data$id_tot), depth); colnames(depth) <- c("id_tot","depth")
depth <- as.data.frame(depth)
depth$depth <- as.numeric(as.character(depth$depth));hist(depth$depth);summary(depth$depth)
depth <- subset(depth, id_tot %in% data$id_tot)

# depth
data <- merge(data, depth, by = "id_tot", all.x = T, all.y = F, suffixes = c("",".x"))
rm(depth);gc()

# profondeur du sol, raison de l'arrêt du sondage, abondance et distribution des racines dans le profil
data <- merge(data, prof.sol,  by = "id_profil", all.x = T, all.y = F, suffixes = c("",".x"))
max.depth.hz <- summaryBy(prof_inf_moy~id_profil, data,FUN=max)  ###For each profile, I extract the maximum lower depth form all its horizons
# max.depth.hz$prof_inf_moy.max <- max.depth.hz$prof_inf_moy.max + 30 ## We add 30 cm to be safer (???)
data <- merge(data, max.depth.hz, by="id_profil", all.x=T, all.y=F)
data <- rename.vars(data, "prof_inf_moy.max", "prof.sol.est")
data$prof.sol.est[!is.na(data$prof_sol_p)] <- data$prof_sol_p[!is.na(data$prof_sol_p)] ## When there is soil depth in the profile description, we assign that one instead. So.... for just 8777 profiles

data.awc <- data

## Clean 
rm(data,eg.est,max.depth.hz, prof.sol)

save(list = ls(), file = "preparation_data_granulo_4.RData")
## end of the script