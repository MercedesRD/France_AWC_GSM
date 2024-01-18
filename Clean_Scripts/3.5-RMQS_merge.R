#################################################################################################################################################

##################### RMQS FOSSES

######## Validation of particle size and coarse elements predictions

### Prepare the validation data

### date:15/01/2018
## Author: Mercedes Roman

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
library(soiltexture)

### Set working directory
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
dir.create(paste0(HomeDir,"Clean_Output/3.5-RMQS_merge"))
setwd(paste0(HomeDir,"Clean_Output/3.5-RMQS_merge"))

### Load coarse elements data
load(paste0(HomeDir,"Clean_Output/3.3-RMQS_EG/RMQS_EG.RData"))

### order
EG_RMQS <-  orderBy(~ id_profil + no_horizon, EG_RMQS)
### eliminate observations with na(abondance_eg)
EG_RMQS <- EG_RMQS[!is.na(EG_RMQS$abondance_eg),] ### 5363 observations
rmqs.eg.profil <- unique(EG_RMQS$id_profil)
rmqs.eg.profil <- sort(decreasing = FALSE, rmqs.eg.profil)

### Load granulo data
load(paste0(HomeDir,"Clean_Output/3.4-RMQS_granulo/RMQS_granulo.RData"))

### Can we join texture and coarse elements data?
granulo_RMQS[granulo_RMQS$id_profil == rmqs.gr.profil[[6]],]
EG_RMQS[EG_RMQS$id_profil == rmqs.gr.profil[[6]],]

### Now that both tables are ordered by id_profil and no_horizon, I create a unique field for merging

conc.champ <- function(x){
    res <- paste(x,collapse=" ")
}

add.id <- function(table, champs, nom.id){
    new.id <- data.matrix(table[,champs],rownames.force=FALSE)
    new.id <- apply(new.id, 1, conc.champ)
    table$new.id  <- new.id 
    table <- rename.vars(table, "new.id",nom.id, info=FALSE)
    return(table)
}

granulo_RMQS <- add.id(granulo_RMQS, c("id_profil","no_horizon"), "id_tot1")
EG_RMQS <- add.id(EG_RMQS, c("id_profil","no_horizon"), "id_tot1")


################################################################################################################################################

### Clean a bit the EG

### eliminate observations from humus horizons
EG_RMQS <- EG_RMQS[!(EG_RMQS$prof_sup_moy<0),]
EG_RMQS <- EG_RMQS[!(EG_RMQS$prof_inf_moy<0),]

### I eliminate some columns I don't need
EG_RMQS <- EG_RMQS[, !colnames(EG_RMQS) %in% c("prof_sup_min", "prof_sup_max", "prof_inf_min", "prof_inf_max")]
EG_RMQS <- EG_RMQS[, -c(9:22)]

### Eliminate those 86 Nas
EG_RMQS <- EG_RMQS[!is.na(EG_RMQS$id_profil),]


#################################################################################################################################################

### CLean the granulo data

### eliminate observations from humus horizons
granulo_RMQS <- granulo_RMQS[!(granulo_RMQS$profondeur_hz_inf <0),]
granulo_RMQS <- granulo_RMQS[!(granulo_RMQS$profondeur_hz_sup<0),]

### For horizons with more than one observation, calculate average
length((granulo_RMQS$id_tot1)) - length(unique(granulo_RMQS$id_tot1)) ### 263 duplicates

granulo_RMQS <- orderBy(~id_tot1, granulo_RMQS)

### check those prelevements that have more than one observation
prel.dupl <- data.matrix(granulo_RMQS[,"id_tot1"],rownames.force=FALSE)
liste.champ.dupl <- duplicated(prel.dupl)
length(liste.champ.dupl[liste.champ.dupl == TRUE])### we have 263 duplicated id_tot1

### What is their ID?
prel.dupl.id <- granulo_RMQS[liste.champ.dupl,"id_tot1"]
###subset
dupl.prel.df <- granulo_RMQS[granulo_RMQS$id_tot1 %in% prel.dupl.id, ]
granulo_RMQS <- granulo_RMQS[!(granulo_RMQS$id_tot1  %in% prel.dupl.id), ]

### Calculate average clay, silt, and sand
clay <- summaryBy(clay ~ id_tot1, dupl.prel.df, FUN = mean, na.rm=TRUE)
silt <- summaryBy(silt ~ id_tot1, dupl.prel.df, FUN = mean,na.rm=TRUE)
sand <- summaryBy(sand ~ id_tot1, dupl.prel.df,FUN = mean,na.rm=TRUE)

dupl.prel.df2 <- merge(clay, silt, by="id_tot1")
dupl.prel.df2 <- merge(dupl.prel.df2, sand, by="id_tot1")
rm(clay, silt, sand)

### Merge to the previous df
dupl.prel.df <- merge(dupl.prel.df, dupl.prel.df2, by="id_tot1")
### and susbtitutte
dupl.prel.df$clay <- dupl.prel.df$clay.mean
dupl.prel.df$silt <- dupl.prel.df$silt.mean
dupl.prel.df$sand <- dupl.prel.df$sand.mean
dupl.prel.df <- dupl.prel.df[,-c(30:32)]

##3 change order of columns
col.1 <- dupl.prel.df[,1]
dupl.prel.df <- dupl.prel.df[, c(2:29)]
dupl.prel.df$id_tot1 <- col.1; rm(col.1)

names(dupl.prel.df)
names(granulo_RMQS)
### bind back and we are all happy
granulo_RMQS <- rbind(granulo_RMQS, dupl.prel.df)
granulo_RMQS <- orderBy(~id_tot1, granulo_RMQS)

### Delete duplicated rows
del.duplicated1 <- function(table, champs, nom.table){ 
    # Permet de supprimer des lignes dont plusieurs champs sont dupliqués. Par exemple : mêmes limites de prelevement et même valeur mesurée
    champ.dupl <- data.matrix(table[,champs],rownames.force=FALSE)
    champ.dupl <- as.matrix(apply(champ.dupl, 1, conc.champ))
    liste.champ.dupl <- duplicated(champ.dupl)
    if (length(liste.champ.dupl[liste.champ.dupl == TRUE]) == 0) {
        print("pas de doublons")} else{
            prof.dupl <- table[liste.champ.dupl,]
            write.table(prof.dupl ,paste("del_duplicated1_",nom.table,".csv",sep=""),sep=";", col.names = T,row.names = F)
            
            while(length(liste.champ.dupl[liste.champ.dupl == TRUE]) != 0){
                champ.dupl <- data.matrix(table[,champs],rownames.force=FALSE)
                champ.dupl <- as.matrix(apply(champ.dupl, 1, conc.champ))
                liste.champ.dupl <- duplicated(champ.dupl)
                table <- subset(table, !liste.champ.dupl)
            }
        }
    return(table)
}

champs <- c("id_tot1","profondeur_hz_sup","profondeur_hz_inf","clay" ,"silt", "sand")
granulo_RMQS <- del.duplicated1(granulo_RMQS, champs, "granulo")

rm(dupl.prel.df, dupl.prel.df2, liste.champ.dupl, prel.dupl.id, prel.dupl, champs)
granulo_RMQS <- granulo_RMQS[!is.na(granulo_RMQS$clay),]

# ### Merge tables --------------------------------------------------------

### Merge (even if there are some rows missing in the other table)
RMQS.data <- merge(granulo_RMQS, EG_RMQS, by="id_tot1", all.x = TRUE, all.y=TRUE)
head(RMQS.data)
setwd(paste0(HomeDir,"Clean_Output/3.5-RMQS_merge"))
write.csv(granulo_RMQS, "granulo_RMQS.csv")
write.csv(EG_RMQS, "EG_RMQS.csv")
RMQS.data.EG.granulo <- RMQS.data[!is.na(RMQS.data$clay) & !is.na(RMQS.data$abondance_eg),]

rm(add.id, conc.champ, del.duplicated1)

##### data is almost ready for validation
save.image("RMQS.validation.RData")
#load("RMQS.validation.RData")
