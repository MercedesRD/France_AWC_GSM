###################################################################################################################################

################ Project: Mapping AWC for metropolitan France

###  Author: Mercedes Roman Dobarco
###  Date 05/03/2018

### Objective: Improving the predictions of coarse elements

#### Task: Try to create new variable for parent material, based on MAT12 and AGLIM1

library(plyr)

#########################################################
### Set the home directory where all files and subfolders are stored
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 

setwd(paste0(HomeDir,"Input/CoarseElements/BDGSF"))
stu <- read.csv("stu.csv", sep=";")

### Check NEW FIELD, MAT_RF
stu[stu$MAT12 %in% c(0,90),]$MAT_RF <- 0
stu[stu$MAT12 %in% c(11,12,15,31,32,41,42,10,40),]$MAT_RF <- 1
stu[stu$MAT12 %in% c(13,14),]$MAT_RF <- 2
stu[stu$MAT12 %in% c(23,24,34,35,63,30),]$MAT_RF <- 3
stu[stu$MAT12 %in% c(21,25,20),]$MAT_RF <-4
stu[stu$MAT12 %in% c(22),]$MAT_RF <- 5
stu[stu$MAT12 %in% c(43,44,52),]$MAT_RF <- 6
stu[stu$MAT12 %in% c(45,53),]$MAT_RF <- 7
stu[stu$MAT12 %in% c(33,51,50),]$MAT_RF <-8
stu[stu$MAT12 %in% c(61,62,64,60),]$MAT_RF <- 9
stu[stu$MAT12 %in% c(70,71,72,73,74,75),]$MAT_RF <-10
stu[stu$MAT12 %in% c(80,81,82,83),]$MAT_RF <- 11
stu[stu$MAT12 %in% c(91),]$MAT_RF <- 12

### Second classification, based on AGLIM1
stu$AGL_RF <- NA
stu[stu$AGLIM1==0,]$AGL_RF <- 0
stu[stu$AGLIM1 %in% c(1,7:31),]$AGL_RF <- 1
stu[stu$AGLIM1 %in% c(2,5),]$AGL_RF <- 3
stu[stu$AGLIM1==3,]$AGL_RF <- 4
stu[stu$AGLIM1==4,]$AGL_RF <- 5
stu[stu$AGLIM1==6,]$AGL_RF <- 6
#stu[stu$AGLIM1 %in% c(7:31),]$AGL_RF <- 7 I joined others with class "No limitation"

### Read SMU table
smu <- read.csv("smu.csv", sep=";")

### Read LinkSTU and SMU
stu_org <- read.csv("stu_org.csv", sep=";")

stu.RF <- stu[,c("STU","MAT12", "MAT_RF")]

### Add the info to the STU_org table
stu_org.2 <- merge(stu_org, stu.RF, by="STU", all.x=TRUE, all.y=FALSE)
length(unique(stu_org.2$SMU)) ### 318 SMU

### Summarize the results by SMY and MAT_RF
stu_org.3 <- ddply(stu_org.2, .(SMU, MAT_RF), summarize, PC.area=sum(PCAREA))
length(unique(stu_org.3$SMU)) ### 318 SMU

### check that for each SMU there is a MAT_RF with at least 50%
require(data.table)
##3 Two ways of doing it
stu_org.4.2 <- setDT(stu_org.3)[,.SD[which.max(PC.area)],by=SMU]
stu_org.3.1 <- as.data.table(stu_org.3)
stu_org.4 <- stu_org.3.1[stu_org.3.1[, .I[which.max(PC.area)],by=SMU]$V1]
all.equal(stu_org.4, stu_org.4.2) 
rm(stu_org.3.1, stu_org.4.2)
stu_org.4 <- as.data.frame(stu_org.4)

### now I can join this to the SMU table
smu <- merge(smu, stu_org.4, by="SMU")

########################################################################################
stu.AGL <- stu[,c("STU","AGLIM1", "AGL_RF")]

### Add the info to the STU_org table
stu_org.5 <- merge(stu_org, stu.AGL, by="STU", all.x=TRUE, all.y=FALSE)
length(unique(stu_org.5$SMU)) ### 318 SMU

### Summarize the results by SMY and MAT_RF
stu_org.6 <- ddply(stu_org.5, .(SMU, AGL_RF), summarize, PC.area=sum(PCAREA))
length(unique(stu_org.6$SMU)) ### 318 SMU

### check that for each SMU there is a MAT_RF with at least 50%
require(data.table)
##3 Two ways of doing it
stu_org.7.2 <- setDT(stu_org.6)[,.SD[which.max(PC.area)],by=SMU]
stu_org.6.1 <- as.data.table(stu_org.6)
stu_org.7 <- stu_org.6.1[stu_org.6.1[, .I[which.max(PC.area)],by=SMU]$V1]
all.equal(stu_org.7, stu_org.7.2) 
rm(stu_org.7.2, stu_org.6.1)
stu_org.7 <- as.data.frame(stu_org.7)

### now I can join this to the SMU table
smu <- merge(smu, stu_org.7, by="SMU")
colnames(smu) <- c("SMU","NB_POLYS" ,"AREA","NB_STU" ,"MAT_RF" ,   "PCarea1", "AGL_RF" ,   "PCarea2")
### Now I can join to shapefile, and rasterize
write.csv(smu, "smu.RF.csv")

rm(stu_org.2,stu_org.3, stu_org.5, stu_org.6)
############################################################################################################################################
### end of the script