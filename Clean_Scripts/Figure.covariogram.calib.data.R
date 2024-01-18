###################################################################################################################################


####### Load packages
library(sp)
library(rgeos)
library(gstat)
library(rgdal)
library(raster)
library(lattice)
library(ggplot2)
library(plyr)
library(ranger)
library(quantregForest)
library(doParallel)
library(foreach)
library(GGally)


##### Figures for article
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/Figures")


### data on granulo
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/4.Cubist/4.3-Cubist_15_30.RData")
save(granulo.data_15_30, file="granulo.data_15_30.RData")
#rm(list=ls())

### data on coarse elements
### Load the calibration datasets
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/3.3-CalDatChelsa/CalibrationDataCoarse.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/8.-CoarseElementsMaps/qrF.models.CoarseElements.RData")
### Predict on igcs locations (to calculate model residuals)
coarse.data_15_30 <- df.output[[3]]
coarse.data_15_30[,c(paste0("Unt.","coarse_15_30"))] <- exp(coarse.data_15_30[,"coarse_15_30"])
coarse.data_15_30 <- coarse.data_15_30[ ,c("Unt.coarse_15_30","x", "y","id_profil","bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                                           "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                                           "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                           "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                                           "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")]
### select only complete cases and drop empty factor levels
coarse.data_15_30 <- coarse.data_15_30[complete.cases(coarse.data_15_30),]
coarse.data_15_30[] <- lapply(coarse.data_15_30, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

### Predict mean
coarse.data_15_30$pred.qrf.mean.coarse_15_30 <- predict(object =  Rocks.qrF[[3]],newdata=coarse.data_15_30, what=mean)
coarse.data_15_30$resid.qrf.coarse_15_30 <- coarse.data_15_30$Unt.coarse_15_30 - coarse.data_15_30$pred.qrf.mean.coarse_15_30

save(coarse.data_15_30, file="coarse.data_15_30.RData")
#rm(list=ls())

## Get only what I need
load("granulo.data_15_30.RData")
load("coarse.data_15_30.RData")

##############################################################################################################################

### Exploratory analysis

## Transform to spatial our texture data
granulo.data_15_30.sp <- granulo.data_15_30
coordinates(granulo.data_15_30.sp) <- ~ x + y
proj4string(granulo.data_15_30.sp) <- CRS("+init=epsg:2154")

coarse.data_15_30.sp <- coarse.data_15_30
coordinates(coarse.data_15_30.sp) <- ~ x + y
proj4string(coarse.data_15_30.sp) <- CRS("+init=epsg:2154")

g_15_30 <- gstat(NULL,id="Silt.alr", formula = alr.limon_15_30 ~ 1, data = granulo.data_15_30.sp[!is.na(granulo.data_15_30.sp@data$alr.limon_15_30),])
g_15_30 <- gstat(g_15_30,id="Clay.alr", formula = alr.argile_15_30 ~ 1, data = granulo.data_15_30.sp[!is.na(granulo.data_15_30.sp@data$alr.argile_15_30),])
g_15_30 <- gstat(g_15_30,id="Coarse elements", formula = Unt.coarse_15_30 ~ 1, data = coarse.data_15_30.sp[!is.na(coarse.data_15_30.sp@data$Unt.coarse_15_30),])

#v.cross <- variogram(g_0_5, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
v.cross <- variogram(g_15_30,width = 12000)
str(v.cross)
plot(v.cross, pl=FALSE)

# ### rescale residuals
# granulo.data_15_30.sp@data$resid.cub.alr.limon_15_30.res <- granulo.data_15_30.sp@data$resid.cub.alr.limon_15_30/sd(granulo.data_15_30.sp@data$resid.cub.alr.limon_15_30)
# granulo.data_15_30.sp@data$resid.cub.alr.argile_15_30.res <- granulo.data_15_30.sp@data$resid.cub.alr.argile_15_30/sd(granulo.data_15_30.sp@data$resid.cub.alr.argile_15_30)
# coarse.data_15_30.sp@data$resid.qrf.coarse_15_30.res <- coarse.data_15_30.sp@data$resid.qrf.coarse_15_30/sd(coarse.data_15_30.sp@data$resid.qrf.coarse_15_30)
# 
# r.g_15_30 <- gstat(NULL,id="Silt.alr", formula = resid.cub.alr.limon_15_30.res ~ 1, data = granulo.data_15_30.sp[!is.na(granulo.data_15_30.sp@data$resid.cub.alr.limon_15_30),])
# r.g_15_30 <- gstat(r.g_15_30,id="Clay.alr", formula = resid.cub.alr.argile_15_30.res ~ 1, data = granulo.data_15_30.sp[!is.na(granulo.data_15_30.sp@data$resid.cub.alr.argile_15_30),])
# r.g_15_30 <- gstat(r.g_15_30,id="Coarse elements", formula = resid.qrf.coarse_15_30.res ~ 1, data = coarse.data_15_30.sp[!is.na(coarse.data_15_30.sp@data$resid.qrf.coarse_15_30),])
# 
# #v.cross <- variogram(g_0_5, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
# r.v.cross <- variogram(r.g_15_30,width = 16000)
# str(r.v.cross)
# plot(r.v.cross, pl=FALSE)
# 
# ###Fit vgm model
# clay_alr.g <- gstat(formula = resid.cub.alr.argile_15_30.res ~ 1, data = granulo.data_15_30.sp[!is.na(granulo.data_15_30.sp@data$resid.cub.alr.argile_15_30),])
# clay_alr.svar <- variogram(clay_alr.g, width = 10000, cutoff=300000)
# plot(clay_alr.svar, plot.nu=FALSE)
# clay_alr.Sph.vgm <- vgm(nugget=0.5, range=200000, psill=0.5,  model="Sph"); plot(clay_alr.svar, clay_alr.Sph.vgm)
# clay_alr.fit.vgm <- fit.variogram(clay_alr.svar, model=clay_alr.Sph.vgm, fit.method=2)
# plot(clay_alr.svar, clay_alr.fit.vgm, pl=FALSE,  col="blue", main="Sph vgm clay-alr (15-30 cm)")
# SSVR(clay_alr.fit.vgm$psill[[1]], clay_alr.fit.vgm$psill[[2]])  ## 52.4 %
# NtoS(clay_alr.fit.vgm$psill[[1]], clay_alr.fit.vgm$psill[[2]])  ## 0.48
# 
# ### use these parameters for LMCR
# r.g_15_30 <- gstat(r.g_15_30, id = "Clay.alr", model = clay_alr.fit.vgm, fill.all=T)
# r.g_15_30 <- fit.lmc(r.v.cross, r.g_15_30, correct.diagonal=1.01)
# 
# ## Save the plot
# #jpeg(filename="ALRSiltClay_LMCR_m7.jpeg", width=700, height=600)
# plot(r.v.cross, model=r.g_15_30$model)
# #dev.off()

######################################################################################################################

### Empirical variogram, witout rescaling
r.g_15_30 <- gstat(NULL,id="Silt.alr", formula = resid.cub.alr.limon_15_30 ~ 1, data = granulo.data_15_30.sp[!is.na(granulo.data_15_30.sp@data$resid.cub.alr.limon_15_30),])
r.g_15_30 <- gstat(r.g_15_30,id="Clay.alr", formula = resid.cub.alr.argile_15_30 ~ 1, data = granulo.data_15_30.sp[!is.na(granulo.data_15_30.sp@data$resid.cub.alr.argile_15_30),])
r.g_15_30 <- gstat(r.g_15_30,id="Coarse elements", formula = resid.qrf.coarse_15_30 ~ 1, data = coarse.data_15_30.sp[!is.na(coarse.data_15_30.sp@data$resid.qrf.coarse_15_30),])

#v.cross <- variogram(g_0_5, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
r.v.cross <- variogram(r.g_15_30, width = 12000)
str(r.v.cross)
plot(r.v.cross, pl=FALSE)

### save to file
stripParams <- list(cex=3, lines=1)
jpeg("Figure1.c.jpeg", width=1772, height= 1575, units="px")
par(mar=c(2,2,0,0), oma=c(4,4,4,4), las=1)
plot(r.v.cross, pl=FALSE, cex=4, par.strip.text = stripParams,cex.axis=2,
     pch=20, col="black",xlab=list(cex=3, font=2), ylab=list(cex=3, font=2) )
dev.off()




#######################################################################################################################

### Scatterplots
data.15_30 <- merge(coarse.data_15_30, granulo.data_15_30, by="id_profil", all=TRUE)
data.15_30 <- data.15_30[,c("resid.cub.alr.argile_15_30","resid.cub.alr.limon_15_30","resid.qrf.coarse_15_30")]
colnames(data.15_30) <- c("Clay.alr residuals", "Silt.alr residuals", "Coarse elements residuals")
pairs(data.15_30, pch=1)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(GGally)

rcorr(as.matrix(data.15_30))
pairs(data.15_30)
M <- cor(data.15_30, use = "na.or.complete")
p <- cor.mtest(data.15_30, use = "na.or.complete")
par(mfrow=c(1,1))


### Option 1: corrplot
corrplot.mixed(M, upper = "ellipse", lower="number", p.mat=p$p, insig="blank", diag = "u")

### Option 2: Base code
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    ### Correlation
    r <- cor(x, y, use = "na.or.complete")
    # p-value calculation
    p <- cor.test(x, y)$p.value
    
    txt <- format(c(r, 0.2), digits = digits, nsmall=2)[1]
    txt <- paste0(prefix, txt)
    
    if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
    
    if(p>0.05) txt <- paste0("")
    else if(p<=0.05)   text(0.5, 0.5, txt, cex = cex.cor)
    
}

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks=10)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "gray70", ...)
}

par(mar=c(40,120,80,120), oma=c(30,120,80,120), las=0, cex=4)
jpeg("Figure2.jpeg", width=1772, height= 1772, units="px")
pairs(data.15_30, upper.panel = panel.cor, diag.panel =panel.hist,
      pch=19, gap=0.5, cex.labels = 5, font.labels=2,las=0,
      cex.axis = 3, mgp=c(0,2,0))
dev.off()

par()



### Option 3: GGally 
# take the inner 5 colors of 7 colors so that the extreme color is not so intense
display.brewer.pal(7, "RdBu")
corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdBu")[2:6]
corColors
#[1] "#FC8D59" "#FEE090" "#FFFFBF" "#E0F3F8" "#91BFDB"

my_custom_cor_color <- function(data, mapping, color = I("black"), sizeRange = c(1, 10), ...) {
    
    # get the x and y data to use the other code
    x <- eval(mapping$x, data)
    y <- eval(mapping$y, data)
    
    ct <- cor.test(x,y)
    
    sig <- symnum(
        ct$p.value, corr = FALSE, na = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("***", "**", "*", ".", " ")
    )
    
    
    r <- unname(ct$estimate)
    rt <- format(r, digits=2)[1]
    tt <- as.character(rt)
    
    # plot the cor value
    p <- ggally_text(
        label = tt, 
        mapping = aes(),
        xP = 0.5, yP = 0.5, 
        size = 20,
        color=color,
        ...
    ) +
        
        # add the sig stars
        geom_text(
            aes_string(
                x = 0.7,
                y = 0.7
            ),
            label = sig, 
            size = 10,
            color = color,
            ...
        ) +
        
        
        
        theme(
            panel.background = element_rect(fill="white", color = "black", linetype = "dashed"),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank()
        ) 
    
    corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdBu")[2:6]
    
    if (r <= -0.8) {
        corCol <- corColors[1]
    } else if (r <= -0.6) {
        corCol <- corColors[2]
    } else if (r < 0.6) {
        corCol <- corColors[3]
    } else if (r < 0.8) {
        corCol <- corColors[4]
    } else {
        corCol <- corColors[5]
    }
    p <- p + ggplot2::theme_bw() + theme(
        panel.background = element_rect(fill= corCol)
    )
    
    p
}


ggpairs(data.15_30, lower=list(continuous="points"),upper = list(continuous = my_custom_cor_color)) + theme_bw()

jpeg("Figure2.jpeg", width=1772, height= 1772, units="px")
par(mar=c(4,4,0,4), oma=c(4,4,0,0), las=1)
ggpairs(data.15_30, lower=list(continuous="points"),upper = list(continuous = my_custom_cor_color)) + theme_bw()
dev.off()

##########################################################################################################################

### check for all other depths

##########################################################################################################################

setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/cross.variograms")
#  ### 0-5 cm -------------------------------------------------------------

### data on granulo
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/4.Cubist/4.1-Cubist_0_5.RData")
save(granulo.data_0_5, file="granulo.data_0_5.RData")
#rm(list=ls())

### data on coarse elements
### Load the calibration datasets
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/3.3-CalDatChelsa/CalibrationDataCoarse.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/8.-CoarseElementsMaps/qrF.models.CoarseElements.RData")
### Predict on igcs locations (to calculate model residuals)
coarse.data_0_5 <- df.output[[1]]
coarse.data_0_5[,c(paste0("Unt.","coarse_0_5"))] <- exp(coarse.data_0_5[,"coarse_0_5"])
coarse.data_0_5 <- coarse.data_0_5[ ,c("Unt.coarse_0_5","x", "y","id_profil","bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                                           "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                                           "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                           "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                                           "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")]
### select only complete cases and drop empty factor levels
coarse.data_0_5 <- coarse.data_0_5[complete.cases(coarse.data_0_5),]
coarse.data_0_5[] <- lapply(coarse.data_0_5, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

### Predict mean
coarse.data_0_5$pred.qrf.mean.coarse_0_5 <- predict(object =  Rocks.qrF[[3]],newdata=coarse.data_0_5, what=mean)
coarse.data_0_5$resid.qrf.coarse_0_5 <- coarse.data_0_5$Unt.coarse_0_5 - coarse.data_0_5$pred.qrf.mean.coarse_0_5

save(coarse.data_0_5, file="coarse.data_0_5.RData")
#rm(list=ls())

## Get only what I need
load("granulo.data_0_5.RData")
load("coarse.data_0_5.RData")

## Transform to spatial our texture data
granulo.data_0_5.sp <- granulo.data_0_5
coordinates(granulo.data_0_5.sp) <- ~ x + y
proj4string(granulo.data_0_5.sp) <- CRS("+init=epsg:2154")

coarse.data_0_5.sp <- coarse.data_0_5
coordinates(coarse.data_0_5.sp) <- ~ x + y
proj4string(coarse.data_0_5.sp) <- CRS("+init=epsg:2154")

### Empirical variogram, witout rescaling
r.g_0_5 <- gstat(NULL,id="Silt.alr", formula = resid.cub.alr.limon_0_5 ~ 1, data = granulo.data_0_5.sp[!is.na(granulo.data_0_5.sp@data$resid.cub.alr.limon_0_5),])
r.g_0_5 <- gstat(r.g_0_5,id="Clay.alr", formula = resid.cub.alr.argile_0_5 ~ 1, data = granulo.data_0_5.sp[!is.na(granulo.data_0_5.sp@data$resid.cub.alr.argile_0_5),])
r.g_0_5 <- gstat(r.g_0_5,id="Coarse elements", formula = resid.qrf.coarse_0_5 ~ 1, data = coarse.data_0_5.sp[!is.na(coarse.data_0_5.sp@data$resid.qrf.coarse_0_5),])

#v.cross <- variogram(g_0_5, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
r.v.cross <- variogram(r.g_0_5, width = 20000)
str(r.v.cross)
plot(r.v.cross, pl=FALSE)

jpeg("cross.vgm_0_5.jpeg", width=500, height= 500, units="px")
plot(r.v.cross, pl=FALSE)
dev.off()


data.0_5 <- merge(coarse.data_0_5, granulo.data_0_5, by="id_profil", all=TRUE)
data.0_5 <- data.0_5[,c("resid.cub.alr.argile_0_5","resid.cub.alr.limon_0_5","resid.qrf.coarse_0_5")]
jpeg("Resid_0_5.jpeg", width=1000, height= 1000, units="px")
ggpairs(data.0_5, lower=list(continuous="points"),upper = list(continuous = my_custom_cor_color)) + theme_bw()
dev.off()


# ### 5-15 cm -------------------------------------------------------------

### data on granulo
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/4.Cubist/4.2-Cubist_5_15.RData")
save(granulo.data_5_15, file="granulo.data_5_15.RData")
#rm(list=ls())

### data on coarse elements
### Load the calibration datasets
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/3.3-CalDatChelsa/CalibrationDataCoarse.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/8.-CoarseElementsMaps/qrF.models.CoarseElements.RData")
### Predict on igcs locations (to calculate model residuals)
coarse.data_5_15 <- df.output[[2]]
coarse.data_5_15[,c(paste0("Unt.","coarse_5_15"))] <- exp(coarse.data_5_15[,"coarse_5_15"])
coarse.data_5_15 <- coarse.data_5_15[ ,c("Unt.coarse_5_15","x", "y","id_profil","bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                                       "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                                       "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                       "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                                       "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")]
### select only complete cases and drop empty factor levels
coarse.data_5_15 <- coarse.data_5_15[complete.cases(coarse.data_5_15),]
coarse.data_5_15[] <- lapply(coarse.data_5_15, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

### Predict mean
coarse.data_5_15$pred.qrf.mean.coarse_5_15 <- predict(object =  Rocks.qrF[[3]],newdata=coarse.data_5_15, what=mean)
coarse.data_5_15$resid.qrf.coarse_5_15 <- coarse.data_5_15$Unt.coarse_5_15 - coarse.data_5_15$pred.qrf.mean.coarse_5_15

save(coarse.data_5_15, file="coarse.data_5_15.RData")
#rm(list=ls())

## Get only what I need
load("granulo.data_5_15.RData")
load("coarse.data_5_15.RData")

## Transform to spatial our texture data
granulo.data_5_15.sp <- granulo.data_5_15
coordinates(granulo.data_5_15.sp) <- ~ x + y
proj4string(granulo.data_5_15.sp) <- CRS("+init=epsg:2154")

coarse.data_5_15.sp <- coarse.data_5_15
coordinates(coarse.data_5_15.sp) <- ~ x + y
proj4string(coarse.data_5_15.sp) <- CRS("+init=epsg:2154")

### Empirical variogram, witout rescaling
r.g_5_15 <- gstat(NULL,id="Silt.alr", formula = resid.cub.alr.limon_5_15 ~ 1, data = granulo.data_5_15.sp[!is.na(granulo.data_5_15.sp@data$resid.cub.alr.limon_5_15),])
r.g_5_15 <- gstat(r.g_5_15,id="Clay.alr", formula = resid.cub.alr.argile_5_15 ~ 1, data = granulo.data_5_15.sp[!is.na(granulo.data_5_15.sp@data$resid.cub.alr.argile_5_15),])
r.g_5_15 <- gstat(r.g_5_15,id="Coarse elements", formula = resid.qrf.coarse_5_15 ~ 1, data = coarse.data_5_15.sp[!is.na(coarse.data_5_15.sp@data$resid.qrf.coarse_5_15),])

#v.cross <- variogram(g_5_15, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
r.v.cross <- variogram(r.g_5_15, width = 20000)
str(r.v.cross)

jpeg("cross.vgm_5_15.jpeg", width=500, height= 500, units="px")
plot(r.v.cross, pl=FALSE)
dev.off()

data.5_15 <- merge(coarse.data_5_15, granulo.data_5_15, by="id_profil", all=TRUE)
data.5_15 <- data.5_15[,c("resid.cub.alr.argile_5_15","resid.cub.alr.limon_5_15","resid.qrf.coarse_5_15")]
jpeg("Resid_5_15.jpeg", width=1000, height= 1000, units="px")
ggpairs(data.5_15, lower=list(continuous="points"),upper = list(continuous = my_custom_cor_color)) + theme_bw()
dev.off()


# ### 30-60 cm -------------------------------------------------------------

### data on granulo
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/4.Cubist/4.4-Cubist_30_60.RData")
save(granulo.data_30_60, file="granulo.data_30_60.RData")
#rm(list=ls())

### data on coarse elements
### Load the calibration datasets
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/3.3-CalDatChelsa/CalibrationDataCoarse.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/8.-CoarseElementsMaps/qrF.models.CoarseElements.RData")
### Predict on igcs locations (to calculate model residuals)
coarse.data_30_60 <- df.output[[4]]
coarse.data_30_60[,c(paste0("Unt.","coarse_30_60"))] <- exp(coarse.data_30_60[,"coarse_30_60"])
coarse.data_30_60 <- coarse.data_30_60[ ,c("Unt.coarse_30_60","x", "y","id_profil","bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                                         "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                                         "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                         "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                                         "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")]
### select only complete cases and drop empty factor levels
coarse.data_30_60 <- coarse.data_30_60[complete.cases(coarse.data_30_60),]
coarse.data_30_60[] <- lapply(coarse.data_30_60, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

### Predict mean
coarse.data_30_60$pred.qrf.mean.coarse_30_60 <- predict(object =  Rocks.qrF[[3]],newdata=coarse.data_30_60, what=mean)
coarse.data_30_60$resid.qrf.coarse_30_60 <- coarse.data_30_60$Unt.coarse_30_60 - coarse.data_30_60$pred.qrf.mean.coarse_30_60

save(coarse.data_30_60, file="coarse.data_30_60.RData")
#rm(list=ls())

## Get only what I need
load("granulo.data_30_60.RData")
load("coarse.data_30_60.RData")

## Transform to spatial our texture data
granulo.data_30_60.sp <- granulo.data_30_60
coordinates(granulo.data_30_60.sp) <- ~ x + y
proj4string(granulo.data_30_60.sp) <- CRS("+init=epsg:2154")

coarse.data_30_60.sp <- coarse.data_30_60
coordinates(coarse.data_30_60.sp) <- ~ x + y
proj4string(coarse.data_30_60.sp) <- CRS("+init=epsg:2154")

### Empirical variogram, witout rescaling
r.g_30_60 <- gstat(NULL,id="Silt.alr", formula = resid.cub.alr.limon_30_60 ~ 1, data = granulo.data_30_60.sp[!is.na(granulo.data_30_60.sp@data$resid.cub.alr.limon_30_60),])
r.g_30_60 <- gstat(r.g_30_60,id="Clay.alr", formula = resid.cub.alr.argile_30_60 ~ 1, data = granulo.data_30_60.sp[!is.na(granulo.data_30_60.sp@data$resid.cub.alr.argile_30_60),])
r.g_30_60 <- gstat(r.g_30_60,id="Coarse elements", formula = resid.qrf.coarse_30_60 ~ 1, data = coarse.data_30_60.sp[!is.na(coarse.data_30_60.sp@data$resid.qrf.coarse_30_60),])

#v.cross <- variogram(g_30_60, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
r.v.cross <- variogram(r.g_30_60, width = 20000)
str(r.v.cross)

jpeg("cross.vgm_30_60.jpeg", width=500, height= 500, units="px")
plot(r.v.cross, pl=FALSE)
dev.off()

data.30_60 <- merge(coarse.data_30_60, granulo.data_30_60, by="id_profil", all=TRUE)
data.30_60 <- data.30_60[,c("resid.cub.alr.argile_30_60","resid.cub.alr.limon_30_60","resid.qrf.coarse_30_60")]
jpeg("Resid_30_60.jpeg", width=1000, height= 1000, units="px")
ggpairs(data.30_60, lower=list(continuous="points"),upper = list(continuous = my_custom_cor_color)) + theme_bw()
dev.off()


# ### 60-100 cm -------------------------------------------------------------

### data on granulo
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/4.Cubist/4.5-Cubist_60_100.RData")
save(granulo.data_60_100, file="granulo.data_60_100.RData")
#rm(list=ls())

### data on coarse elements
### Load the calibration datasets
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/3.3-CalDatChelsa/CalibrationDataCoarse.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/8.-CoarseElementsMaps/qrF.models.CoarseElements.RData")
### Predict on igcs locations (to calculate model residuals)
coarse.data_60_100 <- df.output[[5]]
coarse.data_60_100[,c(paste0("Unt.","coarse_60_100"))] <- exp(coarse.data_60_100[,"coarse_60_100"])
coarse.data_60_100 <- coarse.data_60_100[ ,c("Unt.coarse_60_100","x", "y","id_profil","bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                                           "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                                           "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                           "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                                           "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")]
### select only complete cases and drop empty factor levels
coarse.data_60_100 <- coarse.data_60_100[complete.cases(coarse.data_60_100),]
coarse.data_60_100[] <- lapply(coarse.data_60_100, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

### Predict mean
coarse.data_60_100$pred.qrf.mean.coarse_60_100 <- predict(object =  Rocks.qrF[[3]],newdata=coarse.data_60_100, what=mean)
coarse.data_60_100$resid.qrf.coarse_60_100 <- coarse.data_60_100$Unt.coarse_60_100 - coarse.data_60_100$pred.qrf.mean.coarse_60_100

save(coarse.data_60_100, file="coarse.data_60_100.RData")
#rm(list=ls())

## Get only what I need
load("granulo.data_60_100.RData")
load("coarse.data_60_100.RData")

## Transform to spatial our texture data
granulo.data_60_100.sp <- granulo.data_60_100
coordinates(granulo.data_60_100.sp) <- ~ x + y
proj4string(granulo.data_60_100.sp) <- CRS("+init=epsg:2154")

coarse.data_60_100.sp <- coarse.data_60_100
coordinates(coarse.data_60_100.sp) <- ~ x + y
proj4string(coarse.data_60_100.sp) <- CRS("+init=epsg:2154")

### Empirical variogram, witout rescaling
r.g_60_100 <- gstat(NULL,id="Silt.alr", formula = resid.cub.alr.limon_60_100 ~ 1, data = granulo.data_60_100.sp[!is.na(granulo.data_60_100.sp@data$resid.cub.alr.limon_60_100),])
r.g_60_100 <- gstat(r.g_60_100,id="Clay.alr", formula = resid.cub.alr.argile_60_100 ~ 1, data = granulo.data_60_100.sp[!is.na(granulo.data_60_100.sp@data$resid.cub.alr.argile_60_100),])
r.g_60_100 <- gstat(r.g_60_100,id="Coarse elements", formula = resid.qrf.coarse_60_100 ~ 1, data = coarse.data_60_100.sp[!is.na(coarse.data_60_100.sp@data$resid.qrf.coarse_60_100),])

#v.cross <- variogram(g_60_100, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
r.v.cross <- variogram(r.g_60_100, width = 20000)
str(r.v.cross)

jpeg("cross.vgm_60_100.jpeg", width=500, height= 500, units="px")
plot(r.v.cross, pl=FALSE)
dev.off()

data.60_100 <- merge(coarse.data_60_100, granulo.data_60_100, by="id_profil", all=TRUE)
data.60_100 <- data.60_100[,c("resid.cub.alr.argile_60_100","resid.cub.alr.limon_60_100","resid.qrf.coarse_60_100")]
jpeg("Resid_60_100.jpeg", width=1000, height= 1000, units="px")
ggpairs(data.60_100, lower=list(continuous="points"),upper = list(continuous = my_custom_cor_color)) + theme_bw()
dev.off()

# ### 60-100 cm -------------------------------------------------------------

### data on granulo
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/4.Cubist/4.6-Cubist_100_200.RData")
save(granulo.data_100_200, file="granulo.data_100_200.RData")
#rm(list=ls())

### data on coarse elements
### Load the calibration datasets
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/3.3-CalDatChelsa/CalibrationDataCoarse.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Coarse_elements/Output/RData/8.-CoarseElementsMaps/qrF.models.CoarseElements.RData")
### Predict on igcs locations (to calculate model residuals)
coarse.data_100_200 <- df.output[[6]]
coarse.data_100_200[,c(paste0("Unt.","coarse_100_200"))] <- exp(coarse.data_100_200[,"coarse_100_200"])
coarse.data_100_200 <- coarse.data_100_200[ ,c("Unt.coarse_100_200","x", "y","id_profil","bdforet","cti","curv_long","curv_trans","curvature","ecoclim","eros","etpMax",
                                             "etpMean", "etpMedian","etpMin","EVI_median_june","exposition","graviAL2000","graviBg2000","graviG2000",
                                             "hli","idpr","linear_aspect","modelclc","mrrtf","mrvbf","NDVI_median_jan","NDVI_median_june","roughness",
                                             "rrMax","rrMean","rrMedian","rrMin","sar","scale_pos","slope","slopeascos","slopeassin","slopeastrasp",
                                             "soil1","srr","srtm","tMeanMax","tMeanMean","tMeanMedian","tMeanMin","MAT_RF","AGL_RF")]
### select only complete cases and drop empty factor levels
coarse.data_100_200 <- coarse.data_100_200[complete.cases(coarse.data_100_200),]
coarse.data_100_200[] <- lapply(coarse.data_100_200, function(x) if (is.factor(x)) x[, drop=TRUE] else x)

### Predict mean
coarse.data_100_200$pred.qrf.mean.coarse_100_200 <- predict(object =  Rocks.qrF[[3]],newdata=coarse.data_100_200, what=mean)
coarse.data_100_200$resid.qrf.coarse_100_200 <- coarse.data_100_200$Unt.coarse_100_200 - coarse.data_100_200$pred.qrf.mean.coarse_100_200

save(coarse.data_100_200, file="coarse.data_100_200.RData")
#rm(list=ls())

## Get only what I need
load("granulo.data_100_200.RData")
load("coarse.data_100_200.RData")

## Transform to spatial our texture data
granulo.data_100_200.sp <- granulo.data_100_200
coordinates(granulo.data_100_200.sp) <- ~ x + y
proj4string(granulo.data_100_200.sp) <- CRS("+init=epsg:2154")

coarse.data_100_200.sp <- coarse.data_100_200
coordinates(coarse.data_100_200.sp) <- ~ x + y
proj4string(coarse.data_100_200.sp) <- CRS("+init=epsg:2154")

### Empirical variogram, witout rescaling
r.g_100_200 <- gstat(NULL,id="Silt.alr", formula = resid.cub.alr.limon_100_200 ~ 1, data = granulo.data_100_200.sp[!is.na(granulo.data_100_200.sp@data$resid.cub.alr.limon_100_200),])
r.g_100_200 <- gstat(r.g_100_200,id="Clay.alr", formula = resid.cub.alr.argile_100_200 ~ 1, data = granulo.data_100_200.sp[!is.na(granulo.data_100_200.sp@data$resid.cub.alr.argile_100_200),])
r.g_100_200 <- gstat(r.g_100_200,id="Coarse elements", formula = resid.qrf.coarse_100_200 ~ 1, data = coarse.data_100_200.sp[!is.na(coarse.data_100_200.sp@data$resid.qrf.coarse_100_200),])

#v.cross <- variogram(g_100_200, boundaries = c(1000,1000,2000,2000,2000,2000,3000,3000,3000,5000,5000,rep(6000,240)))
r.v.cross <- variogram(r.g_100_200, width = 20000)
str(r.v.cross)

jpeg("cross.vgm_100_200.jpeg", width=500, height= 500, units="px")
plot(r.v.cross, pl=FALSE)
dev.off()

data.100_200 <- merge(coarse.data_100_200, granulo.data_100_200, by="id_profil", all=TRUE)
data.100_200 <- data.100_200[,c("resid.cub.alr.argile_100_200","resid.cub.alr.limon_100_200","resid.qrf.coarse_100_200")]
jpeg("Resid_100_200.jpeg", width=1000, height= 1000, units="px")
ggpairs(data.100_200, lower=list(continuous="points"),upper = list(continuous = my_custom_cor_color)) + theme_bw()
dev.off()

### end of the script
