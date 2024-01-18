### figures for article on AWC FRANCE
### author: mercedes Roman Dobarco
### Date

### Supplementary material - Sensitivities and variance terms for AWC elementary variance

library(soiltexture)
library(plyr)
library(knitr)
library(wesanderson)
library(ggplot2)
library(ggtern)
library(GGally)
library(raster)
library(gridExtra)
library(sp)
library(rgdal)
library(viridis)
library(rasterVis)
library(viridisLite)
library(RColorBrewer)
library(scales)

### Load the data
HomeDir <- "D:/romandobarco/AWC/AWC_DSM/"  ### Change this with your Home directory 
load(paste0(HomeDir,"Clean_Output/SensitivityAWC_general/AWC.var.decomp.15_30.RData"))

### Load France (even if we dont use it)
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### create raster stack with all my uncertainty terms
unc.terms.stack <- stack(u.Rocks, u.Texture[[2]], u.Texture[[3]], u.Texture[[4]],
                         u.PTF.FC[[2]], u.PTF.FC[[3]], u.PTF.FC[[4]], u.PTF.FC[[5]],
                         u.PTF.FC[[6]], u.PTF.FC[[7]],
                         u.PTF.PWP[[2]], u.PTF.PWP[[3]], u.PTF.PWP[[4]], u.PTF.PWP[[5]],
                         u.PTF.PWP[[6]], u.PTF.PWP[[7]], var.E.AWC.dec)
names(unc.terms.stack[["var.E.AWC.dec.30_60"]]) <- c("var.AWC.30_60")
#levelplot(var.E.AWC.dec)

## Simply using levelplot
levelplot(unc.terms.stack, par.settings = infernoTheme)
levelplot(unc.terms.stack, par.settings = RdBuTheme) +
    layer(sp.polygons(france, lwd=0.7))

## find some colors
my_colors <- viridis_pal(option="A", direction=-1)(25)
display.brewer.pal(11, "RdYlBu")
my_colors2 <- brewer.pal(11, "RdYlBu")

### Define my intervals
# d <- c(seq(-0.0019, -0.0001, by =0.0005),
#        seq(-0.00008, 0.00008, by =0.00005),
#        seq(0.0001, 0.0028, by =0.0005))

d2 <- c(seq(-0.0019, -0.0001, by =0.0005),-0.0001, 0,
        seq(0.0001, 0.0028, by =0.0005))

### Test figure function
levelplot(unc.terms.stack, 
          margin=FALSE,                       
          colorkey=list(                                                    # suppress marginal graphics 
              space='right',                                               # legend ticks and labels                   
              labels=list(at= d2, font=3),  
              axis.line=list(col='black'),
              width=0.75),  
          par.settings=list(
              strip.border=list(col='transparent'),
              strip.background=list(col='transparent'),
              axis.line=list(col='transparent')),
          scales=list(draw=FALSE),   
          col.regions=my_colors2,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))

# ### also input variables and their variances
# levelplot(awc.input, par.settings = viridisTheme) +
#     layer(sp.polygons(france, lwd=0.7))
# 
# ### finally the sensitivities
# ### Create raster stack
# sensitivities.s <- stack(sens.Rocks, sens.Clay, sens.Silt, sens.Beta0, sens.Beta1, sens.Beta2)
# levelplot(sensitivities.s, par.settings = infernoTheme) +
#     layer(sp.polygons(france, lwd=0.7))
# plot(sensitivities.s)
# plot(awc.input)

### Layer 0-5 cm
### Load the data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.0_5/AWC.var.decomp.0_5.RData")
### set the wd
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.0_5")
### Load France (even if we dont use it)
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### create raster stack with all my uncertainty terms
unc.terms.stack <- stack(u.Rocks, u.Texture[[2]], u.Texture[[3]], u.Texture[[4]],
                         u.PTF.FC[[2]], u.PTF.FC[[3]], u.PTF.FC[[4]], u.PTF.FC[[5]],
                         u.PTF.FC[[6]], u.PTF.FC[[7]],
                         u.PTF.PWP[[2]], u.PTF.PWP[[3]], u.PTF.PWP[[4]], u.PTF.PWP[[5]],
                         u.PTF.PWP[[6]], u.PTF.PWP[[7]], var.E.AWC.dec)
names(unc.terms.stack[["var.E.AWC.dec.0_5"]]) <- c("var.AWC.0_5")

levelplot(unc.terms.stack, par.settings = RdBuTheme) +
    layer(sp.polygons(france, lwd=0.7))
my_colors2 <- brewer.pal(11, "RdYlBu")

d2 <- c(seq(-0.0015, -0.0005, by =0.0005), 
        -0.0001, 0, 0.0001,
        seq(0.0005, 0.0025, by =0.0005),0.0032)

### Test figure function
levelplot(unc.terms.stack, 
          margin=FALSE,                       
          colorkey=list(                                                    # suppress marginal graphics 
              space='right',                                               # legend ticks and labels                   
              labels=list(at= d2, font=3),  
              axis.line=list(col='black'),
              width=0.75),  
          par.settings=list(
              strip.border=list(col='transparent'),
              strip.background=list(col='transparent'),
              axis.line=list(col='transparent')),
          scales=list(draw=FALSE),   
          col.regions=my_colors2,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))

### Layer 5-15 cm
### Load the data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.5_15/AWC.var.decomp.5_15.RData")
### set the wd
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.5_15")
### Load France (even if we dont use it)
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### create raster stack with all my uncertainty terms
unc.terms.stack <- stack(u.Rocks, u.Texture[[2]], u.Texture[[3]], u.Texture[[4]],
                         u.PTF.FC[[2]], u.PTF.FC[[3]], u.PTF.FC[[4]], u.PTF.FC[[5]],
                         u.PTF.FC[[6]], u.PTF.FC[[7]],
                         u.PTF.PWP[[2]], u.PTF.PWP[[3]], u.PTF.PWP[[4]], u.PTF.PWP[[5]],
                         u.PTF.PWP[[6]], u.PTF.PWP[[7]], var.E.AWC.dec)
names(unc.terms.stack[["var.E.AWC.dec.5_15"]]) <- c("var.AWC.5_15")

levelplot(unc.terms.stack, par.settings = RdBuTheme) +
    layer(sp.polygons(france, lwd=0.7))
my_colors2 <- brewer.pal(9, "RdYlBu")

d2 <- c(seq(-0.0015, -0.0005, by =0.0005), 
        -0.0001, 0, 0.0001,
        seq(0.0005, 0.0020, by =0.0005))

### Test figure function
levelplot(unc.terms.stack, 
          margin=FALSE,                       
          colorkey=list(                                                    # suppress marginal graphics 
              space='right',                                               # legend ticks and labels                   
              labels=list(at= d2, font=3),  
              axis.line=list(col='black'),
              width=0.75),  
          par.settings=list(
              strip.border=list(col='transparent'),
              strip.background=list(col='transparent'),
              axis.line=list(col='transparent')),
          scales=list(draw=FALSE),   
          col.regions=my_colors2,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))

### Layer 30-60 cm
### Load the data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.30_60/AWC.var.decomp.30_60.RData")
### set the wd
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.30_60")
### Load France (even if we dont use it)
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### create raster stack with all my uncertainty terms
unc.terms.stack <- stack(u.Rocks, u.Texture[[2]], u.Texture[[3]], u.Texture[[4]],
                         u.PTF.FC[[2]], u.PTF.FC[[3]], u.PTF.FC[[4]], u.PTF.FC[[5]],
                         u.PTF.FC[[6]], u.PTF.FC[[7]],
                         u.PTF.PWP[[2]], u.PTF.PWP[[3]], u.PTF.PWP[[4]], u.PTF.PWP[[5]],
                         u.PTF.PWP[[6]], u.PTF.PWP[[7]], var.E.AWC.dec)
names(unc.terms.stack[["var.E.AWC.dec.30_60"]]) <- c("var.AWC.30_60")

levelplot(unc.terms.stack, par.settings = RdBuTheme) +
    layer(sp.polygons(france, lwd=0.7))
my_colors2 <- brewer.pal(11, "RdYlBu")

d2 <- c(seq(-0.002, -0.0005, by =0.0005), 
        -0.0001, 0, 0.0001,
        seq(0.0005, 0.0020, by =0.0005),0.0030)

### Test figure function
levelplot(unc.terms.stack, 
          margin=FALSE,                       
          colorkey=list(                                                    # suppress marginal graphics 
              space='right',                                               # legend ticks and labels                   
              labels=list(at= d2, font=3),  
              axis.line=list(col='black'),
              width=0.75),  
          par.settings=list(
              strip.border=list(col='transparent'),
              strip.background=list(col='transparent'),
              axis.line=list(col='transparent')),
          scales=list(draw=FALSE),   
          col.regions=my_colors2,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))


### Layer 15-30 cm
### Load the data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.15_30/AWC.var.decomp.15_30.RData")
### set the wd
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.15_30")
### Load France (even if we dont use it)
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### create raster stack with all my uncertainty terms
unc.terms.stack <- stack(u.Rocks, u.Texture[[2]], u.Texture[[3]], u.Texture[[4]],
                         u.PTF.FC[[2]], u.PTF.FC[[3]], u.PTF.FC[[4]], u.PTF.FC[[5]],
                         u.PTF.FC[[6]], u.PTF.FC[[7]],
                         u.PTF.PWP[[2]], u.PTF.PWP[[3]], u.PTF.PWP[[4]], u.PTF.PWP[[5]],
                         u.PTF.PWP[[6]], u.PTF.PWP[[7]], var.E.AWC.dec)
names(unc.terms.stack[["var.E.AWC.dec.15_30"]]) <- c("var.AWC.15_30")

levelplot(unc.terms.stack, par.settings = RdBuTheme) +
    layer(sp.polygons(france, lwd=0.7))
my_colors2 <- brewer.pal(9, "RdYlBu")

d2 <- c(seq(-0.0015, -0.0005, by =0.0005), 
        -0.0001, 0, 0.0001,
        seq(0.0005, 0.0020, by =0.0005))

### Test figure function
levelplot(unc.terms.stack, 
          margin=FALSE,                       
          colorkey=list(                                                    # suppress marginal graphics 
              space='right',                                               # legend ticks and labels                   
              labels=list(at= d2, font=3),  
              axis.line=list(col='black'),
              width=0.75),  
          par.settings=list(
              strip.border=list(col='transparent'),
              strip.background=list(col='transparent'),
              axis.line=list(col='transparent')),
          scales=list(draw=FALSE),   
          col.regions=my_colors2,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))


### Layer 60-100 cm
### Load the data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.60_100/AWC.var.decomp.60_100.RData")
### set the wd
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.60_100")
### Load France (even if we dont use it)
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### create raster stack with all my uncertainty terms
unc.terms.stack <- stack(u.Rocks, u.Texture[[2]], u.Texture[[3]], u.Texture[[4]],
                         u.PTF.FC[[2]], u.PTF.FC[[3]], u.PTF.FC[[4]], u.PTF.FC[[5]],
                         u.PTF.FC[[6]], u.PTF.FC[[7]],
                         u.PTF.PWP[[2]], u.PTF.PWP[[3]], u.PTF.PWP[[4]], u.PTF.PWP[[5]],
                         u.PTF.PWP[[6]], u.PTF.PWP[[7]], var.E.AWC.dec)
names(unc.terms.stack[["var.E.AWC.dec.60_100"]]) <- c("var.AWC.60_100")

levelplot(unc.terms.stack, par.settings = RdBuTheme) +
    layer(sp.polygons(france, lwd=0.7))
my_colors2 <- brewer.pal(11, "RdYlBu")

d2 <- c(seq(-0.0045, -0.0005, by =0.001), 
        -0.0001, 0, 0.0001,
        seq(0.0005, 0.0025, by =0.001),0.004)

### Test figure function
levelplot(unc.terms.stack, 
          margin=FALSE,                       
          colorkey=list(                                                    # suppress marginal graphics 
              space='right',                                               # legend ticks and labels                   
              labels=list(at= d2, font=3),  
              axis.line=list(col='black'),
              width=0.75),  
          par.settings=list(
              strip.border=list(col='transparent'),
              strip.background=list(col='transparent'),
              axis.line=list(col='transparent')),
          scales=list(draw=FALSE),   
          col.regions=my_colors2,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))


# layer 100-200 cm

### Load the data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.100_200/AWC.var.decomp.100_200.RData")
### set the wd
setwd("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.100_200")
### Load France (even if we dont use it)
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### create raster stack with all my uncertainty terms
unc.terms.stack <- stack(u.Rocks, u.Texture[[2]], u.Texture[[3]], u.Texture[[4]],
                         u.PTF.FC[[2]], u.PTF.FC[[3]], u.PTF.FC[[4]], u.PTF.FC[[5]],
                         u.PTF.FC[[6]], u.PTF.FC[[7]],
                         u.PTF.PWP[[2]], u.PTF.PWP[[3]], u.PTF.PWP[[4]], u.PTF.PWP[[5]],
                         u.PTF.PWP[[6]], u.PTF.PWP[[7]], var.E.AWC.dec)
names(unc.terms.stack[["var.E.AWC.dec.100_200"]]) <- c("var.AWC.100_200")

levelplot(unc.terms.stack, par.settings = RdBuTheme) +
    layer(sp.polygons(france, lwd=0.7))
my_colors2 <- brewer.pal(11, "RdYlBu")

d2 <- c(seq(-0.007, -0.0005, by =0.002), -0.0005,
        -0.0001, 0, 0.0001,
        seq(0.0005, 0.005, by =0.002), 0.005)

### Test figure function
levelplot(unc.terms.stack, 
          margin=FALSE,                       
          colorkey=list(                                                    # suppress marginal graphics 
              space='right',                                               # legend ticks and labels                   
              labels=list(at= d2, font=3),  
              axis.line=list(col='black'),
              width=0.75),  
          par.settings=list(
              strip.border=list(col='transparent'),
              strip.background=list(col='transparent'),
              axis.line=list(col='transparent')),
          scales=list(draw=FALSE),   
          col.regions=my_colors2,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))

