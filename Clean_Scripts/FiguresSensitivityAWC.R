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
library(gridExtra)


### Load the data
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.0_5/AWC.var.decomp.0_5.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.5_15/AWC.var.decomp.5_15.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.15_30/AWC.var.decomp.15_30.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.30_60/AWC.var.decomp.30_60.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.60_100/AWC.var.decomp.60_100.RData")
load("D:/romandobarco/AWC/AWC_GSM_Dec2017/Output/Recalculation_AWC/SensitivityAWC_general/varAWC.100_200/AWC.var.decomp.100_200.RData")

### Load France
france<- readOGR(dsn="Y:/communs_infosol/Projets/GlobalSoilMap.net/GIS/administratif", layer = "france_L93_ncorse")
france <- spTransform(france, CRS("+init=epsg:2154"))

### does not work because each variable has a different scale. Need to plot one by one
my_rasters <- list(sens.Rocks, sens.Clay, sens.Silt, sens.Beta0, sens.Beta1, sens.Beta2)
plots <- list()
my_colors2 <- brewer.pal(11, "RdYlBu")
for(i in 1:length(my_rasters)){
    
    pi <-  levelplot(my_rasters[[i]],
                     cuts=10,
                     margin =FALSE, 
                     colorkey=list( space='right'),  
                     par.settings=list(
                         strip.border=list(col='transparent'),
                         strip.background=list(col='transparent'),
                         axis.line=list(col='transparent')),
                     col.regions=my_colors2,
                     scales=list(draw=FALSE)) +
        layer(sp.polygons(france, lwd=0.7))
    plots[[i]] <- pi
}


grid.arrange(plots[[1]],plots[[2]],plots[[3]],
             plots[[4]],plots[[5]],plots[[6]],ncol=3)

##################################################################################################################

my_colors <- viridis_pal(option="D", direction=-1)(29)
my_colors2 <- brewer.pal(11, "RdYlBu")
d2 <- c(seq(0, 0.0001, by =0.00001),
        seq(0.0002, 0.002, by =0.0001))
levelplot(u.rast,par.settings = infernoTheme)

u <- stack(u.Rocks, u.Texture[[1]], u.PTF.FC[[1]], u.PTF.PWP[[1]])
levelplot(u.rast, 
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
          col.regions=my_colors,                   
          at= d2)   +                         # colour ramp breaks     
    layer(sp.polygons(france, lwd=0.6))
