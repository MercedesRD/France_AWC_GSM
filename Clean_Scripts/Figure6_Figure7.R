####### Packages
library(soiltexture)
library(plyr)
library(knitr)
library(Hmisc)
library(MASS)
library(wesanderson)
library(ggplot2)
library(fBasics)
library(gdata)
library(ithir)
library(gridExtra)
library(grid)
library(raster)
library(viridis)
library(viridisLite)

### Load packages, and install them if they are not in the library
list.of.packages <- c("soiltexture", "plyr", "MASS","knitr","fBasics", "gdata", "ithir","raster","sp","rgdal",
                      "gridExtra", "grid", "wesanderson","ggplot2", "car", "ggtern",
                      "GGally", "compositions", "rgr"   )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages, repos='http://cran.us.r-project.org')} 

### Set your working directory and load the file:
## setwd("D:/romandobarco/AWC/AWC_DSM/Clean_Output/9-validation.preds.GEVAROVIA")
load("9-Validation.preds.GEVARNOVIA.RData")

### Check we have all variables
names(bdd.validation.s)
names(bdd.validation.PTF)

### calculate upper and lower prediction interval limits
bdd.validation.s$SMFC.pred_UPL <- bdd.validation.s$SMFC.pred + (1.64 * sqrt(bdd.validation.s$SMFC.pred.var))
bdd.validation.s$SMFC.pred_LPL <- bdd.validation.s$SMFC.pred - (1.64 * sqrt(bdd.validation.s$SMFC.pred.var))

bdd.validation.s$SMPWP.pred_UPL <- bdd.validation.s$SMPWP.pred + (1.64 * sqrt(bdd.validation.s$SMPWP.pred.var))
bdd.validation.s$SMPWP.pred_LPL <- bdd.validation.s$SMPWP.pred - (1.64 * sqrt(bdd.validation.s$SMPWP.pred.var))

### calculate upper and lower prediction interval limits
bdd.validation.PTF$smFC.ptf <- predict(object=lm_w20_ClSd, newdata = bdd.validation.PTF,interval = "prediction", level=0.90 )
bdd.validation.PTF$smPWP.ptf <- predict(object=lm_w42_ClSd, newdata = bdd.validation.PTF,interval = "prediction", level=0.90)

bdd.validation.PTF$SMFC.ptf.pred <- bdd.validation.PTF$smFC.ptf[,1]
bdd.validation.PTF$SMFC.ptf_UPL <- bdd.validation.PTF$smFC.ptf[,3]
bdd.validation.PTF$SMFC.ptf_LPL <- bdd.validation.PTF$smFC.ptf[,2]

bdd.validation.PTF$SMPWP.ptf.pred <- bdd.validation.PTF$smPWP.ptf[,1]
bdd.validation.PTF$SMPWP.ptf_UPL <- bdd.validation.PTF$smPWP.ptf[,3]
bdd.validation.PTF$SMPWP.ptf_LPL <- bdd.validation.PTF$smPWP.ptf[,2]

### recalculate HYPRES texture class
gevarnovia.s_tri <- data.frame( "CLAY" = bdd.validation.s$clay, "SILT" = bdd.validation.s$silt,"SAND" = bdd.validation.s$sand)
gevarnovia.PTF_tri <- data.frame( "CLAY" = bdd.validation.PTF$clay, "SILT" = bdd.validation.PTF$silt,"SAND" = bdd.validation.PTF$sand)
bdd.validation.s$Texture <-TT.points.in.classes(tri.data = gevarnovia.s_tri,
                                                class.sys = "HYPRES.TT",
                                                PiC.type = "t", 
                                                tri.sum.tst = FALSE) 
bdd.validation.PTF$Texture <-TT.points.in.classes(tri.data = gevarnovia.PTF_tri,
                                                  class.sys = "HYPRES.TT",
                                                  PiC.type = "t", 
                                                  tri.sum.tst = FALSE)


########################################################################################################################################

library(scales)
library(viridis)
show_col(hue_pal()(5))


par(oma = c(4,4,4,15), mar = c(8,8,8,15))
smfc.dsm.plot <- ggplot(bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),],
                        aes(SMpF2.0, SMFC.pred, color=Texture))
p1 <- smfc.dsm.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
    labs(y = expression("Predicted SMFC (cm"^3*"/cm"^3*")"), 
         x = expression("Measured SMFC (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.6) + ylim(0,0.6)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    scale_colour_manual(values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
    #viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF2.0, y = SMFC.pred_LPL, xend = SMpF2.0, yend = SMFC.pred_UPL),
                 alpha = 5/10,size=3,
                 data = bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),])+
    geom_point(alpha = 8/10, size=9)
    
p1

############################################################################################################

smpwp.dsm.plot <- ggplot(bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),],
                        aes(SMpF4.2, SMPWP.pred, color=Texture))
p2 <- smpwp.dsm.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
    labs(y = expression("Predicted SMPWP (cm"^3*"/cm"^3*")"), 
         x = expression("Measured SMPWP (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.6) + ylim(0,0.6)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
       # viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF4.2, y = SMPWP.pred_LPL, xend = SMpF4.2, yend = SMPWP.pred_UPL),
                 alpha = 5/10,size=3,
                 data = bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),])+
    geom_point(alpha = 8/10, size=9)

p2

#############################################################################################################

smfc.ptf.plot <- ggplot(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),],
                        aes(SMpF2.0, SMFC.ptf.pred, color=Texture))
p3 <- smfc.ptf.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
    labs(y = expression("Predicted SMFC (cm"^3*"/cm"^3*")"), 
         x = expression("Measured SMFC (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.6) + ylim(0,0.6)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    scale_colour_manual(values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
    #viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF2.0, y = SMFC.ptf_UPL, xend = SMpF2.0, yend = SMFC.ptf_LPL),
                 alpha = 5/10,size=3,
                 data = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),])+
    geom_point(alpha = 8/10, size=9)

p3

plot (c(1,2,3,4,5,6,7), c(1,2,3,4,5,6,7),
      xlab=expression("Predicted "*theta[PWP]*" (cm"^3*"/cm"^3*")"), 
      ylab= expression("Measured "*theta[PWP]*" (cm"^3*"/cm"^3*")"))



############################################################################################################

smpwp.ptf.plot <- ggplot(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),],
                         aes(SMpF4.2, SMPWP.ptf.pred, color=Texture))
p4 <- smpwp.ptf.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
    labs(y = expression("Predicted SMPWP (cm"^3*"/cm"^3*")"), 
         x = expression("Measured SMPWP (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.6) + ylim(0,0.6)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
   # viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF4.2, y = SMPWP.ptf_LPL, xend = SMpF4.2, yend = SMPWP.ptf_UPL),
                 alpha = 5/10,size=3,
                 data = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),])+
    geom_point(alpha = 8/10, size=9)

p4

# Create a text
grob1 <- grobTree(textGrob("a)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob2 <- grobTree(textGrob("b)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob3 <- grobTree(textGrob("c)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob4 <- grobTree(textGrob("d)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))

##########################################################################
# Plot
jpeg("Figure6.4.jpeg", width = 3740 , height = 3600,units="px" )
grid.arrange(p1+ annotation_custom(grob1),
             p3+ annotation_custom(grob2),
             p2+ annotation_custom(grob3),
             p4+ annotation_custom(grob4),nrow=2)
dev.off()

bdd.validation.PTF[bdd.validation.PTF$Texture=="VF",]

###########################################################################################################################################
########################################################################################################################################
par(oma = c(4,4,4,15), mar = c(8,8,8,15))
smfc.dsm.plot <- ggplot(bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),],
                        aes(SMpF2.0, SMFC.pred, color=Horizon_profondeur_moyenne))
p1 <- smfc.dsm.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
    labs(y = expression("Predicted "*theta[FC]*" (cm"^3*"/cm"^3*")"), 
         x = expression("Measured "*theta[FC]*" (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.58) + ylim(0,0.58)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF2.0, y = SMFC.pred_LPL, xend = SMpF2.0, yend = SMFC.pred_UPL),
                 alpha = 6/10,size=3,
                 data = bdd.validation.s[!is.na(bdd.validation.s$SMpF2.0),])+
    geom_point(alpha = 8/10, size=9)

p1

############################################################################################################

smpwp.dsm.plot <- ggplot(bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),],
                         aes(SMpF4.2, SMPWP.pred, color=Horizon_profondeur_moyenne))
p2 <- smpwp.dsm.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
  labs(y = expression("Predicted "*theta[PWP]*" (cm"^3*"/cm"^3*")"), 
       x = expression("Measured "*theta[PWP]*" (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.58) + ylim(0,0.58)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF4.2, y = SMPWP.pred_LPL, xend = SMpF4.2, yend = SMPWP.pred_UPL),
                 alpha = 6/10,size=3,
                 data = bdd.validation.s[!is.na(bdd.validation.s$SMpF4.2),])+
    geom_point(alpha = 8/10, size=9)

p2

#############################################################################################################

smfc.ptf.plot <- ggplot(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),],
                        aes(SMpF2.0, SMFC.ptf.pred, color=Horizon_profondeur_moyenne))
p3 <- smfc.ptf.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
  labs(y = expression("Predicted "*theta[FC]*" (cm"^3*"/cm"^3*")"), 
       x = expression("Measured "*theta[FC]*" (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.58) + ylim(0,0.58)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF2.0, y = SMFC.ptf_UPL, xend = SMpF2.0, yend = SMFC.ptf_LPL),
                 alpha = 6/10,size=3,
                 data = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF2.0),])+
    geom_point(alpha = 8/10, size=9)

p3

############################################################################################################

smpwp.ptf.plot <- ggplot(bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),],
                         aes(SMpF4.2, SMPWP.ptf.pred, color=Horizon_profondeur_moyenne))
p4 <- smpwp.ptf.plot + geom_point(alpha = 9/10, size=9) + theme_bw() +
  labs(y = expression("Predicted "*theta[PWP]*" (cm"^3*"/cm"^3*")"), 
       x = expression("Measured "*theta[PWP]*" (cm"^3*"/cm"^3*")"))+
    theme(axis.text=element_text(size=60), 
          axis.title=element_text(size=70,face="bold"),
          legend.text=element_text(size=rel(4)),
          legend.title=element_text(size=rel(5)),
          legend.key.height = unit(30,"points"),
          legend.key.width = unit(100, "points"))+
    xlim(0, 0.58) + ylim(0,0.58)+
    geom_abline(slope=1, intercept = 0, size=2, colour="black")+
    viridis::scale_color_viridis(option = "A",discrete = FALSE ,direction = -1, name="Mean horizon depth") +
    theme(legend.position="bottom") +
    theme(plot.margin = unit(c(0.5,2,0.5,0.5), "cm"))+
    geom_segment(aes(x = SMpF4.2, y = SMPWP.ptf_LPL, xend = SMpF4.2, yend = SMPWP.ptf_UPL),
                 alpha = 6/10,size=3,
                 data = bdd.validation.PTF[!is.na(bdd.validation.PTF$SMpF4.2),])+
    geom_point(alpha = 8/10, size=9)

p4

# Create a text
grob1 <- grobTree(textGrob("a)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob2 <- grobTree(textGrob("b)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob3 <- grobTree(textGrob("c)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))
grob4 <- grobTree(textGrob("d)", x=0.05,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=60)))

##########################################################################
# Plot
library(gridExtra)
jpeg("Figure6.jpeg", width = 3740 , height = 3600,units="px" )
grid.arrange(p1+ annotation_custom(grob1),
             p3+ annotation_custom(grob2),
             p2+ annotation_custom(grob3),
             p4+ annotation_custom(grob4),nrow=2)
dev.off()

#################################################################################################################################

#################################################################################################################################

### calculate the prediction errors
bdd.validation.s$SMFC.GSM.error <- bdd.validation.s$SMFC.pred - bdd.validation.s$SMpF2.0
bdd.validation.s$SMPWP.GSM.error <- bdd.validation.s$SMPWP.pred - bdd.validation.s$SMpF4.2

bdd.validation.PTF$SMFC.PTF.error <- bdd.validation.PTF$SMFC.ptf.pred - bdd.validation.PTF$SMpF2.0
bdd.validation.PTF$SMPWP.PTF.error <- bdd.validation.PTF$SMPWP.ptf.pred - bdd.validation.PTF$SMpF4.2

bdd.validation.s$Texture <- factor(bdd.validation.s$Texture,
                                   levels = c("C", "M", "MF", "F", "VF"))
bdd.validation.PTF$Texture <- factor(bdd.validation.PTF$Texture,
                                     levels = c("C", "M", "MF", "F", "VF"))

### Figure 7
par()
par(oma = c(4,4,4,15), mar = c(8,8,8,15))
jpeg("Figure7.jpeg", width = 2756 , height = 2500,units="px" )
jpeg("Figure7.jpeg", width = 800 , height = 700,units="px" )
par(mfrow=c(2,2))
par(oma = c(2,2,2,0), mar = c(4,2,2,2))
boxplot(bdd.validation.s$SMFC.GSM.error ~ bdd.validation.s$Texture, col="grey", las=1,
        main=expression("GSM "*theta[FC]*" error"), cex.main=2, cex.axis=1.5,
        ylim=c(-0.28,0.28))
abline(0,0, lty=3)
text("a)", x=0.5, y=0.25, cex=1.7)

boxplot(bdd.validation.PTF$SMFC.PTF.error ~ bdd.validation.PTF$Texture, col="grey", las=1,
        main=expression("PTF "*theta[FC]*" error"), cex.main=2, cex.axis=1.5,
        ylim=c(-0.28,0.28))
abline(0,0, lty=3)
text("b)", x=0.5, y=0.25, cex=1.7)

boxplot(bdd.validation.s$SMPWP.GSM.error ~ bdd.validation.s$Texture, col="grey", las=1,
        main=expression("GSM "*theta[PWP]*" error"), cex.main=2, cex.axis=1.5,
        ylim=c(-0.28,0.28))
abline(0,0, lty=3)
text("c)", x=0.5, y=0.25, cex=1.7)

boxplot(bdd.validation.PTF$SMPWP.PTF.error ~ bdd.validation.PTF$Texture, col="grey", las=1,
        main=expression("PTF "*theta[PWP]*" error"), cex.main=2, cex.axis=1.5,
        ylim=c(-0.28,0.28))
abline(0,0, lty=3)
text("d)", x=0.5, y=0.25, cex=1.7)
dev.off()

### end of the script
