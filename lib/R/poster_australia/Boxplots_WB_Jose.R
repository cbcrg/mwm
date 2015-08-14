

library(Hmisc)
library(doBy)
library(nlme)
library(multcomp)
library(Hmisc)
library(beeswarm)

install.packages ('ggplot2')
install.packages ('doBy')
  
install.packages("ggplot2")
library(ggplot2)

setwd("F:/1. Lab_CRG/Statistics Klaus/WB APPCTF/")
APP <- spss.get("F:/1. Lab_CRG/Statistics Klaus/WB APPCTF/Datos_WB_APP_CTF.sav")

app <- read.table("F:/1. Lab_CRG/Statistics Klaus/WB APPCTF/WB_APP_CTF.csv", header=TRUE, 
                     sep=",", row.names="id")
####################
# Code Jose here
####################
head(APP)

#ES este!!!!
bp_APP <- ggplot(APP , aes (Gen.treatment, Ratio.A.B.APPfull, fill = Gen.treatment)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Ratio ??/?? secretase/APP full length\n") + xlab ("\nGroups") + ylab("A.U.\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(1, 5), breaks = seq(1,5, by=0.5)) +
  geom_segment(aes(x = 7.63, y = as.numeric(median(APP [APP$Gen.treatment == "TSEEEGCG", "Ratio_A_B_APPfull"])), 
                   xend = 8.37, yend = as.numeric(median(APP [APP$Gen.treatment == "TSEEEGCG","Ratio_A_B_APPfull"]))), 
               colour="white", size =0.8)

bp_APP 
bp_APP_points <- bp_APP  + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_APP_points

ggsave (bp_APP_points, file="WB_APP_CTF.jpg", width=14, height=7, dpi=900)

#################################################################

#neuros <- spss.get("Datos_neurodegeneracion_colinergica_red.sav")
app <- spss.get("F:/1. Lab_CRG/Statistics Klaus/WB APPCTF/Datos_WB_APP_CTF.sav")
neuros$mouse <- as.character(neuros$mouse)
##
library(beeswarm)
windows(width=14)
par(las=1, font=2, font.lab=2, font.axis=2)
boxplot(Ratio.A.B.APPfull~gentreat, subset(app, reg=='MS' & ncells < 4500), col=rainbow(8), pch=16)
beeswarm(ncells~gentreat, subset(neuros, reg=='MS'), add=T, pch=16)
title('Region: MS; Number of cells')
savePlot('MsNcells', type='png')

windows(width=14)
par(las=1, font=2, font.lab=2, font.axis=2)
boxplot(ncells~gentreat, subset(neuros, reg!='MS'), col=rainbow(8), pch=16)
beeswarm(ncells~gentreat, subset(neuros, reg!='MS'), add=T, pch=16)
title('Region: VDB; Number of cells')
savePlot('VdbNcells', type='png')


####################
# Code figure here
####################

windows(width=14)
par(las=1, font=2, font.lab=2, font.axis=2)
boxplot(Ratio.A.B.APPfull~gentreat, subset(neuros, reg=='MS'), col=rainbow(8), pch=16)
beeswarm(Ratio.A.B.APPfull~gentreat, subset(neuros, reg=='MS'), add=T, pch=16)
title('WB_APP')
savePlot('WB_APP', type='png')

####################
# Code Jose here
####################
head(app)

# only rg=MS
neuros_MS <- neuros [which(neuros$reg == "MS"),]
neuros_MS$gentreat <- gsub("_", "", neuros$gentreat)
app$gentreat <- factor(app$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                                labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

max(neuros_MS$celldens)

# boxplot neuro_MS
bp_neuros_MS <- ggplot(app , aes (gentreat, Ratio.A.B.APPfull, fill = gentreat)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "WB_APP\n") + xlab ("\nGroups") + ylab("Number of neurons\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(4750, 7700), breaks = seq(5000,7500, by=500)) +
  geom_segment(aes(x = 7.63, y = as.numeric(median(app [app$gentreat == "TSEEEGCG", "Ratio.A.B.APPfull"])), 
                   xend = 8.37, yend = as.numeric(median(app [app$gentreat == "TSEEEGCG","Ratio.A.B.APPfull"]))), 
               colour="white", size =0.8)


bp_app
bp_app_points <- bp_app + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_app_points

ggsave (bp_app, file="WB_APP.jpg", width=14, height=7, dpi=900)





windows(width=14)
par(las=1, font=2, font.lab=2, font.axis=2)
boxplot(celldens~gentreat, subset(neuros, reg!='MS'), col=rainbow(8), pch=16)
beeswarm(celldens~gentreat, subset(neuros, reg!='MS'), add=T, pch=16)
title('Region: VDB; Cholinergic neuronal density')
savePlot('VdbCelldens', type='png')

graphics.off()

##
lisi <- list(c(1:4), c(1, 2, 5, 6), c(1, 2, 7, 8))
for(i in 1:3){
  windows(width=10)
  par(las=1, font=2, font.lab=2, font.axis=2) 
  boxplot(ncells~factor(gentreat), subset(neuros, reg=='MS' & ncells < 4500 & as.numeric(gentreat) %in% lisi[[i]]), col=rainbow(4), pch=16)
  beeswarm(ncells~factor(gentreat), subset(neuros, reg=='MS' & ncells < 4500 & as.numeric(gentreat) %in% lisi[[i]]), add=T, pch=16)
  title('Region: MS; Number of cells')
  savePlot(paste0("MsNcellsSubgroups", i))
  dev.off()
}
rm(lisi)

## Some preliminary models
detach(package:Hmisc)
denslmM <- lm(celldens~gentreat, data=subset(neuros, reg=='MS'))
anova(denslmM)
par(mfrow=c(2,2), font=2, font.axis=2, font.lab=2, pch=16)
plot(denslmM)

ncellmM <- lm(ncells~gentreat, data=subset(neuros, reg=='MS' & ncells<4500), weights=1/sqrt(percerro))
anova(ncellmM )
par(mfrow=c(2,2), font=2, font.axis=2, font.lab=2, pch=16)
plot(ncellmM )

