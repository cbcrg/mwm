setwd("/Users/jespinosa/MWM/20150811_graficosPosterSilvina/")

#neuros <- spss.get("Datos_neurodegeneracion_colinergica_red.sav")
neuros <- spss.get("/Users/jespinosa/MWM/20150811_graficosPosterSilvina/Datos_neurodegeneracion_colinergica 18_06_15.sav")
neuros$mouse <- as.character(neuros$mouse)
##
library(beeswarm)
windows(width=14)
par(las=1, font=2, font.lab=2, font.axis=2)
boxplot(ncells~gentreat, subset(neuros, reg=='MS' & ncells < 4500), col=rainbow(8), pch=16)
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
boxplot(celldens~gentreat, subset(neuros, reg=='MS'), col=rainbow(8), pch=16)
beeswarm(celldens~gentreat, subset(neuros, reg=='MS'), add=T, pch=16)
title('Region: MS; Cholinergic neuronal density')
savePlot('MsCelldens', type='png')

####################
# Code Jose here
####################
head(neuros)

# only rg=MS
neuros_MS <- neuros [which(neuros$reg == "MS"),]
neuros_MS$gentreat <- gsub("_", "", neuros$gentreat)
neuros_MS$gentreat <- factor(neuros_MS$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                                labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

max(neuros_MS$celldens)

# boxplot neuro_MS
bp_neuros_MS <- ggplot(neuros_MS , aes (gentreat, celldens, fill = gentreat)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Region: MS; Cholinergic neuronal density\n") + xlab ("\nGroups") + ylab("Number of neurons\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(4750, 7700), breaks = seq(5000,7500, by=500)) +
  geom_segment(aes(x = 7.63, y = as.numeric(median(neuros_MS [neuros_MS$gentreat == "TSEEEGCG", "celldens"])), 
                   xend = 8.37, yend = as.numeric(median(neuros_MS [neuros_MS$gentreat == "TSEEEGCG","celldens"]))), 
               colour="white", size =0.8)


bp_neuros_MS
bp_neuros_MS_points <- bp_neuros_MS + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_neuros_MS_points

ggsave (bp_neuros_MS_points, file="chol_neurons_MS.jpg", width=14, height=7, dpi=900)





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

