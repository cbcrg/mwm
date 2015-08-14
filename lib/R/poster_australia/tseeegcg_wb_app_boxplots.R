library(Hmisc)
library(doBy)
library(nlme)
library(multcomp)
library(Hmisc)
library(beeswarm)
library(ggplot2)

#setwd("F:/1. Lab_CRG/Statistics Klaus/WB APPCTF/")
app <- read.table("/Users/jespinosa/MWM/20150811_graficosPosterSilvina/WB_APPCTF.csv", header=TRUE, sep=";", dec=",")

tail (app)

#Groups inconsistenly named
app$Gen_treatment_label <- gsub("_", "", app$Gen_treatment_label)
app$Gen_treatment_label <- gsub("-", "", app$Gen_treatment_label)
app$Gen_treatment_label <- factor(app$Gen_treatment_label , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                             labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))


# Limits for y axis
min(app$Ratio_A_B_APPfull)
max(app$Ratio_A_B_APPfull)

app [which (app$Gen_treatment_label == "TS" | app$Gen_treatment_label == "WT"),]

# boxplot 
bp_APP <- ggplot(app , aes (Gen_treatment_label, Ratio_A_B_APPfull, fill = Gen_treatment_label)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Ratio ??/?? secretase/APP full length\n") + xlab ("\nGroups") + ylab("A.U.\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(1, 2), breaks = seq(0,2, by=0.2)) +
  geom_segment(aes(x = 7.63, y = as.numeric(median(app [app$Gen_treatment_label == "TSEEEGCG", "Ratio_A_B_APPfull"])), 
                   xend = 8.37, yend = as.numeric(median(app [app$Gen_treatment_label == "TSEEEGCG","Ratio_A_B_APPfull"]))), 
               colour="white", size =0.4)

bp_APP 

bp_APP_points <- bp_APP  + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_APP_points

ggsave (bp_APP_points, file="WB_APP_CTF.jpg", width=14, height=7, dpi=900)
