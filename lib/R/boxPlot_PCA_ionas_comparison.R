#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. April 2013                       ###
#############################################################################
### Phecomp                                                               ###
### Plot weights of mice in boxplot by week                               ### 
#############################################################################

##Loading libraries
library (ggplot2)
library (plyr)
library(reshape)
library (grid) #viewport

##Getting HOME directory
home <- Sys.getenv("HOME")

# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))

# All acquisition variables used to perform the PCA 
# Here I load the first PC
load("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/5setsPC1.R")
# Gallagher used to perform the PCA. Here I load the first PC
#load("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/8setsGall.R")

TS
TSEE
TSEEEGCG
TSEGCG
WT

# TS vs TSEEEGCG
df.TS <- as.data.frame(TS)
id <- c(1:length(df.TS$V1))
colnames (df.TS) <- c("A1", "A2", "A3", "A4", "A5")
df.TS.m <- melt(df.TS)
df.TS.m$group <- "TS"
df.TS.m$id <- id
# TSEEEGCG
df.TSEEEGCG <- as.data.frame(TSEEEGCG)
colnames (df.TSEEEGCG) <- c("A1", "A2", "A3", "A4", "A5")
id_TSEEEGCG <- c((length(df.TS$A1) + 1) : (length(df.TS$A1) + 1 +length(df.TSEEEGCG$A1)-1))
df.TSEEEGCG.m <- melt(df.TSEEEGCG)
df.TSEEEGCG.m$group <- "TSEEEGCG"
df.TSEEEGCG.m$id <- id_TSEEEGCG
df.plot <- rbind(df.TS.m, df.TSEEEGCG.m)
#TSEE
df.TSEE <- as.data.frame(TSEE)
colnames (df.TSEE) <- c("A1", "A2", "A3", "A4", "A5")
id_TSEE <- c(((id_TSEEEGCG[length(id_TSEEEGCG)]+1) : (id_TSEEEGCG[length(id_TSEEEGCG)] + length(df.TSEE$A1))))
df.TSEE.m <- melt(df.TSEE)
df.TSEE.m$group <- "TSEE"
df.TSEE.m$id <- id_TSEE
df.anova <- rbind(df.plot, df.TSEE.m)
#TSEGCG
df.TSEGCG <- as.data.frame(TSEGCG)
colnames (df.TSEGCG) <- c("A1", "A2", "A3", "A4", "A5")
id_TSEGCG <- c(((id_TSEE[length(id_TSEE)]+1) : (id_TSEE[length(id_TSEE)] + length(df.TSEGCG$A1))))
df.TSEGCG.m <- melt(df.TSEGCG)
df.TSEGCG.m$group <- "TSEGCG"
df.TSEGCG.m$id <- id_TSEGCG
df.anova <- rbind(df.anova, df.TSEGCG.m)
#WT
df.WT <- as.data.frame(WT)
colnames (df.WT) <- c("A1", "A2", "A3", "A4", "A5")
id_WT <- c(((id_TSEGCG[length(id_TSEGCG)]+1) : (id_TSEGCG[length(id_TSEGCG)] + length(df.WT$A1))))
df.WT.m <- melt(df.WT)
df.WT.m$group <- "WT"
df.WT.m$id <- id_WT
df.anova <- rbind(df.anova, df.WT.m)

# Adding wt guys that for the anova of the first PC1
load("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/3setsPC1.R")
WTEE
WTEGCG
WTEEEGCG

#WT
df.WTEE <- as.data.frame(WTEE)
colnames (df.WTEE) <- c("A1", "A2", "A3", "A4", "A5")
id_WTEE <- c(((id_WT[length(id_WT)]+1) : (id_WT[length(id_WT)] + length(df.WTEE$A1))))
df.WTEE.m <- melt(df.WTEE)
df.WTEE.m$group <- "WTEE"
df.WTEE.m$id <- id_WTEE
df.anova <- rbind(df.anova, df.WTEE.m)

#WTEGCG
df.WTEGCG <- as.data.frame(WTEGCG)
colnames (df.WTEGCG) <- c("A1", "A2", "A3", "A4", "A5")
id_WTEGCG <- c(((id_WTEE[length(id_WTEE)]+1) : (id_WTEE[length(id_WTEE)] + length(df.WTEGCG$A1))))
df.WTEGCG.m <- melt(df.WTEGCG)
df.WTEGCG.m$group <- "WTEGCG"
df.WTEGCG.m$id <- id_WTEGCG
df.anova <- rbind(df.anova, df.WTEGCG.m)

#WTEEEGCG
df.WTEEEGCG <- as.data.frame(WTEEEGCG)
colnames (df.WTEEEGCG) <- c("A1", "A2", "A3", "A4", "A5")
id_WTEEEGCG <- c(((id_WTEGCG[length(id_WTEGCG)]+1) : (id_WTEGCG[length(id_WTEGCG)] + length(df.WTEEEGCG$A1))))
df.WTEEEGCG.m <- melt(df.WTEEEGCG)
df.WTEEEGCG.m$group <- "WTEEEGCG"
df.WTEEEGCG.m$id <- id_WTEEEGCG
df.anova <- rbind(df.anova, df.WTEEEGCG.m)

### PLOT 2 groups

boxPlots <- ggplot(df.plot , aes (variable, value, fill = group, color=group)) + 
  geom_boxplot(show_guide=FALSE) + 
  scale_color_manual(values = c("black", "brown")) +
  scale_fill_manual(name = "Group", values = c("green", "black"), labels = c ("TS", "TSEEEGCG")) +
  #scale_fill_manual (name = "Group", values = c ("green", "brown")) +
  labs(title = "") + xlab ("\nAcquisition days") + ylab("PC1\n")

boxPlots + annotate("text", x=4, y=6, label="*", size=10) + 
           annotate("text", x=5, y=7.6, label="*", size=10)



## ANOVA
demo1 <- read.csv("http://www.ats.ucla.edu/stat/data/demo1.csv")
class(demo1$time)
# df.plot$time <- as.integer(df.plot$variable)
# 
# df.plot [with(df.plot, order(id)), ]
## Convert variables to factor
df.anova <- within(df.anova, {
  group <- factor(group)
  time <- factor(variable)
  id <- factor(id)
})

df.anova

demo1.aov <- aov(value ~ group * time + Error(id), data = df.anova)
# demo1.aov <- aov(pulse ~ group * time + Error(id), data = demo1)
# I set the interaction to perform the ttest
df.anova$interaction <- paste(df.anova$group, df.anova$time, sep="_")
pairwise.t.test(df.anova$value, df.anova$interaction, , p.adj="hochberg", paired=F)


# Comparison only TS vs TSEEEGCG
df_TS_TSEEEGCG <- rbind(df.TS.m, df.TSEEEGCG.m)
df_TS_TSEEEGCG$interaction <- paste (df_TS_TSEEEGCG$variable, df_TS_TSEEEGCG$group,sep="_")
pairwise.t.test (df_TS_TSEEEGCG$value, df_TS_TSEEEGCG$interaction, , p.adj="bonferroni", paired=F)

length(df.anova$value)
length(df.anova$interaction)
?pairwise.t.test

summary(demo1.aov)
class(df.plot$value)
length(Group)
length(Value)

# Example
# Group <- c("A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B", "C","C","C","C","C","C","C","C") 
# Value <- c(1,2,4,1,1,2,2,3,3,4,4,2,3,4,4,3,4,5,3,5,5,3,4,6) 
# Participant <- c("1","2","3","4","5","6","7","8","1","2","3","4","5","6","7","8", "1","2","3","4","5","6","7","8") 
# pairwise.t.test(Value, Group, p.adj="bonferroni", paired=T)
# 
# data <- data.frame(Participant, Group, Value) 
# aov <- aov(Value ~ factor(Group) + Error(factor(Participant)/factor(Group)), data) 
# summary(aov)

# demo1 <- read.csv("http://www.ats.ucla.edu/stat/data/demo1.csv")
# ## Convert variables to factor
# demo1 <- within(demo1, {
#   group <- factor(group)
#   time <- factor(time)
#   id <- factor(id)
# })

# Exporting the data to excel
df.WT$group <- "WT"
df.WTEE$group <- "WTEE"
df.WTEGCG$group <- "WTEGCG"
df.WTEEEGCG$group <- "WTEEEGCG"
df.TS$group <- "TS"
df.TSEE$group <- "TSEE"
df.TSEGCG$group <- "TSEGCG"
df.TSEEEGCG$group <- "TSEEEGCG"

df.spss <- rbind (df.WT, df.WTEE, df.WTEGCG, df.WTEEEGCG, df.TS, df.TSEE, df.TSEGCG, df.TSEEEGCG)

library(xlsx)
# write.xlsx (df.spss, paste(home, "/20150515_PCA_old_frotiersPaper/data/anova_PC1.xlsx", sep="")) 
# write.xlsx (df.spss, paste(home, "/sharedWin/anova_PC1.xlsx", sep=""))

boxPlots <- ggplot(df.anova , aes (variable, value, fill = group, color=group)) + 
#   geom_boxplot(show_guide=FALSE) + 
  geom_boxplot() + 
  scale_color_manual(values = c("black", "brown","black", "black","black","black", "black","black")) +
  scale_fill_manual(name = "Group", values = c("green", "black","red","blue","pink","yellow", "purple", "gray", "brown")) +
  #scale_fill_manual (name = "Group", values = c ("green", "brown")) +
  labs(title = "") + xlab ("\nAcquisition days") + ylab("PC1\n")
boxPlots

#### 
# Comparison of only day 5 and TS
df.anova

df.anova.ts<-subset(df.anova, grepl("TS", df.anova$group))
df.anova.ts.a5 <- subset(df.anova.ts, grepl("A5", df.anova.ts$variable))

ts_a5.aov <- aov(value ~ group  + Error(id), data = df.anova.ts.a5)
summary(ts_a5.aov)

pairwise.t.test (df.anova.ts.a5$value, df.anova.ts.a5$group, p.adj="bonferroni")
pairwise.t.test (df.anova.ts.a5$value, df.anova.ts.a5$group, p.adj="hochberg")


uniqueInitials <- c("red", "lightblue","lightgreen", "darkgreen")
initialShapes <- unlist(lapply(uniqueInitials, utf8ToInt))

# Setting order for ploting
df.anova.ts.a5$group <- factor(df.anova.ts.a5$group , levels=c("TS","TSEE", "TSEGCG", "TSEEEGCG"), 
                             labels=c("TS","TSEE", "TSEGCG", "TSEEEGCG"))



# dailyInt_theme <- theme_update (panel.border = element_rect(colour = "black"))
# boxPlots <- ggplot(df.anova.ts.a5 , aes (variable, value, fill = group, color=group)) +

boxPlots <- ggplot(df.anova.ts.a5 , aes (group, value, fill = group)) + 
#                    geom_boxplot() +
                   geom_boxplot(show_guide=FALSE) +
#             guides(color=guide_legend('Model',override.aes=list(shape=c(1,1,6,6))))
  
  scale_fill_manual(name = "Group", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Session 5 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-5.5, 9.5)) +
  geom_segment(aes(x = 3.63, y = median(df.anova.ts.a5[df.anova.ts.a5$group == "TSEEEGCG","value"]), 
                 xend = 4.37, yend = median(df.anova.ts.a5[df.anova.ts.a5$group == "TSEEEGCG","value"])), colour="white")

boxPlots 

#PLOT_paper
ggsave (boxPlots, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_ts_a5.jpg", sep=""), dpi=900)

###########################################
# Comparison of only day 5 and TS
df.anova.ts

df.anova.ts.a1 <- subset(df.anova.ts, grepl("A1", df.anova.ts$variable))

ts_a1.aov <- aov(value ~ group  + Error(id), data = df.anova.ts.a1)
summary(ts_a1.aov)

# pairwise.t.test (df.anova.ts.a1$value, df.anova.ts.a1$group, p.adj="bonferroni")
# pairwise.t.test (df.anova.ts.a1$value, df.anova.ts.a1$group,p.adj="hochberg")

uniqueInitials <- c("red", "lightblue","lightgreen", "darkgreen")
initialShapes <- unlist(lapply(uniqueInitials, utf8ToInt))

# Setting order for ploting
df.anova.ts.a1$group <- factor(df.anova.ts.a5$group , levels=c("TS","TSEE", "TSEGCG", "TSEEEGCG"), 
                               labels=c("TS","TSEE", "TSEGCG", "TSEEEGCG"))



# dailyInt_theme <- theme_update (panel.border = element_rect(colour = "black"))
# boxPlots <- ggplot(df.anova.ts.a5 , aes (variable, value, fill = group, color=group)) +
med_tseeegcg<-median(df.anova.ts.a1[df.anova.ts.a1$group == "TSEEEGCG","value"])
boxPlots <- ggplot(df.anova.ts.a1 , aes (group, value, fill = group)) + 
#                      geom_boxplot() +
   geom_boxplot(show_guide=FALSE) +
  #             guides(color=guide_legend('Model',override.aes=list(shape=c(1,1,6,6))))
  
  scale_fill_manual(name = "Group", values=c("darkgreen", "lightblue", "orange", "black")) +
#   scale_colour_manual (name="Group",values = c(rep("gray",4))) +
  labs(title = "Session 1 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
#   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-5.5, 9.5)) +
  geom_segment(aes(x = 3.63, y = median(df.anova.ts.a1[df.anova.ts.a1$group == "TSEEEGCG","value"]), xend = 4.37, yend = median(df.anova.ts.a1[df.anova.ts.a1$group == "TSEEEGCG","value"])), colour="white")
                                 
boxPlots

#PLOT_paper
# ggsave (boxPlots, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_ts_a1.jpg", sep=""), dpi=900)

# Plotting a legend with squares
l <- ggplot() + geom_point(data=df.anova.ts.a1 , aes (x=group, y=value, colour = group), shape=15, size=5) +
     scale_colour_manual (values=c("darkgreen", "lightblue", "orange", "black"))
l <- l + guides(color=guide_legend(title=NULL)) 
l <- l + theme(legend.key = element_blank())
l
# ggsave (l, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_Legend_TS.jpg", sep=""), dpi=900)

################
# Boxplots of wt
################

df.anova.wt <- subset(df.anova, grepl("WT", df.anova$group))
df.anova.wt.a5 <- subset(df.anova.wt, grepl("A5", df.anova.ts$variable))
df.anova.wt.a1 <- subset(df.anova.wt, grepl("A1", df.anova.ts$variable))

boxPlots_wt_a1 <- ggplot(df.anova.wt.a1 , aes (group, value, fill = group)) + 
  #                      geom_boxplot() +
  geom_boxplot(show_guide=FALSE) +
  #             guides(color=guide_legend('Model',override.aes=list(shape=c(1,1,6,6))))
  
  scale_fill_manual(name = "Group", values = c("red", "blue", "magenta",  "yellow")) + 
#                     values=c("darkgreen", "lightblue", "orange", "black")) +
  #   scale_colour_manual (name="Group",values = c(rep("gray",4))) +
  labs(title = "Session 1 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6,8), limits=c(-6.5, 9.5))
#   geom_segment(aes(x = 3.63, y = median(df.anova.wt.a1[df.anova.wt.a1$group == "TSEEEGCG","value"]), xend = 4.37, yend = median(df.anova.wt.a1[df.anova.wt.a1$group == "WTEEEGCG","value"])), colour="white")

boxPlots_wt_a1

boxPlots_wt_a5 <- ggplot(df.anova.wt.a5 , aes (group, value, fill = group)) + 
  #                      geom_boxplot() +
  geom_boxplot(show_guide=FALSE) +
  #             guides(color=guide_legend('Model',override.aes=list(shape=c(1,1,6,6))))
  
  scale_fill_manual(name = "Group", values = c("red", "blue", "magenta",  "yellow")) + 
  #                     values=c("darkgreen", "lightblue", "orange", "black")) +
  #   scale_colour_manual (name="Group",values = c(rep("gray",4))) +
  labs(title = "Session 5 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6,8), limits=c(-6.5, 9.5))
#   geom_segment(aes(x = 3.63, y = median(df.anova.wt.a1[df.anova.wt.a1$group == "TSEEEGCG","value"]), xend = 4.37, yend = median(df.anova.wt.a1[df.anova.wt.a1$group == "WTEEEGCG","value"])), colour="white")

boxPlots_wt_a5

#PLOT_paper

# ggsave (boxPlots_wt_a1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_wt_a1.jpg", sep=""), dpi=900)
# ggsave (boxPlots_wt_a5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_wt_a5.jpg", sep=""), dpi=900)

# Plotting a legend with squares
l <- ggplot() + geom_point(data=df.anova.wt.a1 , aes (x=group, y=value, colour = group), shape=15, size=5) +
  scale_colour_manual (values=c("red", "blue", "magenta",  "yellow"))
l <- l + guides(color=guide_legend(title=NULL)) 
l <- l + theme(legend.key = element_blank())
l
# ggsave (l, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_Legend_WT.jpg", sep=""), dpi=900)
















##################################



#                     guides(colour = guide_legend(override.aes = list(shape = 11)))

#scale_color_manual(values = c("green", "black","red", "blue")) +
#              scale_fill_manual(name = "Group", values = c("red", "lightblue","lightgreen", "darkgreen"))+
#                     , labels = c ("TS", "TSEEEGCG")) +
  #scale_fill_manual (name = "Group", values = c ("green", "brown")) +

#             theme(legend.text = element_text(colour=c("red", "lightblue","lightgreen", "darkgreen"), size = 16, face = "bold"))
#             theme(element_text(colour=c("red", "lightblue","lightgreen", "darkgreen"), size = 16, face = "bold"))
            guides(colour = guide_legend(override.aes = list(values=initialShapes)))+
        
#               theme (legend.text = element_text(colour = c('red', "lightblue","lightgreen", "darkgreen")))+
            theme(legend.key=element_rect(fill=NA))






          scale_shape_manual(values = initialShapes) 
       
              guides(fill = guide_legend(override.aes = list(values=initialShapes)))
              
guides(fill = guide_legend(override.aes = list(linetype = 0, shape='')), 
                   colour = guide_legend(override.aes = list(linetype=c(0,1,1)
                                                     , shape=c(16,NA,NA))))

#   theme (legend.text = element_text(colour=c("red", "lightblue","lightgreen", "darkgreen"), size = 16, face = "bold"))

boxPlots + scale_shape_manual(values = initialShapes) 
boxPlots$theme$legend.text$colour<-c("red", "lightblue","lightgreen", "darkgreen")
boxPlots + annotate("text", x=4, y=6, label="*", size=10) + 
  annotate("text", x=5, y=7.6, label="*", size=10)
                 