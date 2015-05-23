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
source (paste (home, "/git/phecomp/lib/R/plotParamPublication.R", sep=""))

load("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/5setsPC1.R")
load("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/8setsGall.R")
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

with(df.anova, interaction.plot(df.anova, group, time,
                              ylim = c(5, 20), lty= c(1, 12), lwd = 3,
                             ylab = "mean of pulse", xlab = "time", trace.label = "group"))

demo1.aov <- aov(value ~ group * time + Error(id), data = df.anova)
demo1.aov <- aov(pulse ~ group * time + Error(id), data = demo1)
df.anova$interaction <- paste(df.anova$group, df.anova$time, sep="_")
  pairwise.t.test(df.anova$value, df.anova$interaction, , p.adj="bonferroni", paired=T)
length(df.anova$value)
length(df.anova$interaction)
?pairwise.t.test

summary(demo1.aov)
class(df.plot$value)
length(Group)
length(Value)

pairwise.t.test(Value, Group, p.adj="bonferroni", paired=T)
Group <- c("A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B", "C","C","C","C","C","C","C","C") 
Value <- c(1,2,4,1,1,2,2,3,3,4,4,2,3,4,4,3,4,5,3,5,5,3,4,6) 
Participant <- c("1","2","3","4","5","6","7","8","1","2","3","4","5","6","7","8", "1","2","3","4","5","6","7","8") 
data <- data.frame(Participant, Group, Value) 
aov <- aov(Value ~ factor(Group) + Error(factor(Participant)/factor(Group)), data) 
summary(aov)

demo1 <- read.csv("http://www.ats.ucla.edu/stat/data/demo1.csv")
## Convert variables to factor
demo1 <- within(demo1, {
  group <- factor(group)
  time <- factor(time)
  id <- factor(id)
})


boxPlots <- ggplot(df.anova , aes (variable, value, fill = group, color=group)) + 
#   geom_boxplot(show_guide=FALSE) + 
  geom_boxplot() + 
  scale_color_manual(values = c("black", "brown","black", "black","black","black", "black","black")) +
  scale_fill_manual(name = "Group", values = c("green", "black","red","blue","pink","yellow", "purple", "gray", "brown")) +
  #scale_fill_manual (name = "Group", values = c ("green", "brown")) +
  labs(title = "") + xlab ("\nAcquisition days") + ylab("PC1\n")
boxPlots
