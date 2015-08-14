# Analysis of the PAT data
# Lectura de datos

library(Hmisc)
library(doBy)
library(nlme)
library(multcomp)

#setwd("C:/Users/scatuara/Desktop/Statistics Klaus/PAT/SPSS")
setwd("/Users/jespinosa/MWM/20150811_graficosPosterSilvina/")

datos <- spss.get("PA hasta G14 grupos hembras_LongFormat.sav", lowernames=T)

# no original data in mail, I load the Rdata there the table is the correct one
load ("/Users/jespinosa/MWM/20150811_graficosPosterSilvina/PAdata_May29.RData")

datos <- orderBy(~id+day, datos)
rownames(datos) <- 1:nrow(datos)
dim(datos)
View(datos)
summary(datos)
# Diferencias resp. valor basal
datos$pa.diff <- with(datos, pa.shock - rep(pa.shock[seq(1, nrow(datos), 3)], each=3))

# Modelos para analizar los datos
# 8 groups, 44 contrasts
#################################
# modeloPAdiff <- lme(pa.diff~gentreatment*as.factor(day), datos, random=~1|id)
# summary(modeloPAdiff )
# anova(modeloPAdiff )

# Matriz de contrastes para modeloPAdiff
#lev <- levels(datos$gentreatment)
#ctmat <- matrix(0, nr = 44, nc = 24)
#colnames(ctmat) <- c('Icept', lev[-1], paste('Time', 2:3), paste(rep(lev[-1], 2), 'T', rep(2:3, each=7), sep=''))
#rownames(ctmat)[44] <- 'Day 6: WTEEEGCG vs. TSEEEGCG'
#rownames(ctmat)[1:7] <- paste("Day 1: WTNEH20", lev[-1], sep=' vs. ')
#rownames(ctmat)[8:11] <- paste("Day 1: TSNEH20", lev[c(4,5,6,8)], sep=' vs. ')
#rownames(ctmat)[12:14] <- paste("Day 1: WTEEH20", lev[c(4,5,7)], sep=' vs. ')
#rownames(ctmat)[15:16] <- paste("Day 1: TSEEH20", lev[c(6,8)], sep=' vs. ')
#rownames(ctmat)[17:19] <- paste("Day 1: WTNEEGCG", lev[6:8], sep=' vs. ')
#rownames(ctmat)[20:21] <- paste("Day 1: TSNEEGCG", lev[7:8], sep=' vs. ')
#rownames(ctmat)[22] <- 'Day 1: WTEEEGCG vs. TSEEEGCG'
#rownames(ctmat)[23:29] <- paste("Day 6: WTNEH20", lev[-1], sep=' vs. ')
#rownames(ctmat)[30:33] <- paste("Day 6: TSNEH20", lev[c(4,5,6,8)], sep=' vs. ')
#rownames(ctmat)[34:36] <- paste("Day 6: WTEEH20", lev[c(4,5,7)], sep=' vs. ')
#rownames(ctmat)[37:38] <- paste("Day 6: TSEEH20", lev[c(6,8)], sep=' vs. ')
#rownames(ctmat)[39:41] <- paste("Day 6: WTNEEGCG", lev[6:8], sep=' vs. ')
#rownames(ctmat)[42:43] <- paste("Day 6: TSNEEGCG", lev[7:8], sep=' vs. ')
## A continuaci?n se han introducido los valores a mano: fix(ctmat)
#comment(ctmat) <- 'Matrix for 44 contrasts'

#comps <- glht(modeloPAdiff, linfct = ctmat)
#summary(comps)
#confint(comps)
#cons <- round(cbind(confint(comps)$confint, summary(comps)$test$pva), 4)
#colnames(cons)[4] <- 'p value'
#cons

## Less contrasts
#ctmat.red <- ctmat[c(3, 20, 23, 28, 31, 33, 36, 40, 42, 44), ]
#comps <- glht(update(modeloPAdiff, log(pa.diff)~.), linfct = ctmat.red)
#summary(comps)
#confint(comps)
#cons <- round(cbind(confint(comps)$confint, summary(comps)$test$pva), 4)
#colnames(cons)[4] <- 'p value'
#cons

## Figures of all 8 groups
##########################
mean1 <- with(subset(datos, day==1), tapply(pa.diff, gentreatment, mean, na.rm=T))
mean6 <- with(subset(datos, day==6), tapply(pa.diff, gentreatment, mean, na.rm=T))

install.packages("beeswarm")


library(beeswarm)
jpeg('BoxplotsD1.jpg', width=800, height=400)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0))
boxplot(pa.diff~gentreatment, subset(datos, day==1), xlab='Treatment group', main='After 1 day', pch=16, ylim=c(-100, 350), xaxt='n')
axis(1, 1:8, labels=c('Name1',....., 'Name8'))
beeswarm(pa.diff~gentreatment, subset(datos, day==1), add=T, pch=16)
title('Differences from baseline', outer=T, cex.main=1.7)
points(1:8, mean1, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
dev.off()

####################
# Code figure here
####################

jpeg('BoxplotsD6.jpg', width=800, height=400)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0))
boxplot(pa.diff~gentreatment, subset(datos, day==6), xlab='Treatment group', main='After 6 days', pch=16, ylim=c(-100, 350), xaxt='n')
axis(1, 1:8, labels=c('Name1',....., 'Name8'))
beeswarm(pa.diff~gentreatment, subset(datos, day==6), add=T, pch=16)
#points(1:8, mean6, type='b', pch=17, col=2, lwd=2)
#legend('topleft', 'Mean', lwd=2, pch=17, col=2)
title('Differences from baseline', outer=T, cex.main=1.7)
dev.off()

####################
# Code Jose here
####################
head (rem_data_all_var)
head (datos)

# Solo coge el dia 6 por eso me sale diferente
datos_d6 <- datos [which(datos$day == 6),]
datos_d6$gentreatment <- gsub("H20", "", datos_d6$gentreatment)
datos_d6$gentreatment <- gsub("NE", "", datos_d6$gentreatment)
datos_d6$gentreatment <- factor(datos_d6$gentreatment , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                                    labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

# boxplot pa.diff
bp_pa.diff <- ggplot(datos_d6 , aes (gentreatment, pa.diff, fill = gentreatment)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Title here\n") + xlab ("\nGroups") + ylab("Your axis title here\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(-102, 352), breaks = seq(-100,350, by=50)) +
  
  geom_segment(aes(x = 7.63, y = median(datos_d6 [datos_d6$gentreatment == "TSEEEGCG", "pa.diff"]), 
                   xend = 8.37, yend = median(datos_d6 [datos_d6$gentreatment == "TSEEEGCG","pa.diff"])), 
                   colour="white", size =0.8)

bp_pa.diff
bp_pa.diff_points <- bp_pa.diff + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_pa.diff_points

ggsave (bp_pa.diff_points, file="pa_diff_d6.jpg", width=14, height=7, dpi=900)
# 4 groups, 2*6 contrasts
#########################
#> levels(datos$gentreatment)
#[1] "WTNEH20"  "TSNEH20"  "WTEEH20"  "TSEEH20"  "WTNEEGCG" "TSNEEGCG" "WTEEEGCG" "TSEEEGCG"
subdats <- subset(datos, substr(gentreatment, 1, 1)=='T')
subdats$genotype <- NULL
subdats$gentreatment <- factor(subdats$gentreatment)

submodPAdiff <- lme(pa.diff~gentreatment*as.factor(day), subdats, random=~1|id)
summary(submodPAdiff)
anova(submodPAdiff)

# Pairwise comparisons on both days
#lev <- levels(subdats$gentreatment)
#ctmat12 <- matrix(0, nr = 12, nc = 12)
#colnames(ctmat12) <- c('Icept', lev[-1], paste('Time', 2:3), paste(rep(lev[-1], 2), 'T', rep(2:3, each=3), sep=''))
#rownames(ctmat12)[12] <- 'Day 6: TSNEEGCG vs. TSEEEGCG'
#rownames(ctmat12)[1:3] <- paste("Day 1: TSNEH20", lev[-1], sep=' vs. ')
#rownames(ctmat12)[4:5] <- paste("Day 1: TSEEH20", lev[c(3, 4)], sep=' vs. ')
#rownames(ctmat12)[6] <- 'Day 1: TSNEEGCG vs. TSEEEGCG'
#rownames(ctmat12)[7:9] <- paste("Day 6: TSNEH20", lev[-1], sep=' vs. ')
#rownames(ctmat12)[10:11] <- paste("Day 6: TSEEH20", lev[c(3, 4)], sep=' vs. ')
#rm(lev)
## A continuaci?n se han introducido los valores a mano: fix(ctmat12)
#comment(ctmat12) <- 'Matrix for 12 contrasts'

comps <- glht(submodPAdiff, linfct = ctmat12)
summary(comps)
confint(comps)
cons <- round(cbind(confint(comps)$confint, summary(comps)$test$pva), 4)
colnames(cons)[4] <- 'p value'
cons

## Figure
#########
mean1 <- with(subset(subdats, day==1), tapply(pa.diff, gentreatment, mean, na.rm=T))
mean6 <- with(subset(subdats, day==6), tapply(pa.diff, gentreatment, mean, na.rm=T))

jpeg('Boxplots4groups.jpg', width=800, height=800)
par(mfrow=c(2, 1), font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0))
boxplot(pa.diff~gentreatment, subset(subdats, day==1), xlab='Treatment group', main='After 1 day', pch=16, ylim=c(-100, 350), xaxt='n')
axis(1, at=1:4, c('Name1', 'Name2', 'Name3', 'name4'))
beeswarm(pa.diff~gentreatment, subset(subdats, day==1), add=T, pch=16)
points(1:4, mean1, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
boxplot(pa.diff~gentreatment, subset(subdats, day==6), xlab='Treatment group', main='After 6 days', pch=16, ylim=c(-100, 350), xaxt='n')
axis(1, at=1:4, c('Name1', 'Name2', 'Name3', 'name4'))
beeswarm(pa.diff~gentreatment, subset(subdats, day==6), add=T, pch=16)
points(1:4, mean6, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
title('Differences from baseline', outer=T, cex.main=2)
dev.off()

# 4 groups (now wild types), 2*6 contrasts
##########################################
#> levels(datos$gentreatment)
#[1] "WTNEH20"  "TSNEH20"  "WTEEH20"  "TSEEH20"  "WTNEEGCG" "TSNEEGCG" "WTEEEGCG" "TSEEEGCG"
subdatsW <- subset(datos, substr(gentreatment, 1, 1)=='W')
subdatsW$genotype <- NULL
subdatsW$gentreatment <- factor(subdatsW$gentreatment)

submodPAdiffW <- lme(pa.diff~gentreatment*as.factor(day), subdatsW, random=~1|id)
summary(submodPAdiffW)
anova(submodPAdiffW)

# Pairwise comparisons on both days
lev <- levels(subdatsW$gentreatment)
ctmat12W <- matrix(0, nr = 12, nc = 12)
colnames(ctmat12W) <- c('Icept', lev[-1], paste('Time', 2:3), paste(rep(lev[-1], 2), 'T', rep(2:3, each=3), sep=''))
rownames(ctmat12W)[12] <- 'Day 6: WTNEEGCG vs. WTEEEGCG'
rownames(ctmat12W)[1:3] <- paste("Day 1: WTNEH20", lev[-1], sep=' vs. ')
rownames(ctmat12W)[4:5] <- paste("Day 1: WTEEH20", lev[c(3, 4)], sep=' vs. ')
rownames(ctmat12W)[6] <- 'Day 1: WTNEEGCG vs. WTEEEGCG'
rownames(ctmat12W)[7:9] <- paste("Day 6: WTNEH20", lev[-1], sep=' vs. ')
rownames(ctmat12W)[10:11] <- paste("Day 6: WTEEH20", lev[c(3, 4)], sep=' vs. ')
rm(lev)
ctmat12W[1:12, 1:12] <- ctmat12[1:12, 1:12]
comment(ctmat12W) <- 'Matrix for 12 contrasts (Wild type)'

compsW <- glht(submodPAdiffW, linfct = ctmat12W)
summary(compsW)
confint(compsW)
consW <- round(cbind(confint(compsW)$confint, summary(compsW)$test$pva), 4)
colnames(consW)[4] <- 'p value'
consW

## Figures
##########
mean1 <- with(subset(subdatsW, day==1), tapply(pa.diff, gentreatment, mean, na.rm=T))
mean6 <- with(subset(subdatsW, day==6), tapply(pa.diff, gentreatment, mean, na.rm=T))

jpeg('Boxplots4groupsD1W.jpg', width=800, height=900)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0))
boxplot(pa.diff~gentreatment, subset(subdatsW, day==1), xlab='Treatment group', main='After 1 day', pch=16, ylim=c(-75, 350), xaxt='n')
axis(1, at=1:4, labels=c('H1', 'h2', 'h3', 'h4')) 
axis(2, at=seq(-50, 350, 100))
beeswarm(pa.diff~gentreatment, subset(subdatsW, day==1), add=T, pch=16)
points(1:4, mean1, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
dev.off()
jpeg('Boxplots4groupsD6W.jpg', width=800, height=900)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0))
boxplot(pa.diff~gentreatment, subset(subdatsW, day==6), xlab='Treatment group', main='After 6 days', pch=16, ylim=c(-75, 350), xaxt='n')
axis(1, at=1:4, labels=c('H1', 'h2', 'h3', 'h4')) 
axis(2, at=seq(-50, 350, 100))
beeswarm(pa.diff~gentreatment, subset(subdatsW, day==6), add=T, pch=16)
points(1:4, mean6, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
title('Differences from baseline', outer=T, cex.main=2)
dev.off()

### Extra comparisons
#####################

subsd1 <- subset(datos, day==1 & as.numeric(gentreatment) <=2)
t.test(pa.diff~gentreatment, subsd1, var.equal=T)
subsd6 <- subset(datos, day==6 & as.numeric(gentreatment) <=2)
t.test(pa.diff~gentreatment, subsd6, var.equal=T)
rm(subsd1, subsd6)


############################
subdats <- subset(datos, gentreatment %in% levels(gentreatment)[c(1, 2, 7, 8)])
subdats$gentreatment <- factor(subdats$gentreatment)
mean0 <- with(subset(subdats, day==0), tapply(pa.shock, gentreatment, mean, na.rm=T))
mean1 <- with(subset(subdats, day==1), tapply(pa.diff, gentreatment, mean, na.rm=T))
mean6 <- with(subset(subdats, day==6), tapply(pa.diff, gentreatment, mean, na.rm=T))

jpeg('Boxplots4groupsBaseline_NEW.jpg', width=800, height=900)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0), cex.lab=1.5, cex.axis=1.4)
boxplot(pa.shock~gentreatment, subset(subdats, day==0), xlab='Treatment group', main='Baseline', pch=16, ylim=c(0, 350), cex.main=2, ylab='Seconds')
beeswarm(pa.shock~gentreatment, subset(subdats, day==0), add=T, pch=16, col=1:4)
points(1:4, mean0, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
dev.off()

jpeg('Boxplots4groupsDay1OrigData.jpg', width=800, height=900)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0), cex.lab=1.5, cex.axis=1.4)
boxplot(pa.shock~gentreatment, subset(subdats, day==1), xlab='Treatment group', main='Day 1', pch=16, ylim=c(0, 350), cex.main=2, ylab='Seconds')
beeswarm(pa.shock~gentreatment, subset(subdats, day==1), add=T, pch=16, col=1:4)
points(1:4, mean1, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
dev.off()

jpeg('Boxplots4groupsDay6OrigData.jpg', width=800, height=900)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0), cex.lab=1.5, cex.axis=1.4)
boxplot(pa.shock~gentreatment, subset(subdats, day==6), xlab='Treatment group', main='Day 6', pch=16, ylim=c(0, 350), cex.main=2, ylab='Seconds')
beeswarm(pa.shock~gentreatment, subset(subdats, day==6), add=T, pch=16, col=1:4)
points(1:4, mean6, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
dev.off()

jpeg('Boxplots4groupsD1_NEW.jpg', width=800, height=900)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0))
boxplot(pa.diff~gentreatment, subset(subdats, day==1), xlab='Treatment group', main='After 1 day', pch=16, ylim=c(-75, 350))
axis(2, at=seq(-50, 350, 100))
beeswarm(pa.diff~gentreatment, subset(subdats, day==1), add=T, pch=16)
points(1:4, mean1, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
dev.off()
jpeg('Boxplots4groupsD6_NEW.jpg', width=800, height=900)
par(font=2, font.lab=2, font.axis=2, oma=c(0, 0, 1, 0))
boxplot(pa.diff~gentreatment, subset(subdats, day==6), xlab='Treatment group', main='After 6 days', pch=16, ylim=c(-75, 350))
axis(2, at=seq(-50, 350, 100))
beeswarm(pa.diff~gentreatment, subset(subdats, day==6), add=T, pch=16)
points(1:4, mean6, type='b', pch=17, col=2, lwd=2)
legend('topleft', 'Mean', lwd=2, pch=17, col=2)
title('Differences from baseline', outer=T, cex.main=2)
dev.off()
