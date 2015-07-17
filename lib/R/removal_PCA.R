#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. May 2015                         ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### Modified from Ionas script called PCA_clean for PCA analysis          ### 
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))

# Calling libraries
library(Hmisc)
library(calibrate)
  
# Loading functions:
source (paste (home, "/git/phecomp/lib/R/plotParamPublication.R", sep=""))

rem_data = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav", sep=""))
ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters except latency.sav", sep=""))
head (ma3)
# Last 5 rows are empty 
ma3 <- head(ma3,-5)

ma3_filt_rem_data <- ma3[ , grepl( "REM" , names( ma3 ) ) ]
ma3_filt_rem_data$ID <- ma3 [ , grepl( "ID" , names( ma3 ) ) ]

tail (rem_data)
tail (ma3_filt_rem_data)


# All variables in rem_data are not present in ma3_filt_rem_data
# We need to add them: NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET
# We keep IDs to perform the joining
rem_data_var = rem_data [ , c(1,7:10)]
rem_data_all_var <- merge(rem_data, ma3_filt_rem_data, all=TRUE)

# Guys starting with 1400277xx are missing in the second table, thus when merging they appear as NA, I delete them
rem_data_all_var <- head (rem_data_all_var, -5)

head (rem_data_all_var)

# Median with all available rem variables
# tbl_stat_median <-with (rem_data_all_var, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET), list (GENTREAT), FUN=function (x) median=median(x)))
tbl_stat_median <-with (rem_data_all_var, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET, DIST.REM, GALLINDEX.REM, GALLDIST.REM, SPEED.REM, PERC.NE.REM, PERC.CENTRE.REM, PERC.PERI.REM, WISHAW.REM), list (GENTREAT), FUN=function (x) median=median(x)))
genotype_tt <- as.factor(tbl_stat_median$Group.1)
genotype_tt_ionas <- as.factor(c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
tbl_stat_median <- tbl_stat_median [,-1]
PCA_rem <- prcomp(tbl_stat_median, scale=TRUE)
summary(PCA_rem)

# Plot PCA color by genotype
pca_rem_2plot <- as.data.frame (PCA_rem$x)
pca_rem_2plot$PC1_neg <- -pca_rem_2plot$PC1
g_genotype_tt <- ggplot(pca_rem_2plot, aes(PC1_neg, PC2)) + geom_point(aes(colour=genotype_tt_ionas), size=4) +                                                           
  labs (title = "PCA") +
  #   scale_color_manual(values=cols, labels=c("1", "2", "3")) +
  #   geom_text (aes (label=genotype_tt), hjust=0, vjust=-0.5)
  geom_text (aes (label=genotype_tt_ionas), hjust=0.5, vjust=-0.5)+
  xlim (c(-5, 5)) + ylim (c(-1.5,2))+
  labs (title = "PCA removal") +  
  labs (x = "\nPC1", y="PC2\n") +
  theme (legend.key=element_rect(fill=NA), legend.title=element_blank())
#, position=position_jitter(h=0), alpha = 1)

g_genotype_tt

#################################
# Individual values for boxPlot
# install.packages("FactoMineR")
library(FactoMineR)
tbl_ind_rem <- rem_data_all_var[,c(7:length(rem_data_all_var[1,]))]
tbl_ind_rem_med<- rbind (tbl_stat_median, tbl_ind_rem)

res_pca_rem = PCA(tbl_ind_rem_med, scale.unit=TRUE, ind.sup=c(9:91))

new_coord_rem <- cbind(res_pca_rem$ind.sup$coord[,1],res_pca_rem$ind.sup$coord[,2])
new_coord_rem <- as.data.frame(new_coord_rem)
new_coord_rem$genotype <- rem_data_all_var[,c(6)]
new_coord_rem$genotype <- gsub("H20", "", new_coord_rem$genotype)
new_coord_rem$genotype <- gsub("NE", "", new_coord_rem$genotype)
new_coord_rem$genotype <- factor(new_coord_rem$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                             labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_plot_individuals_rem <- ggplot (data=new_coord_rem, aes (V1, V2)) + 
#   geom_text (aes(label=day, colour = genotype), size=3, show_guide = FALSE) +
  geom_point () 

#PLOT_paper
# ggsave (pca_plot_individuals_rem, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_individuals_rem.jpg", sep=""), dpi=900)

# ANOVA
new_coord_rem_TS <- subset(new_coord_rem, grepl("TS", new_coord_rem$genotype))
new_coord_rem_TS$id <- c(1:length(new_coord_rem_TS$V1))
rem.aov <- aov(V1 ~ genotype  + Error(id), data = new_coord_rem_TS)
summary(rem.aov)

# Post-hoc
# df.anova$interaction <- paste(df.anova$group, df.anova$time, sep="_")
pairwise.t.test(new_coord_rem$V1, new_coord_rem$genotype, , p.adj="hochberg", paired=F)

# TS vs TSEEEGCG 
new_coord_rem_old_var_TS_TSEEEGCG <- subset(new_coord_rem, genotype == "TS" | genotype == "TSEEEGCG")
new_coord_rem_old_var_TS_TSEEEGCG$id <- c(1:length(new_coord_rem_old_var_TS_TSEEEGCG$V1))
new_coord_rem_old_var_TS_TSEEEGCG.aov <- aov(V1 ~ genotype  + Error(id), data = new_coord_rem_old_var_TS_TSEEEGCG)
summary(new_coord_rem_old_var_TS_TSEEEGCG.aov)


# boxplot
boxPlots_rem <- ggplot(new_coord_rem_TS , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values=c("green", "lightblue", "orange", "black")) +
  labs(title = "Removal PC1\n") + xlab ("\ngentreat") + ylab("PC1\n") +
  theme (legend.title=element_blank())

boxPlots_rem

boxPlots_rem
+ 
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-5.5, 9.5))
  
  scale_color_manual(values = c("black", "brown")) +
  scale_fill_manual(name = "Group", values = c("green", "black"), labels = c ("TS", "TSEEEGCG")) +
  #scale_fill_manual (name = "Group", values = c ("green", "brown")) +
  labs(title = "") + xlab ("\nAcquisition days") + ylab("PC1\n")

plot(res,choix="var")
cols=c("green","lightblue","black","orange","red","blue","darkgrey","magenta")




# OLD code with only data from removal data from file TS_old_removal.sav
# tbl_stat_mean <-with (df.act_sum, aggregate (cbind (V6), list (index=index, V4=V4), FUN=function (x) c (mean=mean(x), std.error=std.error(x))))
# tbl_stat_mean$mean <- tbl_stat_mean$V18 [,1]
# tbl_stat_mean$std.error <- tbl_stat_mean$V18 [,2]
tbl_stat_median <-with (rem_data, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET), list (GENTREAT), FUN=function (x) median=median(x)))

# First column are labels, transform to row.names
rownames(tbl_stat_median) <- tbl_stat_median$Group.1
tbl_stat_median <- tbl_stat_median [,-1]
res_rem <- prcomp(tbl_stat_median, scale=TRUE)
summary(res_rem)
res_rem$x
# png("figures/PCAmedians.png",res=300,width=15,height=15,unit="cm")
# dev.off()
# plot(res_rem$x[,1],res_rem$x[,2],type="n",main="PCA of group medians",xlab="PA1 (97% of variance)",ylab="PA2 (2% of variance)")
# points(res_rem$x[1,1],res_rem$x[1,2],col="red", pch=21, bg="red")
# textxy(res_rem$x[1,1],res_rem$x[1,2], row.names(res_rem$x)[1], cex=0.7,col="black")
# points(res_rem$x[2,1],res_rem$x[2,2],col="green", pch=21, bg="green")
# textxy(res_rem$x[2,1]+2,res_rem$x[2,2]+0.1, row.names(res_rem$x)[2], cex=0.7,col="black")
# points(res_rem$x[3,1],res_rem$x[3,2],col="yellow", pch=21, bg="yellow")
# textxy(res_rem$x[3,1]-0.7,res_rem$x[3,2]+0.1, row.names(res_rem$x)[3], cex=0.7,col="black")
# points(res_rem$x[4,1],res_rem$x[4,2],col="yellow", pch=21, bg="yellow")
# textxy(res_rem$x[4,1],res_rem$x[4,2], row.names(res_rem$x)[4], cex=0.7,col="black")
# points(res_rem$x[5,1],res_rem$x[5,2],col="magenta", pch=21, bg="magenta")
# textxy(res_rem$x[5,1]-0.4,res_rem$x[5,2], row.names(res_rem$x)[5], cex=0.7,col="black")
# points(res_rem$x[6,1],res_rem$x[6,2],col="magenta", pch=21, bg="lightblue")
# textxy(res_rem$x[6,1]+0.5,res_rem$x[6,2], row.names(res_rem$x)[6], cex=0.7,col="black")
# points(res_rem$x[7,1],res_rem$x[7,2],col="blue", pch=21, bg="blue")
# textxy(res_rem$x[7,1]-3.5,res_rem$x[7,2]+0.05, row.names(res_rem$x)[7], cex=0.7,col="black")
# points(res_rem$x[8,1],res_rem$x[8,2],col="gray", pch=21, bg="gray")
# textxy(res_rem$x[8,1]+1,res_rem$x[8,2]-0.05, row.names(res_rem$x)[8], cex=0.7,col="black")

# Plot PCA color by genotype
pca_rem_2plot <- as.data.frame (res_rem$x)
genotype_tt <- as.factor(row.names(res_rem$x))
length (genotype_tt)

g_genotype_tt <- ggplot(pca_rem_2plot, aes(PC1, PC2)) + geom_point(aes(colour=genotype_tt), size=4) +                                                           
  labs (title = "PCA") +
#   scale_color_manual(values=cols, labels=c("1", "2", "3")) +
#   geom_text (aes (label=genotype_tt), hjust=0, vjust=-0.5)
  geom_text (aes (label=genotype_tt), hjust=0.3, vjust=-0.5)
#, position=position_jitter(h=0), alpha = 1)
g_genotype_tt

# Biplot
PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}

PCbiplot(res_rem)


#########################
# I select for the pca the same variables than in the acquisition data of Ionas
#########################

selected_var_rem <- rem_data_all_var [,c(1:8,10,12,14,15,17,18)] 
head (selected_var_rem)
# Median with all available rem variables
# tbl_stat_median <-with (rem_data_all_var, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET), list (GENTREAT), FUN=function (x) median=median(x)))
tbl_stat_median <-with (selected_var_rem, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, LATENCY.TARGET, GALLINDEX.REM, SPEED.REM, PERC.NE.REM, PERC.PERI.REM, WISHAW.REM), list (GENTREAT), FUN=function (x) median=median(x)))                                                            
genotype_tt <- as.factor(tbl_stat_median$Group.1)
genotype_tt_ionas <- as.factor(c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
tbl_stat_median <- tbl_stat_median [,-1]
PCA_rem <- prcomp(tbl_stat_median, scale=TRUE)
summary(PCA_rem)

# Plot PCA color by genotype
pca_rem_2plot <- as.data.frame (PCA_rem$x)
pca_rem_2plot$PC1_neg <- -pca_rem_2plot$PC1
g_genotype_tt <- ggplot(pca_rem_2plot, aes(PC1_neg, PC2)) + geom_point(aes(colour=genotype_tt_ionas), size=4) +                                                           
  labs (title = "PCA") +
  #   scale_color_manual(values=cols, labels=c("1", "2", "3")) +
  #   geom_text (aes (label=genotype_tt), hjust=0, vjust=-0.5)
  geom_text (aes (label=genotype_tt_ionas), hjust=0.5, vjust=-0.5)+
  xlim (c(-5, 5)) + ylim (c(-1.5,2))+
  labs (title = "PCA removal") +  
  labs (x = "\nPC1", y="PC2\n") +
  theme (legend.key=element_rect(fill=NA), legend.title=element_blank())
#, position=position_jitter(h=0), alpha = 1)

g_genotype_tt

#################################
# Individual values for boxPlot
# install.packages("FactoMineR")
library(FactoMineR)
selected_var_rem
tbl_ind_rem <- selected_var_rem[,c(7:length(selected_var_rem[1,]))]
tbl_ind_rem_med<- rbind (tbl_stat_median, tbl_ind_rem)
res_pca_rem = PCA(tbl_ind_rem_med, scale.unit=TRUE, ind.sup=c(9:91))

new_coord_rem <- cbind(res_pca_rem$ind.sup$coord[,1],res_pca_rem$ind.sup$coord[,2])
new_coord_rem <- as.data.frame(new_coord_rem)
# length(new_coord_rem$V1)
# length(selected_var_rem$ID)

new_coord_rem$genotype <- selected_var_rem[,c(6)]
new_coord_rem$genotype <- gsub("H20", "", new_coord_rem$genotype)
new_coord_rem$genotype <- gsub("NE", "", new_coord_rem$genotype)
new_coord_rem$genotype <- factor(new_coord_rem$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                                 labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_plot_individuals_rem <- ggplot (data=new_coord_rem, aes (V1, V2)) + 
  #   geom_text (aes(label=day, colour = genotype), size=3, show_guide = FALSE) +
  geom_point () 

#PLOT_paper
# ggsave (pca_plot_individuals, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_individuals.jpg", sep=""), dpi=900)

### Circle with PCA of the variables in ggplot

row.names(res_pca_rem$var$coord)
circle_plot <- as.data.frame (res_pca_rem$var$coord)
# Changing the sign of the first dimension and second dimension so that is the same direction as in the Acq
circle_plot
circle_plot[,1] <- -circle_plot[,1]
circle_plot[,2] <- -circle_plot[,2]

labels_v <- row.names(res_pca_rem$var$coord)
labels_v
# [1] "NUMBER.ENTRIES" "PERM.TIME"      "LATENCY.TARGET" "GALLINDEX.REM"  "SPEED.REM"      "PERC.NE.REM"    "PERC.PERI.REM" 
# [8] "WISHAW.REM"    
labels_v <- c("Nentries", "permtime", "latency", "gallindex", "speed", "percentne", "percentperi", "whishaw")

# circle_plot$labels <- lab_names

pos_labels <- labels_v [c(3,4,7)]
pos_positions <- circle_plot [c(3,4,7), c(1,2)]
pos_positions[3,2] <- pos_positions[3,2] - 0.04

neg_labels <- labels_v [c(1,2,5,6,8)]
neg_positions <- circle_plot [c(1,2,5,6,8), c(1,2)]

angle <- seq(-pi, pi, length = 50)
df.circle <- data.frame(x = sin(angle), y = cos(angle))

#aes(x=PC1, y=PC2, colour=gentreat )) 

p_circle_plot <- ggplot(circle_plot) + 
  geom_segment (data=circle_plot, aes(x=0, y=0, xend=-Dim.1, yend=-Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, color="red") +
  xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
  geom_text (data=neg_positions, aes (x=-Dim.1, y=-Dim.2, label=neg_labels, hjust=1.2), show_guide = FALSE, size=5) + 
  geom_text (data=pos_positions, aes (x=-Dim.1, y=-Dim.2, label=pos_labels, hjust=-0.3), show_guide = FALSE, size=5) +
  geom_vline (xintercept = 0, linetype="dotted") +
  geom_hline (yintercept=0, linetype="dotted") +
  labs (title = "PCA of the variables\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
  #        geom_polygon(aes(x, y), data = df, inherit.aes = F, Fill=NA)
  #                         scale_x_continuous(breaks=1:10)  
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)



p_circle_plot <- ggplot(circle_plot) + 
  geom_segment (data=circle_plot, aes(x=0, y=0, xend=-Dim.1, yend=-Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, color="red") +
#   xlim (c(-1.29, 1.29)) + ylim (c(-1.4, 1.4)) +
  xlim (c(-1.3, 1.3)) + ylim (c(-1.3, 1.3)) +
  geom_text (data=neg_positions, aes (x=-Dim.1+0.25, y=-Dim.2, label=neg_labels, hjust=1.2), show_guide = FALSE, size=5) + 
  geom_text (data=pos_positions, aes (x=-Dim.1-0.3, y=-Dim.2, label=pos_labels, hjust=-0.3), show_guide = FALSE, size=5) +
  geom_vline (xintercept = 0, linetype="dotted") +
  geom_hline (yintercept=0, linetype="dotted") +
  labs (title = "PCA of the variables\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
  #        geom_polygon(aes(x, y), data = df, inherit.aes = F, Fill=NA)
  #                         scale_x_continuous(breaks=1:10)  
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)
p_circle_plot
#ggsave (p_circle_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "circle_plot_rem.jpg", sep=""), width = 10, height = 10, dpi=900)

# ANOVA
# TS all treatments
new_coord_rem_TS <- subset(new_coord_rem, grepl("TS", new_coord_rem$genotype))
new_coord_rem_TS$id <- c(1:length(new_coord_rem_TS$V1))
rem.aov <- aov(V1 ~ genotype  + Error(id), data = new_coord_rem_TS)
summary(rem.aov)

### Box plot for supplementary figure
boxPlots_rem <- ggplot(new_coord_rem_TS , aes (genotype, V1, fill = genotype)) + 
  #                    geom_boxplot() +
#   geom_boxplot(show_guide=TRUE) +
  geom_boxplot(show_guide=FALSE) +
  #             guides(color=guide_legend('Model',override.aes=list(shape=c(1,1,6,6))))
  
  scale_fill_manual(name = "Group", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Probe trial PC1\n") + xlab ("\ngentreat") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  scale_y_continuous(breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10,12), limits=c(-10.5, 12.5)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_rem_TS[new_coord_rem_TS$genotype == "TSEEEGCG","V1"]), 
                   xend = 4.37, yend = median(new_coord_rem_TS[new_coord_rem_TS$genotype == "TSEEEGCG","V1"])), colour="white")

boxPlots_rem 

#PLOT_paper
ggsave (boxPlots_rem, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_2_PCA_sup/", "boxPlot_rem.jpg", sep=""), dpi=900)
# TS vs TSEEEGCG 
new_coord_rem_TS_TSEEEGCG <- subset(new_coord_rem_TS, genotype == "TS" | genotype == "TSEEEGCG")
new_coord_rem_TS_TSEEEGCG$id <- c(1:length(new_coord_rem_TS_TSEEEGCG$V1))
rem_TS_TSEEEGCG.aov <- aov(V1 ~ genotype  + Error(id), data = new_coord_rem_TS_TSEEEGCG)
summary(rem_TS_TSEEEGCG.aov)

pairwise.t.test(new_coord_rem_TS_TSEEEGCG$V1, new_coord_rem_TS_TSEEEGCG$genotype, , p.adj="hochberg", paired=F)

# Comprobacion de las medias para ver porque no es significativo
mean (new_coord_rem_TS_TSEEEGCG [ which(new_coord_rem_TS_TSEEEGCG$genotype=="TS"),"V1"])
mean (new_coord_rem_TS_TSEEEGCG [ which(new_coord_rem_TS_TSEEEGCG$genotype=="TSEEEGCG"),"V1"])
library (plotrix)
std.error(new_coord_rem_TS_TSEEEGCG [ which(new_coord_rem_TS_TSEEEGCG$genotype=="TS"),"V1"])
std.error (new_coord_rem_TS_TSEEEGCG [ which(new_coord_rem_TS_TSEEEGCG$genotype=="TSEEEGCG"),"V1"])

stderr(new_coord_rem_TS_TSEEEGCG [ which(new_coord_rem_TS_TSEEEGCG$genotype=="TSEEEGCG"),"V1"])
boxPlots_rem <- ggplot(new_coord_rem_TS_TSEEEGCG , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values=c("green", "black")) +
  labs(title = "Removal PC1\n") + xlab ("\ngentreat") + ylab("PC1\n") +
  theme (legend.title=element_blank())

boxPlots_rem

# textxy(res_rem$x[2,1],res_rem$x[2,2], row.names(res_rem$x)[1], cex=0.7,col="green")
# lines(res_rem$x[2,1], res_rem$x[2,2],col="black")
# 
# textxy(res_rem$x[c(5,13,21,29,37),1],res_rem$x[c(5,13,21,29,37),2],c(1:5),cex=0.7,col="green")
# lines(res_rem$x[c(5,13,21,29,37),1],res_rem$x[c(5,13,21,29,37),2],col="green")
# textxy(res_rem$x[c(7,15,23,31,39),1],res_rem$x[c(7,15,23,31,39),2],c(1:5),cex=0.7,col="yellow")
# lines(res_rem$x[c(7,15,23,31,39),1],res_rem$x[c(7,15,23,31,39),2],col="yellow")
# textxy(res_rem$x[c(3,11,19,27,35),1],res_rem$x[c(3,11,19,27,35),2],c(1:5),cex=0.7,col="magenta")
# lines(res_rem$x[c(3,11,19,27,35),1],res_rem$x[c(3,11,19,27,35),2],col="magenta")
# textxy(res_rem$x[c(6,14,22,30,38),1],res_rem$x[c(6,14,22,30,38),2],c(1:5),cex=0.7,col="lightblue")
# lines(res_rem$x[c(6,14,22,30,38),1],res_rem$x[c(6,14,22,30,38),2],col="lightblue")
# textxy(res_rem$x[c(2,10,18,26,34),1],res_rem$x[c(2,10,18,26,34),2],c(1:5),cex=0.7,col="blue")
# lines(res_rem$x[c(2,10,18,26,34),1],res_rem$x[c(2,10,18,26,34),2],col="blue")
# textxy(res_rem$x[c(4,12,20,28,36),1],res_rem$x[c(4,12,20,28,36),2],c(1:5),cex=0.7,col="blue")
# lines(res_rem$x[c(4,12,20,28,36),1],res_rem$x[c(4,12,20,28,36),2],col="grey")
# legend(x="bottomleft",c("WT","TS","WTEE","TSEE","WTEGCG","TSEGCG","WTEEEGCG","TSEEEGCG"),col=c("red","green","blue","lightblue","magenta","yellow","grey","black"),lty=1,cex=0.6)
# altnamesM=c("dist","","latency\ngallindex","speed","percentne","percenter","wishaw","percentperi")
# dev.off()
# png("figures/PCAmedians_variables.png",res=300,width=15,height=15,unit="cm")
# plot(res_rem$rotation[,1],res_rem$rotation[,2],type="n",xlim=c(-1,1),asp=1,main="Variable contributions to principal axes",xlab="PA1 (80% of variance)",ylab="PA2 (11% of variance)")
# arrows(0,0,res_rem$rotation[,1],res_rem$rotation[,2],angle=5,length=0.1)
# textxy(res_rem$rotation[,1],res_rem$rotation[,2],altnamesM,cex=1)
# dev.off()
# 
# 
# 
# 
# 
# ResMMALL=prcomp(MMall,scale=TRUE)
# 
# MM5[1,]=apply(M[treat=="WT",],2,median)
# 
# M=as.matrix(ma2[rem_data$day=="Day 5",c(4,5,7,8,9,10,11,12)])
# treat=as.matrix(ma2[ma2$day=="Day 5",2])
# M4=as.matrix(ma2[ma2$day=="Day 4",c(4,5,7,8,9,10,11,12)])
# treat4=as.matrix(ma2[ma2$day=="Day 4",2])
# M3=as.matrix(ma2[ma2$day=="Day 3",c(4,5,7,8,9,10,11,12)])
# treat3=as.matrix(ma2[ma2$day=="Day 3",2])
# M2=as.matrix(ma2[ma2$day=="Day 2",c(4,5,7,8,9,10,11,12)])
# treat2=as.matrix(ma2[ma2$day=="Day 2",2])
# M1=as.matrix(ma2[ma2$day=="Day 1",c(4,5,7,8,9,10,11,12)])
# treat1=as.matrix(ma2[ma2$day=="Day 1",2])
# 
# 
# 
