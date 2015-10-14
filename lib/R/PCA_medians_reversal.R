#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Oct 2015                         ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### Reversal session PCA                                                  ### 
#############################################################################

# Calling libraries
library(Hmisc)
library(calibrate)
library(multcomp)
library(ggplot2)

##Getting HOME directory
home <- Sys.getenv("HOME") 

# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))
rev_data <- spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters_all_including_3_extra_mice_reestrustured.sav", sep=""))
colnames(rev_data)[1] <- "id"
rev_data <- subset(rev_data, grepl("Rev", day))
head(rev_data,20)

# Latencies are in a different table
latencies.df=spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_OLD_LATENCIES_A_MANO_ALL_3_mice.sav", sep=""))
head(latencies.df)
colnames(latencies.df)[1] <- "id"

# I only need REVERSAL data
LAT.REV1 <- subset(latencies.df, select=c(id, LAT.REV1))
LAT.REV2 <- subset(latencies.df, select=c(id, LAT.REV2))
LAT.REV3 <- subset(latencies.df, select=c(id, LAT.REV3))
colnames(LAT.REV1)[2]<-"latency"
colnames(LAT.REV2)[2]<-"latency"
colnames(LAT.REV3)[2]<-"latency"
LAT.REV1$day <- "Rev 1"
LAT.REV2$day <- "Rev 2"
LAT.REV3$day <- "Rev 3"

latencies.rev<-rbind (LAT.REV1, LAT.REV2, LAT.REV3)

# I generate a table for each reversal session with all the other variables, the simplest way to mix everything
reversal_var_lat <- merge(rev_data, latencies.rev, by=c("id","day"), all = TRUE)
head (reversal_var_lat)

#### WARNING animal 130054726 has any value for removal regardless which table I am using
# I delete it
reversal_var_lat_f <- subset (reversal_var_lat, id != 130054726) 
head (reversal_var_lat_f)
n_animals <- length(reversal_var_lat_f[,1])
rownames(reversal_var_lat_f) <- c(1:n_animals)

# I remove the variable percentage center because is exactly equal to 1-percentperi
rev_data_f_6v <- subset(reversal_var_lat_f, select=-c(percenter))
n_col <- dim (rev_data_f_6v)[2]

# Setting the same name for the variables such as in the other scripts
colnames (rev_data_f_6v) [4] <- "distance"
colnames (rev_data_f_6v) [3] <- "gentreat"
colnames (rev_data_f_6v) [7] <- "percentsw"

#Removing Rev from days colunm
rev_data_f_6v$day <- gsub ("Rev ", "", rev_data_f_6v$day)
head (rev_data_f_6v)
# Saving for permutation test
# write.table(rev_data_f_6v, "/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/rev_data_f_6v.csv", sep="\t")

variables_list <- colnames(rev_data_f_6v) [4:n_col] 

tbl_median <- with (rev_data_f_6v, aggregate (cbind (distance, gallindex, speed, percentsw, perperi, wishaw, latency), 
                                    list (gentreat=gentreat, day=day), FUN=median))

rownames (tbl_median) <- paste (tbl_median[,1], gsub ("Day ", "", tbl_median[,2]), sep="")

tbl_ind <- rev_data_f_6v
tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_median_plus1 <- n_median + 1
max_var_i <- dim(tbl_med_ind)[2]
# velocity vs density
# plot(sort (tbl_med_ind$speed), sort(tbl_med_ind$distance))

res = PCA(tbl_med_ind[,(3:max_var_i)], scale.unit=TRUE, ind.sup=c(n_median_plus1:length(tbl_med_ind[,1]))) 

# summary_resPCA<- summary(res)

# Variance of PC1 and PC2
var_PC1 <- round (res$eig [1,2])
var_PC2 <- round (res$eig [2,2])

# Coordinates are store here
# res$ind$coord --- rownames(res$ind$coord)
pca2plot <- as.data.frame (res$ind$coord[1:n_median,])
pca2plot$gen_day <- row.names(pca2plot)
pca2plot$days <-  as.factor(as.numeric (gsub(".*([0-9]+)$", "\\1", pca2plot$gen_day)))
pca2plot$gentreat <-  as.factor(gsub("([A-Z]+).*$", "\\1", pca2plot$gen_day))

pca2plot$gentreat <- factor(pca2plot$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                            labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_medians_rev <- ggplot(pca2plot, aes(x=-Dim.1, y=-Dim.2, colour=gentreat )) + 
  #                           geom_path (size = 1,show_guide = T) + 
  geom_path (size = 1,show_guide = F) + 
  scale_color_manual(values=c("red", "darkgreen", "blue", "lightblue", 
                              "magenta", "orange", "gray", "black")) +
  #                           geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show_guide = T)+
  geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show_guide = F)+
  theme(legend.key=element_rect(fill=NA)) +
  labs(title = "PCA of group medians\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
       y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  #                           guides(colour = guide_legend(override.aes = list(size = 10)))+
  guides(colour = guide_legend(override.aes = list(size = 1)))+
  theme(legend.key=element_rect(fill=NA))

#PLOT_paper
pca_medians_rev

# keeping aspect ratio
pca_medians_rev_aspect_ratio <- pca_medians_rev + coord_fixed() + 
  scale_x_continuous (limits=c(-4, 4), breaks=-4:4) + 
  scale_y_continuous (limits=c(-3, 3), breaks=-3:3)

pca_medians_rev_aspect_ratio
# ggsave (pca_medians_rev_aspect_ratio, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", 
#           "PCA_medians_reversal_legend.jpg", sep=""), width = 10, height = 6, dpi=900)
# ggsave (pca_medians_rev_aspect_ratio, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", 
#            "PCA_medians_reversal_NO_legend.jpg", sep=""), width = 9, height = 6, dpi=900)

### Circle Plot
circle_plot <- as.data.frame (res$var$coord)
labels_v <- row.names(res$var$coord)

neg_labels <- labels_v [c(1,2,3,7)]
neg_positions <- circle_plot [c(1,2,3,7), c(1,2)]

# change positions for labels
# neg_positions [2,2] <- neg_positions [2,2] - 0.03 
# neg_positions [3,2] <- neg_positions [3,2] + 0
neg_positions [4,2] <- neg_positions [4,2] - 0.02

pos_labels <- labels_v [c(4,5,6)]
pos_positions <- circle_plot [c(4,5,6), c(1,2)]

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
  labs (title = "PCA of the variables\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  #        geom_polygon(aes(x, y), data = df, inherit.aes = F, Fill=NA)
  #                         scale_x_continuous(breaks=1:10)  
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)
base_size <- 12
dailyInt_theme <- theme_update (axis.title.x = element_text (size=base_size * 2, face="bold"),
                                axis.title.y = element_text (size=base_size * 2, angle = 90, face="bold"),
                                plot.title = element_text (size=base_size * 2, face="bold"))

p_circle_plot

# ggsave (p_circle_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig6_Reversal_PCA/", "circle_plot_reversal.jpg", sep=""), width = 10, height = 10, dpi=900)

############
## BARPLOT
df.bars <- cbind (as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE)), names(res$var$coord[,1])[order(res$var$coord[,1]^2,decreasing=TRUE)])
df.bars_to_plot <- as.data.frame(df.bars)
df.bars_to_plot$index <- as.factor (df.bars_to_plot$V2)
class (df.bars_to_plot$V1)
df.bars_to_plot$value <- as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE))
df.bars_to_plot$index <- factor(df.bars_to_plot$index, levels = df.bars_to_plot$index[order(df.bars_to_plot$value, decreasing=TRUE)])

bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
  ylim (c(0, 71)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC1\n", x = "", y="Contribution in %\n") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot

#PLOT_paper
# ggsave (bars_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "bar_contribution.jpg", sep=""), dpi=900)


df.bars_PC2 <- cbind (as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE)), names(res$var$coord[,2])[order(res$var$coord[,2]^2,decreasing=TRUE)])
df.bars_to_plot_PC2 <- as.data.frame(df.bars_PC2)
df.bars_to_plot_PC2$index <- as.factor (df.bars_to_plot_PC2$V2)
# class (df.bars_to_plot_PC2$V1)
# df.bars_to_plot_PC2$value <- as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE))
df.bars_to_plot_PC2$value <- as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE))

df.bars_to_plot_PC2$index
df.bars_to_plot_PC2$index <- factor(df.bars_to_plot_PC2$index, levels = df.bars_to_plot_PC2$index[order(df.bars_to_plot_PC2$value, decreasing=TRUE)])

# df.bars_to_plot_PC2$value <- rev(df.bars_to_plot_PC2$value)
bars_plot_PC2 <- ggplot (data=df.bars_to_plot_PC2, aes(x=index, y=value)) + 
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC2\n", x = "", y="Contribution in %\n") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot_PC2

#PLOT_paper
# Final version
# ggsave (bars_plot_PC2, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "bar_contribution_PC2.jpg", sep=""), dpi=900)










###################################
# Plot of supplementary individuals
day <- c()
genotype <- c()
new_coord <- as.data.frame(cbind(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2]))

# the info of the genotype and the tt was in this table
#tbl_ind
new_coord$day <- gsub ("Day ", "", tbl_ind$day)
new_coord$genotype <- tbl_ind$gentreat

new_coord$genotype <- factor(new_coord$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                             labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_plot_individuals <- ggplot (data=new_coord, aes (V1, V2)) + 
  geom_text (aes(label=day, colour = genotype), size=5, show_guide = FALSE) +
  scale_color_manual(values=c("red", "darkgreen", "blue", "lightblue", 
                              "magenta", "orange", "gray", "black")) +
  xlim (c(-6, 6.5)) + ylim (c(-8, 6.5)) +
  geom_path (data=pca2plot, aes(x=-Dim.1, y=-Dim.2, colour=gentreat),size = 1,show_guide = FALSE) +
  labs(title = "Individual as supplementary points\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
       y=paste("PC2 (", var_PC2, "% of variance)\n", sep = ""))  +
  coord_fixed()

pca_plot_individuals

#PLOT_paper
# ggsave (pca_plot_individuals, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "PCA_individuals_review.jpg", sep=""),
#         height = 10, width = 10, dpi=900)


##############
# Adding the plot of PCA cloud for day A1 and day A5 only for trisomics
new_coord
# indiv <- rownames(pca_indiv_2plot)
# new_coord_indiv <- rbind(indiv, new_coord)
PC1_TS  <- subset(new_coord, grepl("TS", genotype))
# PC1_TS_rev1 <- subset(PC1_TS, day=="1", c("V1", "V2", "day","genotype", "id"))
PC1_TS_rev1 <- subset(PC1_TS, day=="1", c("V1", "V2", "day","genotype"))


colnames (PC1_TS_rev1) <- c("PC1","PC2", "day", "genotype_tt")
# Link to the explanation of how to perform density plots with ggplot2
# http://stackoverflow.com/questions/19791181/density-shadow-around-the-data-with-ggplot2-r
# Quick facts
# n controls the smoothness of the density polygon.
# h is the bandwidth of the density estimation.
# bins controls the number of density levels.
# alpha transparency allows to make the cloud more or less transparent depending on the level
# level, Computed density, is the amount you have to increase to change from one level to the next

p_cloud_ts_rev1 <- ggplot(PC1_TS_rev1, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 #                  geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) +
                 geom="polygon", color=NA, h=5, n=100, bins=6, show_guide = FALSE) +
  #   scale_alpha(range = c(0.00, 1)) +
  #   scale_size(range = c(0, 0.01), guide = "none")  
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values=c("darkgreen", "lightblue", "darkorange", "black"), 
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values=c("darkgreen", "lightblue", "orange", "black"),
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
  #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  #   coord_cartesian(xlim=c(-7, 9),
  #                   ylim=c(-10, 10)) +
  coord_cartesian(xlim=c(-8, 9),
                  ylim=c(-10.5, 8)) +
  labs(title = "PCA coordinates density, trisomic groups, session 1\n", x = "\nPC1", y="PC2\n")

p_cloud_ts_rev1 <- p_cloud_ts_rev1 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_ts_rev1

PC1_TS_rev3 <- subset(PC1_TS, day=="3", c("V1", "V2", "day","genotype"))

colnames (PC1_TS_rev3) <- c("PC1","PC2", "day", "genotype_tt") 
p_cloud_ts_rev3 <- ggplot(PC1_TS_rev3, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=5, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("darkgreen", "lightblue", "darkorange", "black"), 
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("darkgreen", "lightblue","orange", "black"),
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
  #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  #   coord_cartesian(xlim=c(-6, 8.5),
  #                   ylim=c(-9, 6)) +
  coord_cartesian(xlim=c(-8, 9),
                  ylim=c(-10.5, 8)) +
  labs(title = "PCA coordinates density, trisomic groups, session 3\n", x = "\nPC1", y="PC2\n")

# p_cloud_ts_rev1 + scale_alpha_continuous(range=c(0.3,0.5))
# p_cloud_ts_rev3 + scale_alpha_continuous(range=c(0.3,0.5))
p_cloud_ts_rev3 <- p_cloud_ts_rev3 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_ts_rev3

#PLOT_paper
#
setwd("/Users/jespinosa/20150515_PCA_old_frotiersPaper/figures/fig6_Reversal_PCA/")

# ggsave (p_cloud_ts_rev1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig6_Reversal_PCA/", "PCA_rev1_ts_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)
# ggsave (p_cloud_ts_rev3, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig6_Reversal_PCA/", "PCA_rev3_ts_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)

#### BOXPLOTS of PC1
new_coord

new_coord_TS.rev1 <- subset(new_coord, grepl("TS", new_coord$genotype) & day==1)
new_coord_TS.rev3 <- subset(new_coord, grepl("TS", new_coord$genotype) & day==3)

boxPlots.TS.rev1 <- ggplot(new_coord_TS.rev1 , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Reversal session 1 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4), limits=c(-6, 4.5)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_TS.rev1[new_coord_TS.rev1$genotype == "TSEEEGCG","V1"]), xend = 4.37, yend = median(new_coord_TS.rev1[new_coord_TS.rev1$genotype == "TSEEEGCG","V1"])), colour="white")

boxPlots.TS.rev1

boxPlots.TS.rev3 <- ggplot(new_coord_TS.rev3 , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Reversal session 3 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4), limits=c(-6, 5)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_TS.rev3[new_coord_TS.rev3$genotype == "TSEEEGCG","V1"]), xend = 4.37, yend = median(new_coord_TS.rev3[new_coord_TS.rev3$genotype == "TSEEEGCG","V1"])), colour="white")

boxPlots.TS.rev3

boxPlots.TS.rev1.line <- boxPlots.TS.rev1 + geom_hline(yintercept = 0, colour="gray")
boxPlots.TS.rev3.line <- boxPlots.TS.rev3 + geom_hline(yintercept = 0, colour="gray")

#PLOT_paper
ggsave (boxPlots.TS.rev1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_ts_rev1.jpg", sep=""),width = 10, height = 7, dpi=900)
ggsave (boxPlots.TS.rev3, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_ts_rev3.jpg", sep=""),width = 10, height = 7, dpi=900)

# ggsave (boxPlots.TS.rev1.line, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_ts_rev1.jpg", sep=""), dpi=900)
# ggsave (boxPlots.TS.rev3.line, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_ts_rev3.jpg", sep=""), dpi=900)

### DENSITY PLOTS WT
#################
# Plots of the wt
##############
new_coord_WT.rev1 <- subset(new_coord, grepl("WT", new_coord$genotype) & day==1)
colnames (new_coord_WT.rev1) <- c("PC1","PC2", "day", "genotype_tt")

p_cloud_wt_rev1 <- ggplot(new_coord_WT.rev1, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 #                  geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) +
                 geom="polygon", color=NA, h=5, n=100, bins=6, show_guide = FALSE) +
  #   scale_alpha(range = c(0.00, 1)) +
  #   scale_size(range = c(0, 0.01), guide = "none")  
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values=c("red", "blue", "magenta",  "yellow"),
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values=c("red", "blue", "magenta",  "yellow"),
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-8, 9),
                  ylim=c(-9.5, 8)) +
  labs(title = "PCA coordinates density, wt groups, session 1\n", x = "\nPC1", y="PC2\n")

p_cloud_wt_rev1 <- p_cloud_wt_rev1 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_wt_rev1

ggsave (p_cloud_wt_rev1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "p_cloud_wt_rev1.jpg", sep=""), dpi=900)

new_coord_WT.rev3 <- subset(new_coord, grepl("WT", new_coord$genotype) & day==3)
colnames (new_coord_WT.rev3) <- c("PC1","PC2", "day", "genotype_tt")

p_cloud_wt_rev3 <- ggplot(new_coord_WT.rev3, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=5, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("red", "blue", "magenta",  "yellow"), 
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("red", "blue", "magenta",  "yellow"),
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  #   coord_cartesian(xlim=c(-6, 8.5),
  #                   ylim=c(-9, 6)) +
  coord_cartesian(xlim=c(-8, 9),
                  ylim=c(-9.5, 8)) +
  labs(title = "PCA coordinates density, wt groups, session 3\n", x = "\nPC1", y="PC2\n")

# p_cloud_wt_rev1 + scale_alpha_continuous(range=c(0.3,0.5))
# p_cloud_WT_rev3 + scale_alpha_continuous(range=c(0.3,0.5))
p_cloud_wt_rev3 <- p_cloud_wt_rev3 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_wt_rev3

ggsave (p_cloud_wt_rev3, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "p_cloud_wt_rev3.jpg", sep=""), dpi=900)
                                                 





#### BOXPLOTS of PC1

new_coord_WT.rev1 <- subset(new_coord, grepl("WT", new_coord$genotype) & day==1)
new_coord_WT.rev3 <- subset(new_coord, grepl("WT", new_coord$genotype) & day==3)

boxPlots.WT.rev1 <- ggplot(new_coord_WT.rev1 , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "blue", "magenta",  "yellow")) + 
  labs(title = "Reversal session 1 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6), limits=c(-6, 6)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_TS.rev1[new_coord_TS.rev1$genotype == "TSEEEGCG","V1"]), xend = 4.37, yend = median(new_coord_TS.rev1[new_coord_TS.rev1$genotype == "TSEEEGCG","V1"])), colour="white")

boxPlots.WT.rev1

boxPlots.WT.rev3 <- ggplot(new_coord_WT.rev3 , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "blue", "magenta",  "yellow")) +
  labs(title = "Reversal session 3 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6), limits=c(-6, 6)) 
# +
#   geom_segment(aes(x = 3.63, y = median(new_coord_WT.rev3[new_coord_TS.rev3$genotype == "TSEEEGCG","V1"]), xend = 4.37, yend = median(new_coord_TS.rev3[new_coord_TS.rev3$genotype == "TSEEEGCG","V1"])), colour="white")

boxPlots.WT.rev3

# boxPlots.WT.rev1.line <- boxPlots.TS.rev1 + geom_hline(yintercept = 0, colour="gray")
# boxPlots.WT.rev3.line <- boxPlots.TS.rev3 + geom_hline(yintercept = 0, colour="gray")

#PLOT_paper
# ggsave (boxPlots.WT.rev1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_wt_rev1.jpg", sep=""),width = 10, height = 7, dpi=900)
# ggsave (boxPlots.WT.rev3, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_wt_rev3.jpg", sep=""),width = 10, height = 7, dpi=900)

# ggsave (boxPlots.TS.rev1.line, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_ts_rev1.jpg", sep=""), dpi=900)
# ggsave (boxPlots.TS.rev3.line, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/figReversal_PCA/", "boxPlot_ts_rev3.jpg", sep=""), dpi=900)
















