#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. May 2015                         ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### PCA analysis figures for frontiers paper                              ### 
#############################################################################

library("ggplot2")
library("Hmisc")
library("FactoMineR") #PCA
 
##Getting HOME directory
home <- Sys.getenv("HOME")

#path of R session is /users/cn/ierb/work/MaraDierssen/Silvina/

# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))

#install.packages("calibrate")
library(calibrate)

#old version data
# load (paste (home, "/20150515_PCA_old_frotiersPaper/data/", "ma2.RData", sep=""))
# data with only 3 new individuals
dir_data <- paste (home, "/Dropbox (CRG)/backups_old_mac", sep="")
load (paste (home, "/20150515_PCA_old_frotiersPaper/data/", "ma2_3sup_TSEGCG.RData", sep=""))
# load (paste (home, "/Dropbox (CRG)/backups_old_mac/20150515_PCA_old_frotiersPaper/data/", "ma2_3sup_TSEGCG.RData", sep=""))

####
# My own way
# Data is in ma2

# ma2_TSDoubleTT <-ma2[ma2$gentreat =="TSEEEGCG" ,]
# 
# ma2_TSDoubleTT$factorID <- as.factor(ma2_TSDoubleTT$ID) 
# 
# pca_latencies_doubleTT <- ggplot(ma2_TSDoubleTT, aes(x=day, y=latency, colour=factorID, group=factorID)) + 
#   #                           geom_path (size = 1,show_guide = T) + 
#   geom_path (size = 1,show_guide = T)  
# 
# pca_latencies_doubleTT  
##   130050203 Day 1	130050203	TS	Day 1	1239.75780	72.12244	47.8600	20.748321	24.129086	7.759941	74.07075
#47.8600  
head (ma2)
n_animals <- length(ma2[,1])
rownames(ma2) <- c(1:n_animals)
colnames(ma2)[1] <- "id"
n_col <- dim (ma2)[2]
# Saving for permutation test
# write.table(ma2, "/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/ts65_old_3sup_tsegcg_rev.csv", sep="\t")

variables_list <- colnames(ma2) [4:n_col] 

tbl_median <- with (ma2, aggregate (cbind (distance, gallindex, latency, speed, percentne, whishaw, percentperi), 
                                    list (gentreat=gentreat, day=day), FUN=median))

rownames (tbl_median) <- paste (tbl_median[,1], gsub ("Day ", "", tbl_median[,2]), sep="")

tbl_ind <- ma2

# velocity versus distance
plot(sort (tbl_ind$speed), sort(tbl_ind$distance))

tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_median_plus1 <- n_median + 1

## boxplot for thesis presentation
# boxplot_thesis <- subset(ma2, gentreat=="WT" & day=="Day 1")
# boxplot_thesis <- subset(ma2, gentreat=="WT" & day=="Day 5")
# boxplot_thesis$Value <- boxplot_thesis$latency
# 
boxplot_thesis_WT <- subset(ma2, gentreat=="WT") 
#                             & day %in% c("Day 1", "Day 5"))
# boxplot_thesis_WT$Value <- boxplot_thesis_WT$latency

boxplot_thesis_TS <- subset(ma2, gentreat=="TS") 
#                             & day %in% c("Day 1", "Day 5"))
# boxplot_thesis_TS$Value <- boxplot_thesis_TS$latency
# # boxplot_thesis <- subset(ma2, gentreat=="TS" & day=="Day 5")
# 
# boxplot_thesis_plot <- ggplot(boxplot_thesis, aes (gentreat, latency, fill = gentreat)) +
#   geom_boxplot(show.legend=FALSE, lwd=4) +
#   geom_point (position = position_jitter(width = 0.2), colour="black", show.legend=FALSE, size=8)
#   theme (axis.text=element_blank()) +
#   scale_y_continuous (limits=c(0, 60)) +
#   labs (x = "\nWT acq 5", y="Latency\n") +
#   #        labs (x = "", y="") +
#   scale_fill_manual(values=c("red")) +
#   theme(axis.title.y = element_text(size=90), axis.title.x = element_text(size=90))
# 
# boxplot_thesis_plot
boxplot_thesis_TS_5 <- subset(ma2, gentreat=="TS" & day %in% c("Day 5"))
mean(boxplot_thesis_TS_5$latency)


boxplot_thesis_plot_WT <- ggplot(boxplot_thesis_WT, aes (gentreat, latency, fill = gentreat)) +
  geom_boxplot(show.legend=FALSE, lwd=4) +
  geom_point (position = position_jitter(width = 0.2), colour="black", show.legend=FALSE, size=8) +
#   theme (axis.text=element_blank()) +
  stat_summary(fun.y=mean, colour="white", geom="point", size=10) +
  scale_y_continuous (limits=c(8, 61), breaks=c(10,60)) +
  labs (x = "\nDays", y="Latency\n") +
  #        labs (x = "", y="") +
  scale_fill_manual(values=c("red")) +
  scale_y_continuous (limits=c(8, 61), breaks=c(10,20,30,40,50,60)) +
  facet_wrap (~day, ncol = 5) +
  theme(axis.title.y = element_text(size=50), axis.title.x = element_text(size=50)) + 
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(size = 30))
boxplot_thesis_plot_WT
# ggsave (boxplot_thesis_plot_WT, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "boxplot_thesis_lat_WT.png", sep=""), width = 15, height = 10, dpi=300)  

boxplot_thesis_plot_TS <- ggplot(boxplot_thesis_TS, aes (gentreat, latency, fill = gentreat)) +
  geom_boxplot(show.legend=FALSE, lwd=4) +
#   geom_point (position = position_jitter(width = 0.2), colour="black", show.legend=FALSE, size=8) +
#   theme (axis.text=element_blank()) +
  scale_y_continuous (limits=c(8, 61), breaks=c(10,20,30,40,50,60)) +
  stat_summary(fun.y=mean, colour="red", geom="point", size=10) +
  labs (x = "\nDays", y="Latency\n") +
  #        labs (x = "", y="") +
  scale_fill_manual(values=c("limegreen")) +
  facet_wrap (~day, ncol = 5) +
  theme(axis.title.y = element_text(size=50), axis.title.x = element_text(size=50)) + 
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(size = 30))
boxplot_thesis_plot_TS

boxplot_thesis_TS

# ggsave (boxplot_thesis_plot_TS, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "boxplot_thesis_lat_TS.png", sep=""), width = 15, height = 10, dpi=300)  


# ggsave (boxplot_thesis_plot, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "boxplot_thesis_lat_acq1.png", sep=""), width = 10, height = 20, dpi=300)  
# ggsave (boxplot_thesis_plot, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "boxplot_thesis_lat_acq5.png", sep=""), width = 10, height = 20, dpi=300)  
subset(tbl_median, gentreat=="WTEE")
##################

res = PCA(tbl_med_ind[,(3:9)], scale.unit=TRUE, ind.sup=c(n_median_plus1:length(tbl_med_ind[,1]))) 

# res$var$coord[,1],res$var$coord[,2]
summary_resPCA<- summary(res)

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

pca_medians_acq <- ggplot(pca2plot, aes(x=-Dim.1, y=-Dim.2, colour=gentreat )) + 
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
pca_medians_acq

# keeping aspect ratio
pca_medians_acq_aspect_ratio <- pca_medians_acq + coord_fixed() + 
  scale_x_continuous (limits=c(-4, 4), breaks=-4:4) + 
  scale_y_continuous (limits=c(-3, 3), breaks=-3:3)

pca_medians_acq_aspect_ratio
# ggsave (pca_medians_acq_aspect_ratio, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig1_PCA/", 
#           "PCA_medians_legend.jpg", sep=""), width = 10, height = 6, dpi=900)
# ggsave (pca_medians_acq_aspect_ratio, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig1_PCA/", 
#            "PCA_medians_NO_legend.jpg", sep=""), width = 9, height = 6, dpi=900)

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
# base_size <- 12
# dailyInt_theme <- theme_update (axis.title.x = element_text (size=base_size * 2, face="bold"),
#                                 axis.title.y = element_text (size=base_size * 2, angle = 90, face="bold"),
#                                 plot.title = element_text (size=base_size * 2, face="bold"))

## Plots circle for slides thesis
size_titles <- 34
size_axis <- 28
p_circle_big_title <- p_circle_plot + coord_fixed() +
                      theme(plot.title = element_text(size=size_titles)) + 
                      theme(axis.title.x = element_text(size=size_axis)) +
                      theme(axis.title.y = element_text(size=size_axis))
p_circle_big_title
# ggsave (p_circle_big_title, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)
## Plots circle for slides thesis
################ Plots circle for slides thesis

# ggsave (p_circle_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)
# ggsave (p_circle_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig1_PCA/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)

############
## BARPLOT
df.bars <- cbind (as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE)), names(res$var$coord[,1])[order(res$var$coord[,1]^2,decreasing=TRUE)])
df.bars_to_plot <- as.data.frame(df.bars)
df.bars_to_plot$index <- as.factor (df.bars_to_plot$V2)
class (df.bars_to_plot$V1)
df.bars_to_plot$value <- as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE))
df.bars_to_plot$index <- factor(df.bars_to_plot$index, levels = df.bars_to_plot$index[order(df.bars_to_plot$value, decreasing=TRUE)])

bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
  ylim (c(0, 83)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC1\n", x = "", y="Contribution in %\n") +
  ## Plots circle for slides thesis
#   scale_x_discrete (labels=c("Gallagher","% NE","Distance","% periphery", "Whishaw","Latency","Speed")) +
  ## Plots circle for slides thesis
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot

#PLOT_paper
# ggsave (bars_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "bar_contribution.jpg", sep=""), dpi=900, height=5, width=10)

# ggsave (bars_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig4_PCA/", "bar_contribution.jpg", sep=""), dpi=900)

## SfN poster
# ggsave (bars_plot, file=paste(home, "/Dropbox (Personal)/2015_SfN/figures/", "bar_contribution.jpg", sep=""), dpi=900)

## Plots circle for slides thesis
# ggsave (bars_plot, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "bar_contribution.jpg", sep=""), dpi=900)
           
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
  ylim (c(0, 83)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC2\n", x = "", y="Contribution in %\n") +
  ## Plots circle for slides thesis
#   scale_x_discrete (labels=c("Speed", "Whishaw", "Distance", "Latency", "% periphery", "% NE", "Gallagher")) +
  ## Plots circle for slides thesis
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot_PC2

#PLOT_paper
# ggsave (bars_plot_PC2, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig4_PCA/", "bar_contribution_PC2.jpg", sep=""), dpi=900, height=5, width=10)
# Final version
# ggsave (bars_plot_PC2, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig4_PCA/", "bar_contribution_PC2.jpg", sep=""), dpi=900)

## SfN poster
# ggsave (bars_plot_PC2, file=paste(home, "/Dropbox (Personal)/2015_SfN/figures/", "bar_contribution_PC2.jpg", sep=""), dpi=900)

## Plots circle for slides thesis
# ggsave (bars_plot_PC2, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "bar_contribution_PC2.jpg", sep=""), dpi=900)

#PC3
df.bars_PC3 <- cbind (as.numeric(sort(res$var$coord[,3]^2/sum(res$var$coord[,3]^2)*100,decreasing=TRUE)), names(res$var$coord[,3])[order(res$var$coord[,3]^2,decreasing=TRUE)])
df.bars_to_plot_PC3 <- as.data.frame(df.bars_PC3)
df.bars_to_plot_PC3$index <- as.factor (df.bars_to_plot_PC3$V2)
# class (df.bars_to_plot_PC3$V1)
# df.bars_to_plot_PC2$value <- as.numeric(sort(res$var$coord[,3]^2/sum(res$var$coord[,3]^2)*100,decreasing=TRUE))
df.bars_to_plot_PC3$value <- as.numeric(sort(res$var$coord[,3]^2/sum(res$var$coord[,3]^2)*100,decreasing=TRUE))

df.bars_to_plot_PC3$index
df.bars_to_plot_PC3$index <- factor(df.bars_to_plot_PC3$index, levels = df.bars_to_plot_PC3$index[order(df.bars_to_plot_PC3$value, decreasing=TRUE)])

# df.bars_to_plot_PC2$value <- rev(df.bars_to_plot_PC2$value)
bars_plot_PC3 <- ggplot (data=df.bars_to_plot_PC3, aes(x=index, y=value)) + 
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC3\n", x = "", y="Contribution in %\n") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot_PC3


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
# ggsave (pca_plot_individuals, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_1_PCA_sup/", "PCA_individuals_review.jpg", sep=""),
#         height = 10, width = 10, dpi=900)

## Plots individuals for slides thesis
pca_plot_individuals_thesis <- ggplot (data=new_coord, aes (V1, V2)) + 
  geom_text (aes(label=day, colour = genotype), size=5, show.legend = FALSE) +
  scale_color_manual(values=c("red", "darkgreen", "blue", "lightblue", 
                              "magenta", "orange", "gray", "black")) +
#   xlim (c(-6, 6.5)) + ylim (c(-8, 6.5)) +
  xlim (c(-8, 8)) + ylim (c(-8, 6.5)) +
  #   geom_path (data=pca2plot, aes(x=-Dim.1, y=-Dim.2, colour=gentreat),size = 1,show_guide = FALSE) +
  labs(title = "Individual as supplementary points\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
       y=paste("PC2 (", var_PC2, "% of variance)\n", sep = ""))  +
  coord_fixed()

pca_plot_individuals_thesis <- pca_plot_individuals_thesis + panel_border() + theme(panel.border = element_rect(colour = "black")) +
  theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis)) +
  coord_fixed()

pca_plot_individuals_thesis

# ggsave (pca_plot_individuals_thesis, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "PCA_individuals_acq_only_points", ".png", sep=""),
#         height = 10, width = 10, dpi=300)

wt_animals <- subset (new_coord, genotype=="WT")
pca_plot_individuals_thesis_WT <- ggplot (data=wt_animals, aes (V1, V2)) + 
  geom_text (aes(label=day, colour = genotype), size=8, show.legend = FALSE) +
  scale_color_manual(values=c("red", "darkgreen", "blue", "lightblue", 
                              "magenta", "orange", "gray", "black")) +
  #   xlim (c(-6, 6.5)) + ylim (c(-8, 6.5)) +
  xlim (c(-8, 8)) + ylim (c(-8, 6.5)) +
  #   geom_path (data=pca2plot, aes(x=-Dim.1, y=-Dim.2, colour=gentreat),size = 1,show_guide = FALSE) +
  labs(title = "WT individuals as supplementary points\n", x = paste("\nPC1","", sep=""), 
       y=paste("PC2", "\n", sep = ""))  +
  coord_fixed()
# size_titles <- 50
# size_axis <- 40
pca_plot_individuals_thesis_WT <- pca_plot_individuals_thesis_WT + theme(panel.border = element_rect(colour = "black")) +
  theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis)) +
  coord_fixed()

pca_plot_individuals_thesis_WT

# ggsave (pca_plot_individuals_thesis_WT, file=paste("/Users/jespinosa/Dropbox (CRG)/thesis_presentation/figures/PCA_frontiers/", "PCA_individuals_acq_only_points_WT", ".png", sep=""),
#         height = 10, width = 10, dpi=300)

##############
# Adding the plot of PCA cloud for day A1 and day A5 only for trisomics
new_coord
# indiv <- rownames(pca_indiv_2plot)
# new_coord_indiv <- rbind(indiv, new_coord)
PC1_TS  <- subset(new_coord, grepl("TS", genotype))
# PC1_TS_acq1 <- subset(PC1_TS, day=="1", c("V1", "V2", "day","genotype", "id"))
PC1_TS_acq1 <- subset(PC1_TS, day=="1", c("V1", "V2", "day","genotype"))


colnames (PC1_TS_acq1) <- c("PC1","PC2", "day", "genotype_tt")
# Link to the explanation of how to perform density plots with ggplot2
# http://stackoverflow.com/questions/19791181/density-shadow-around-the-data-with-ggplot2-r
# Quick facts
# n controls the smoothness of the density polygon.
# h is the bandwidth of the density estimation.
# bins controls the number of density levels.
# alpha transparency allows to make the cloud more or less transparent depending on the level
# level, Computed density, is the amount you have to increase to change from one level to the next

p_cloud_ts_acq1 <- ggplot(PC1_TS_acq1, aes(PC1, PC2, color=genotype_tt)) + 
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

p_cloud_ts_acq1 <- p_cloud_ts_acq1 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_ts_acq1

PC1_TS_acq5 <- subset(PC1_TS, day=="5", c("V1", "V2", "day","genotype"))

colnames (PC1_TS_acq5) <- c("PC1","PC2", "day", "genotype_tt") 
p_cloud_ts_acq5 <- ggplot(PC1_TS_acq5, aes(PC1, PC2, color=genotype_tt)) + 
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
  labs(title = "PCA coordinates density, trisomic groups, session 5\n", x = "\nPC1", y="PC2\n")

# p_cloud_ts_acq1 + scale_alpha_continuous(range=c(0.3,0.5))
# p_cloud_ts_acq5 + scale_alpha_continuous(range=c(0.3,0.5))
p_cloud_ts_acq5 <- p_cloud_ts_acq5 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_ts_acq5

#PLOT_paper
#
setwd("/Users/jespinosa/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/")

# ggsave (p_cloud_ts_acq1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "PCA_acq1_ts_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)
# ggsave (p_cloud_ts_acq5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "PCA_acq5_ts_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)

#### BOXPLOTS of PC1
new_coord

new_coord_TS.a1 <- subset(new_coord, grepl("TS", new_coord$genotype) & day==1)
new_coord_TS.a5 <- subset(new_coord, grepl("TS", new_coord$genotype) & day==5)

boxPlots.TS.a1 <- ggplot(new_coord_TS.a1 , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
#   scale_fill_manual(name = "Genotype", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Session 1 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-5.5, 9.5)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_TS.a1[new_coord_TS.a1$genotype == "TSEEEGCG","V1"]), xend = 4.37, yend = median(new_coord_TS.a1[new_coord_TS.a1$genotype == "TSEEEGCG","V1"])), colour="white")

boxPlots.TS.a5 <- ggplot(new_coord_TS.a5 , aes (genotype, V1, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Session 5 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-5.5, 9.5)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_TS.a5[new_coord_TS.a5$genotype == "TSEEEGCG","V1"]), xend = 4.37, yend = median(new_coord_TS.a5[new_coord_TS.a5$genotype == "TSEEEGCG","V1"])), colour="white")

boxPlots.TS.a5

#PLOT_paper
# ggsave (boxPlots.TS.a1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_ts_a1.jpg", sep=""), dpi=900)
# ggsave (boxPlots.TS.a5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "boxPlot_ts_a5.jpg", sep=""), dpi=900)

### DENSITY PLOTS WT
#################
# Plots of the wt
# This plots do not change after the revision because we do not add any wt animal!!! YUPPIiiii
##############
new_coord_WT.a1 <- subset(new_coord, grepl("WT", new_coord$genotype) & day==1)
colnames (new_coord_WT.a1) <- c("PC1","PC2", "day", "genotype_tt")

p_cloud_wt_acq1 <- ggplot(new_coord_WT.a1, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 #                  geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) +
                 geom="polygon", color=NA, h=5, n=100, bins=6, show_guide = FALSE) +
  #   scale_alpha(range = c(0.00, 1)) +
  #   scale_size(range = c(0, 0.01), guide = "none")  
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values=c("red", "blue", "magenta",  "gray"),
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values=c("red", "blue", "magenta",  "gray"),
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-8, 9),
                  ylim=c(-9.5, 8)) +
  labs(title = "PCA coordinates density, wt groups, session 1\n", x = "\nPC1", y="PC2\n")

p_cloud_wt_acq1 <- p_cloud_wt_acq1 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_wt_acq1

new_coord_WT.a5 <- subset(new_coord, grepl("WT", new_coord$genotype) & day==5)
colnames (new_coord_WT.a5) <- c("PC1","PC2", "day", "genotype_tt")

p_cloud_wt_acq5 <- ggplot(new_coord_WT.a5, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=5, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("red", "blue", "magenta",  "gray"), 
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("red", "blue", "magenta",  "gray"),
                     labels = c("WT", "WTEE", "WTEGCG", "WTEEEGCG")) + 
  #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  #   coord_cartesian(xlim=c(-6, 8.5),
  #                   ylim=c(-9, 6)) +
  coord_cartesian(xlim=c(-8, 9),
                  ylim=c(-9.5, 8)) +
  labs(title = "PCA coordinates density, wt groups, session 5\n", x = "\nPC1", y="PC2\n")

# p_cloud_wt_acq1 + scale_alpha_continuous(range=c(0.3,0.5))
# p_cloud_WT_acq5 + scale_alpha_continuous(range=c(0.3,0.5))
p_cloud_wt_acq5 <- p_cloud_wt_acq5 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_wt_acq5

# Plot for paper
setwd("/Users/jespinosa/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/")

#PLOT_paper
# ggsave (p_cloud_wt_acq1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig5_PCA/", "PCA_acq1_wt_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)
# ggsave (p_cloud_wt_acq5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig5_PCA/", "PCA_acq5_wt_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)





#######################
# PC2 of WT group
# Tables with data are:
# new_coord_WT.a1
# new_coord_WT.a5

boxPlots.PC1.WT.a1 <- ggplot(new_coord_WT.a1 , aes (genotype_tt, PC1, fill = genotype_tt)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "blue", "magenta",  "gray")) +
  labs(title = "Session 1 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6), limits=c(-9, 7))

boxPlots.PC1.WT.a1 
# + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

boxPlots.PC1.WT.a5 <- ggplot(new_coord_WT.a5 , aes (genotype_tt, PC1, fill = genotype_tt)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "blue", "magenta",  "gray")) +
  labs(title = "Session 5 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6), limits=c(-9, 7)) 

#PLOT_paper
# ggsave (boxPlots.PC1.WT.a1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig5_PCA/", "boxPlot_wt_a1.jpg", sep=""), dpi=900)
# ggsave (boxPlots.PC1.WT.a5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig5_PCA/", "boxPlot_wt_a5.jpg", sep=""), dpi=900)

# Plotting a legend with scuares and colors
l <- ggplot() + geom_point(data=new_coord_WT.a1, aes (x=PC1, y=PC2, colour = genotype_tt), shape=15, size=5) +
  scale_colour_manual (values=c("red", "blue", "magenta",  "gray"))
l <- l + guides(color=guide_legend(title=NULL)) 
l <- l + theme(legend.key = element_blank())
l
# /20150515_PCA_old_frotiersPaper/figures/fig5_PCA

# ggsave (l, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig5_PCA/", "legend_squares_wt.jpg", sep=""), dpi=900)


new_coord_all.a5



# Both WT and TS PC2

new_coord_all.a1 <- subset(new_coord, day==1)
new_coord_all.a5 <- subset(new_coord, day==5)








####################
#### BOXPLOTS of PC2
new_coord

# new_coord_TS.a1 <- subset(new_coord, grepl("TS", new_coord$genotype) & day==1)
# new_coord_TS.a5 <- subset(new_coord, grepl("TS", new_coord$genotype) & day==5)

boxPlots.PC2.TS.a1 <- ggplot(new_coord_TS.a1 , aes (genotype, V2, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Session 1 PC2\n") + xlab ("\nGroups") + ylab("PC2\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6,8), limits=c(-9, 8)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_TS.a1[new_coord_TS.a1$genotype == "TSEEEGCG","V2"]), xend = 4.37, yend = median(new_coord_TS.a1[new_coord_TS.a1$genotype == "TSEEEGCG","V2"])), colour="white")

boxPlots.PC2.TS.a1

boxPlots.PC2.TS.a5 <- ggplot(new_coord_TS.a5 , aes (genotype, V2, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("darkgreen", "lightblue", "orange", "black")) +
  labs(title = "Session 5 PC2\n") + xlab ("\nGroups") + ylab("PC2\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6,8), limits=c(-9, 8)) +
  geom_segment(aes(x = 3.63, y = median(new_coord_TS.a5[new_coord_TS.a5$genotype == "TSEEEGCG","V2"]), xend = 4.37, yend = median(new_coord_TS.a5[new_coord_TS.a5$genotype == "TSEEEGCG","V2"])), colour="white")

boxPlots.PC2.TS.a5

#PLOT_paper
# ggsave (boxPlots.PC2.TS.a1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_PC2_boxplot/", "boxPlot_ts_PC2_a1.jpg", sep=""), dpi=900)
# ggsave (boxPlots.PC2.TS.a5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_PC2_boxplot/", "boxPlot_ts_PC2_a5.jpg", sep=""), dpi=900)

#######################
# PC2 of WT group
# Tables with data are:
# new_coord_WT.a1
# new_coord_WT.a5

boxPlots.PC2.WT.a1 <- ggplot(new_coord_WT.a1 , aes (genotype_tt, PC2, fill = genotype_tt)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "blue", "magenta",  "yellow")) +
  labs(title = "Session 1 PC2\n") + xlab ("\nGroups") + ylab("PC2\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6), limits=c(-9, 7))

boxPlots.PC2.WT.a1

boxPlots.PC2.WT.a5 <- ggplot(new_coord_WT.a5 , aes (genotype_tt, PC2, fill = genotype_tt)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "blue", "magenta",  "yellow")) +
  labs(title = "Session 5 PC2\n") + xlab ("\nGroups") + ylab("PC2\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6), limits=c(-9, 7)) 

boxPlots.PC2.WT.a5

#PLOT_paper
# ggsave (boxPlots.PC2.WT.a1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_PC2_boxplot/", "boxPlot_wt_PC2_a1.jpg", sep=""), dpi=900)
# ggsave (boxPlots.PC2.WT.a5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_PC2_boxplot/", "boxPlot_wt_PC2_a5.jpg", sep=""), dpi=900)

# Both WT and TS PC2

new_coord_all.a1 <- subset(new_coord, day==1)
new_coord_all.a5 <- subset(new_coord, day==5)

boxPlots.PC2.all.a1 <- ggplot(new_coord_all.a1, aes (genotype, V2, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Session 1 PC2\n") + xlab ("\nGroups") + ylab("PC2\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6,8), limits=c(-9, 9)) +
  geom_segment(aes(x = 7.63, y = median(new_coord_all.a1 [new_coord_all.a1$genotype == "TSEEEGCG","V2"]), 
                   xend = 8.37, yend = median(new_coord_all.a1 [new_coord_all.a1$genotype == "TSEEEGCG","V2"])), 
               colour="white", size =0.8)
boxPlots.PC2.all.a1

boxPlots.PC2.all.a5 <- ggplot(new_coord_all.a5, aes (genotype, V2, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "Genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Session 5 PC2\n") + xlab ("\nGroups") + ylab("PC2\n") +
  theme (legend.title=element_blank())+ 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-8,-6,-4,-2,0,2,4,6,8), limits=c(-9, 9)) +
  geom_segment(aes(x = 7.63, y = median(new_coord_all.a5 [new_coord_all.a5$genotype == "TSEEEGCG","V2"]), 
               xend = 8.37, yend = median(new_coord_all.a5 [new_coord_all.a5$genotype == "TSEEEGCG","V2"])), 
              colour="white", size =0.8)
boxPlots.PC2.all.a5

#PLOT_paper
# ggsave (boxPlots.PC2.all.a1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_PC2_boxplot/", "boxPlot_all_PC2_a1.jpg", sep=""), width=14, height=7, dpi=900)
# ggsave (boxPlots.PC2.all.a5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_PC2_boxplot/", "boxPlot_all_PC2_a5.jpg", sep=""), width=14, height=7, dpi=900)

### Plot of 3D variables 
library("plot3D")

x<- tbl_median$gallindex
y<-tbl_median$percentperi
z<-tbl_median$latency
col <- as.integer(tbl_median$gentreat)
labels_3d <- gsub("Day ", "",tbl_median$day)

plot (tbl_median$gallindex, tbl_median$speed)
plot (tbl_median$distance, tbl_median$percentperi)
plot (tbl_median$gallindex, tbl_median$whishaw)
plot (tbl_median$percentne, tbl_median$speed)

plot (tbl_median$percentperi, tbl_median$whishaw) #ok
plot (tbl_median$gallindex, tbl_median$latency) #ok
plot (tbl_median$whishaw, tbl_median$latency) #ok

?scatter3D
# scatter3D(x,y,z, bty = "g", pch = 18,
#           colvar=col, col = c("red", "darkgreen", "blue", "lightblue", 
#                               "magenta", "orange", "gray", "black"),
#           ticktype = "detailed",
#           colkey = list(at = c(1, 2, 3,4,5,6,7,8), side = 1, 
#                         addlines = TRUE, length = 0.5, width = 0.5,
#                         labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG","TSEEEGCG")),
#           xlab = "Gallagher",
#           ylab ="% periphery", zlab = "Latency")

text3D(x,y,z, bty = "g", 
       colvar=col, col = c("red", "darkgreen", "blue", "lightblue", 
                           "magenta", "orange", "gray", "black"),
       ticktype = "detailed",
       labels = labels_3d,
#        colkey = list(at = c(1, 2, 3,4,5,6,7,8), side = 1, 
#                      addlines = TRUE, length = 0.5, width = 0.5,
#                      labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG","TSEEEGCG")),
      colkey = FALSE,
       xlab = "\nGallagher",
       ylab ="\n% periphery", zlab = "\nLatency")
