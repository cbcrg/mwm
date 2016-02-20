#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Oct 2015                         ###
#############################################################################
### MWM young ts65dn mice                                                 ###
### PCA analysis                                                          ### 
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")
library("ggplot2")
library("Hmisc")
library("FactoMineR") #PCA
# install.packages("cowplot")
library("cowplot")

jtracks_data = spss.get(paste (home, "/20151001_ts65_young_MWM/data/Jtracks parameters Young TS_SUBCONJ_REV_R_FORMAT.sav", sep=""))
#smart_data = spss.get(paste (home, "/20151001_ts65_young_MWM/data/Jtracks parameters Young TS_SUBCONJ.sav", sep=""))

# data_nxf = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Ts65Dn_OLD_ACQ1_ACQ5_SUBCONJ.sav", sep=""))
# head (data_nxf)
# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))
source (paste (home, "/git/mwm/lib/R/stat_density_2d_function.R", sep=""))

# Parameter to set plot qualities
dpi_q <- 300
# dpi_q <- 900
# img_format <- ".jpg"
img_format <- ".tiff"
size_titles <- 18
size_axis <- 14

# The last records are empty
tail(jtracks_data, 50)
# tail(smart_data)

head (jtracks_data)
# head(smart_data)

# Percentage in the center is the inverse to time in the periphery, we removed it 
jtracks_data

#####################
# Last 50 rows are empty 
# Parecen que los datos son iguales elimino de la tabla el reversal 
# Me quedo con la tabla jtracks que tiene la misma estructura que la que hemos utilizado en el frontiers
jtracks_data_filt <- head(jtracks_data, -50)
tail(jtracks_data_filt)
dim (jtracks_data_filt)
young_acq <- subset(jtracks_data_filt, grepl("ACQ", DAY))
dim (young_acq)
young_acq_7var <- subset(young_acq, select=-c(GALL.DIST, PER.CENTER))
# I set same labels than in the other script (frontiers)
colnames(young_acq_7var) <- c("id", "gentreat","day", "distance", "gallindex", "latency", "speed", "percentne", "percentperi", "whishaw")
head (young_acq_7var)

# Checking that actually speed decreases from day 4 to day 5
# wt_5<- subset(young_acq_7var, day==5 & gentreat=="WT")
# wt_4<-subset(young_acq_7var, day==4 & gentreat=="WT")
# mean(wt_5$speed)
# mean(wt_4$speed)

young_acq_7var$gentreat <- gsub("H20", "", young_acq_7var$gentreat)
young_acq_7var$gentreat <- gsub("NE", "", young_acq_7var$gentreat)
young_acq_7var$day <- gsub("ACQ", "", young_acq_7var$day)
tbl4permutation <- young_acq_7var 
tbl4permutation$day <- paste("Day", young_acq_7var$day)
# Saving for permutation test
# write.table(tbl4permutation, "/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv", sep="\t")

# young_acq_7var$gentreat <- factor(young_acq_7var$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
#                                   labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
young_acq_7var$gentreat <- factor(young_acq_7var$gentreat , levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"), 
                                  labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"))

tbl_median <- with (young_acq_7var, aggregate (cbind (distance, gallindex, latency, speed, percentne, percentperi, whishaw), 
                    list (gentreat=gentreat, day=day), FUN=median))

rownames (tbl_median) <- paste (tbl_median[,1], gsub ("Day ", "", tbl_median[,2]), sep="")
tbl_ind <- young_acq_7var
tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_median_plus1 <- n_median + 1

res = PCA(tbl_med_ind[,(3:9)], scale.unit=TRUE, ind.sup=c(n_median_plus1:length(tbl_med_ind[,1]))) 

## Performing the pca of only the medians to see the directionality of variables directly from the plot of the 
## in order to avoid mistakes when I am twisting the directions to be the same as in the old paper
PCA(tbl_median[,(3:9)], scale.unit=TRUE)

tbl_med_ind[,(3:9)]
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

# pca2plot$gentreat <- factor(pca2plot$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
#                             labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
pca2plot$gentreat <- factor(pca2plot$gentreat , levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"), 
                                  labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"))

pca_medians_acq <- ggplot(pca2plot, aes(x=-Dim.1, y=-Dim.2, colour=gentreat )) + 
  geom_path (size = 1,show.legend = F) + 
  scale_color_manual(values=c("red", "darkgreen", "magenta", "black")) +
  
  #                           geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show.legend = T)+
  geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show.legend = F)+
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
  scale_x_continuous (limits=c(-4, 5), breaks=-4:5) + 
  scale_y_continuous (limits=c(-2, 3), breaks=-2:3)


pca_medians_acq_aspect_ratio_big_title <- pca_medians_acq_aspect_ratio +  theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis))

pca_medians_acq_aspect_ratio_leg <- pca_medians_acq_aspect_ratio  + geom_path (size = 1,show.legend = T)

pca_medians_acq_aspect_ratio_leg

# ggsave (pca_medians_acq_aspect_ratio_leg, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", 
#          "PCA_medians_legend.jpg", sep=""), width = 10, height = 6, dpi=900)
# ggsave (pca_medians_acq_aspect_ratio, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", 
#         "PCA_medians_NO_legend.jpg", sep=""), width = 9, height = 6, dpi=900)
# ggsave (pca_medians_acq_aspect_ratio_big_title, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", 
#                                                            "PCA_medians_NO_legend.jpg", sep=""), width = 15, height = 10, dpi=900)
        
### Circle Plot
circle_plot <- as.data.frame (res$var$coord)
labels_v <- row.names(res$var$coord)

neg_labels <- labels_v [c(1,2,3,6)]
neg_positions <- circle_plot [c(1,2,3,6), c(1,2)]

# change positions for labels
# neg_positions [2,2] <- neg_positions [2,2] - 0.03 
# neg_positions [3,2] <- neg_positions [3,2] + 0
# neg_positions [4,2] <- neg_positions [4,2] - 0.02

pos_labels <- labels_v [c(4,5,7)]
pos_positions <- circle_plot [c(4,5,7), c(1,2)]

angle <- seq(-pi, pi, length = 50)
df.circle <- data.frame(x = sin(angle), y = cos(angle))

p_circle_plot <- ggplot(circle_plot) + 
  geom_segment (data=circle_plot, aes(x=0, y=0, xend=-Dim.1, yend=-Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, color="red") +
#   xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
  scale_x_continuous(limits=c(-1.3, 1.3), breaks=(c(-1,0,1))) +
  scale_y_continuous(limits=c(-1.3, 1.3), breaks=(c(-1,0,1))) +
  geom_text (data=neg_positions, aes (x=-Dim.1, y=-Dim.2, label=neg_labels, hjust=1.2), show.legend = FALSE, size=5) + 
  geom_text (data=pos_positions, aes (x=-Dim.1, y=-Dim.2, label=pos_labels, hjust=-0.3), show.legend = FALSE, size=5) +
  geom_vline (xintercept = 0, linetype="dotted") +
  geom_hline (yintercept=0, linetype="dotted") +
  labs (title = "PCA of the variables\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)

p_circle_big_title <- p_circle_plot + coord_fixed() +
                      theme(plot.title = element_text(size=size_titles)) + 
                      theme(axis.title.x = element_text(size=size_axis)) +
                      theme(axis.title.y = element_text(size=size_axis))
# No axis
#                       theme(panel.border = element_blank(), axis.line = element_blank())
p_circle_big_title
# ggsave (p_circle_plot, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", "circle_plot.jpg", sep=""), width = 6, height = 6, dpi=900)
# ggsave (p_circle_big_title, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)

############
## BARPLOT
df.bars <- cbind (as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE)), names(res$var$coord[,1])[order(res$var$coord[,1]^2,decreasing=TRUE)])
df.bars_to_plot <- as.data.frame(df.bars)
df.bars_to_plot$index <- as.factor (df.bars_to_plot$V2)
class (df.bars_to_plot$V1)
df.bars_to_plot$value <- as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE))
df.bars_to_plot$index <- factor(df.bars_to_plot$index, levels = df.bars_to_plot$index[order(df.bars_to_plot$value, decreasing=TRUE)])

bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
  ylim (c(0, 80)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC1\n", x = "", y="Contribution in %\n") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot

bar_plot_big_title <- bars_plot + theme(plot.title = element_text(size=size_titles)) + 
                      theme(axis.title.x = element_text(size=size_axis)) +
                      theme(axis.title.y = element_text(size=size_axis))
#PLOT_paper
# ggsave (bars_plot, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq", "bar_contribution.jpg", sep=""), dpi=900)
# ggsave (bar_plot_big_title, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", "bar_contribution.jpg", sep=""), dpi=900)

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
  ylim (c(0, 80)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC2\n", x = "", y="Contribution in %\n") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot_PC2

bars_plot_PC2_big_title <- bars_plot_PC2 + theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis))

#PLOT_paper
# Final version
# ggsave (bars_plot_PC2, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", "bar_contribution_PC2.jpg", sep=""), dpi=900)
# ggsave (bars_plot_PC2_big_title, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_PCA_acq/", "bar_contribution_PC2.jpg", sep=""), dpi=900)

###################################
# Plot of supplementary individuals
day <- c()
genotype <- c()
new_coord <- as.data.frame(cbind(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2]))
dim(new_coord)

# the info of the genotype and the tt was in this table
new_coord$day <- gsub ("Day ", "", tbl_ind$day)
new_coord$genotype <- tbl_ind$gentreat

new_coord$genotype <- factor(new_coord$genotype , levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"), 
                             labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"))
  
pca_plot_individuals <- ggplot (data=new_coord, aes (V1, V2)) + 
  geom_text (aes(label=day, colour = genotype), size=5, show.legend = FALSE) +
  scale_color_manual(values=c("red", "darkgreen", "magenta", "black")) +
#   xlim (c(-5, 10)) + ylim (c(-5, 5)) +
  scale_x_continuous (limits=c(-5, 10), breaks=seq(-5, 10, by=5)) +
  scale_y_continuous (limits=c(-5, 5), breaks=seq(-5, 5, by=5)) +
  geom_path (data=pca2plot, aes(x=-Dim.1, y=-Dim.2, colour=gentreat),size = 1, show.legend = TRUE) +
  guides(color=guide_legend(guide_legend(title = "Group"))) +
  labs(title = "Individual as supplementary points\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
       y=paste("PC2 (", var_PC2, "% of variance)\n", sep = ""))  +
  coord_fixed()

pca_plot_individuals

#PLOT_paper
# ggsave (pca_plot_individuals, file=paste(home, "/20151001_ts65_young_MWM/figures/", "PCA_individuals.jpg", sep=""),
#         height = 10, width = 10, dpi=900)

head(new_coord)

## Individual variation as density plots
p_cloud_indiv_by_day <- ggplot(new_coord, aes(V1, V2, color=genotype, label=day)) + 
  stat_density2d(aes(fill=factor(genotype), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=4, bins=6, show.legend = FALSE) + 
  geom_point(show.legend = FALSE) + 
  scale_color_manual(name='genotype', 
                     values = c("red", "darkgreen", "magenta", "black"),                                            
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
#                      values = c("red","gray", "green", "black"),
#                      labels = c("WT","WTEEEGCG", "TS",  "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
#                      values = c("red","gray", "green", "black"), 
#                      labels = c("WT","WTEEEGCG", "TS",  "TSEEEGCG")) + 
                     values = c("red", "darkgreen", "magenta", "black"),                                            
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
  geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "Density plot of individual variation\n", x = "\nPC1", y="PC2\n") +
  scale_alpha_continuous(range=c(0.3,0.5)) 
p_cloud_indiv_by_day
# +
#   geom_path (data=pca2plot, aes(x=Dim.1, y=Dim.2, colour=gentreat), size = 0.5, linetype = 2, show.legend = FALSE)

p_cloud_indiv_by_day
p_cloud_indiv_by_day_facet <- p_cloud_indiv_by_day + facet_wrap(~genotype, ncol = 2)
                  
p_cloud_indiv_by_day_facet_lines <- p_cloud_indiv_by_day_facet + geom_vline(xintercept = 0, colour="gray") + 
                                    geom_hline(yintercept = 0, colour="gray")

p_cloud_indiv_by_day_facet_lines
#PLOT_presentation
# ggsave (p_cloud_indiv_by_day_facet_lines, file=paste(home, "/20151001_ts65_young_MWM/figures/", "PCA_density_indiv_byDay.jpg", sep=""), 
#          width = 10, height = 10, dpi=900)

# As I have only four groups I can plot all the grous and day 1 and 5 in the same plot
PC1_acq1_5 <- subset(new_coord, day %in% c("1", "5"), c("V1", "V2", "day","genotype"))
colnames (PC1_acq1_5) <- c("PC1","PC2", "day", "genotype_tt") 
PC1_acq1_5$day <- as.factor(PC1_acq1_5$day)

# p_cloud_acq1_5 <- ggplot(PC1_acq1_5, aes(PC1, PC2, color=genotype_tt)) + 
#                           stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
#                                          geom="polygon", color=NA, n=100, h=5, bins=6, show.legend = F) +
#                           geom_point(show.legend = F) + 
#                           scale_color_manual(name='genotype_tt', 
#                                              values = c("red", "darkgreen", "magenta", "black"),                                            
#                                              labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
#                           scale_fill_manual( name='gentreat', 
#                                              values = c("red", "darkgreen", "magenta", "black"),
#                                              labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) + 
#                           #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
#                           scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
#                           scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
#                           coord_cartesian(xlim=c(-8, 11),
#                                           ylim=c(-5, 6.5)) +
#                           labs(title = "PCA coordinates density, session 1 and 5\n", x = "\nPC1", y="PC2\n")

# p_cloud_acq1_5_facet <- p_cloud_acq1_5 + facet_grid(day ~ genotype_tt, margins=FALSE) + geom_vline(xintercept = 0, colour="gray") +
#                         geom_hline(yintercept = 0, colour="gray")
# p_cloud_acq1_5_facet

## Labels facet
levels(PC1_acq1_5$day) <- c("Session 1","Session 5")

# Trying to do density plot again, the changed the function arggggg
p_cloud_acq1_5 <- ggplot(PC1_acq1_5, aes(PC1, PC2, color=genotype_tt)) + 
#   stat_density_2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
#                   geom="polygon"), color=NA, n=100, h=4, bins=6, binwidth=6, show.legend = F) + 
  stat_density_2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=5, bins=6, show.legend = F) +  
  geom_point(show.legend = F) + 
  scale_y_continuous(limits = c(-5, 6), breaks=seq(-5,5,5)) +
  scale_x_continuous(limits = c(-6, 11), breaks=seq(-5,10,5)) +
  scale_color_manual(name='genotype_tt', 
                     values = c("red", "darkgreen", "magenta", "black"),                                            
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
  scale_fill_manual( name='gentreat', 
                     values = c("red", "darkgreen", "magenta", "black"),
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
#   scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
#   scale_y_continuous(expand=c(0.3, 0)) #+ # don't reach edges of plot.
  labs(title = "PCA coordinates density, session 1 and 5\n", x = "\nPC1", y="PC2\n")
 
# p_cloud_acq1_5
size_strips <- 12
p_cloud_acq1_5_facet <- p_cloud_acq1_5 + facet_grid(day ~ genotype_tt, margins=FALSE) + geom_vline(xintercept = 0, colour="gray") +
  geom_hline(yintercept = 0, colour="gray") +
  theme(strip.text.x = element_text(size=size_strips, face="bold"), strip.text.y = element_text(size=size_strips, face="bold", angle=90)) +
  theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis)) +
  theme(strip.background = element_blank()) + 
  panel_border() + theme(panel.border = element_rect(colour = "black"))

#   theme(panel.background = element_rect(fill=NA, col="black"))+
#   theme(panel.border = element_rect(colour = "black"))
#   theme(panel.grid.major = element_line(colour = "black"))
#   theme(strip.background = element_rect(fill="black"))
#   theme(strip.background = element_blank())

p_cloud_acq1_5_facet_coord <- p_cloud_acq1_5_facet + coord_fixed()
p_cloud_acq1_5_facet_coord

#PLOT_presentation
# ggsave (p_cloud_acq1_5_facet, file=paste(home, "/20151001_ts65_young_MWM/figures/", "PCA_density_byGroupAndDay.jpg", sep=""), 
#         width = 10, height = 6, dpi=900)

#PLOT_paper
# ggsave (p_cloud_acq1_5_facet_coord, file=paste(home, "/20151001_ts65_young_MWM/figures/", "PCA_density_byGroupAndDay_coordFixed", img_format, sep=""), 
#         dpi=dpi_q)

#### BOXPLOTS of PC1
# I can use the same table
PC1_acq1_5$genotype <- PC1_acq1_5$genotype_tt

# Changing again the labels of acquisiton
levels(PC1_acq1_5$day) <- c("1","5")
PC1.a1 <- subset(PC1_acq1_5,  day==1)
PC1.a5 <- subset(PC1_acq1_5,  day==5)

boxPlots.PC1.a1 <- ggplot(PC1.a1, aes (genotype, PC1, fill = genotype)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "black")) +
  labs(title = "Session 1 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank()) + 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8,10,12), limits=c(-4.5, 13)) +
  geom_segment(aes(x = 3.63, y = median(PC1.a1[PC1.a1$genotype == "TSEEEGCG","PC1"]), xend = 4.37, yend = median(PC1.a1[PC1.a1$genotype == "TSEEEGCG","PC1"])), colour="white")

# p_cloud_indiv_by_day_facet <- boxPlots.PC1 + facet_wrap(~genotype, ncol = 2)
boxPlots.PC1.a1.line <- boxPlots.PC1.a1 + geom_hline(yintercept = 0, colour="gray") + 
                        theme(plot.title = element_text(size=size_titles)) + 
                        theme(axis.title.x = element_text(size=size_axis)) +
                        theme(axis.title.y = element_text(size=size_axis)) +
                        panel_border() + theme(panel.border = element_rect(colour = "black"))
boxPlots.PC1.a1.line

## Session 5
boxPlots.PC1.a5 <- ggplot(PC1.a5, aes (genotype, PC1, fill = genotype)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "black")) +
  labs(title = "Session 5 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank()) + 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8,10,12), limits=c(-4.5, 13)) +
  geom_segment(aes(x = 3.63, y = median(PC1.a5[PC1.a1$genotype == "TSEEEGCG","PC1"]), xend = 4.37, yend = median(PC1.a5[PC1.a5$genotype == "TSEEEGCG","PC1"])), colour="white")

# p_cloud_indiv_by_day_facet <- boxPlots.PC1 + facet_wrap(~genotype, ncol = 2)
boxPlots.PC1.a5.line <- boxPlots.PC1.a5 + geom_hline(yintercept = 0, colour="gray") +
                        theme(plot.title = element_text(size=size_titles)) + 
                        theme(axis.title.x = element_text(size=size_axis)) +
                        theme(axis.title.y = element_text(size=size_axis)) +
                        panel_border() + theme(panel.border = element_rect(colour = "black"))


## This work but is a line not a bracket
# boxPlots.PC1.a1.line + geom_segment(aes(x = 1, y = 7.2, xend = 2, yend = 7.2), colour="black")

# I need the factor genotype, to make it work, I just add a fake genotype
sl_1 <- data.frame(x = c(1, 1, 2, 2), y = c(3.2, 3.5, 3.5, 3.2), genotype=rep("TS",4))
sl_2 <- data.frame(x = c(2, 2, 3, 3), y = c(4.4, 4.7, 4.7, 4.4), genotype=rep("TS",4))
sl_3 <- data.frame(x = c(2, 2, 4, 4), y = c(5.6, 5.9, 5.9, 5.6), genotype=rep("TS",4))

boxPlots.PC1.a1.line.sig <- boxPlots.PC1.a1.line + geom_path(data = sl_1, aes(x = x, y = y)) +
                                                   geom_path(data = sl_2, aes(x = x, y = y)) +
                                                   geom_path(data = sl_3, aes(x = x, y = y)) +
                                                   annotate("text", x=1.5, y=3.6,label="***", size=10) +
                                                   annotate("text", x=2.5, y=4.8,label="***", size=10) +
                                                   annotate("text", x=3, y=6,label="**", size=10)

boxPlots.PC1.a1.line.sig

# I need the factor genotype, to make it work, I just add a fake genotype
sl_1 <- data.frame(x = c(1, 1, 2, 2), y = c(7.2, 7.5, 7.5, 7.2), genotype=rep("TS",4))
sl_2 <- data.frame(x = c(2, 2, 3, 3), y = c(8.4, 8.7, 8.7, 8.4), genotype=rep("TS",4))
sl_3 <- data.frame(x = c(3, 3, 4, 4), y = c(9.6, 9.9, 9.9, 9.6), genotype=rep("TS",4))
sl_4 <- data.frame(x = c(2, 2, 4, 4), y = c(10.8, 11.1, 11.1, 10.8), genotype=rep("TS",4))
sl_5 <- data.frame(x = c(1, 1, 4, 4), y = c(12, 12.3, 12.3, 12), genotype=rep("TS",4))

boxPlots.PC1.a5.line.sig <- boxPlots.PC1.a5.line + geom_path(data = sl_1, aes(x = x, y = y)) +
                                                   geom_path(data = sl_2, aes(x = x, y = y)) +
                                                   geom_path(data = sl_3, aes(x = x, y = y)) +
                                                   geom_path(data = sl_4, aes(x = x, y = y)) +
                                                   geom_path(data = sl_5, aes(x = x, y = y)) +
                                                   annotate("text", x=1.5, y=7.6,label="***", size=10) +
                                                   annotate("text", x=3.5, y=10,label="*", size=10) +
                                                   annotate("text", x=2.5, y=8.8,label="***", size=10) +
                                                   annotate("text", x=3, y=11.2,label="*", size=10) +
                                                   annotate("text", x=2.5, y=12.4,label="**", size=10)

###################
###################
# boxplot but facet
####

PC1_acq1_5

## Labels facet
levels(PC1_acq1_5$day) <- c("Session 1", "Session 5")

# I need the factor genotype, to make it work, I just add a fake genotype
sl_1 <- data.frame(x = c(1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 4, 4, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 2, 2, 4, 4, 1, 1, 4, 4), 
                   y = c(3.2, 3.5, 3.5, 3.2, 4.4, 4.7, 4.7, 4.4, 5.6, 5.9, 5.9, 5.6, 7.2, 7.5, 
                         7.5, 7.2, 8.4, 8.7, 8.7, 8.4,9.6, 9.9, 9.9, 9.6, 10.8, 11.1, 11.1, 10.8, 12, 12.3, 12.3, 12), 
                   genotype=c(rep("WT",4),rep("TS",4), rep("TSEEEGCG",4),rep("WT",4),rep("TS",4), rep("TSEEEGCG",4),
                              rep("WTEEEGCG",4), rep("TS_fake",4)), 
                   day=c(rep("Session 1", 12), rep("Session 5", 20)))
sl_1$genotype <- factor(sl_1$genotype, levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG", "TS_fake"), 
                              labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG", "TS_fake"))

stars_plot <- data.frame(x_pos = c(1.5, 2.5, 3, 1.5, 3.5, 2.5, 3, 2.5), y_pos = c(3.6, 4.8, 6, 7.6, 10, 8.8, 11.2, 12.4), 
                         label=c("***","***","**", "***", "*", "***", "*", "**"), 
                         day=c(rep("Session 1", 3), 
                               rep("Session 5", 5)),
                         genotype=c(rep("WT",8)))                       

median_line_5 <- data.frame(x=3.63,y=median(PC1.a5[PC1.a5$genotype == "TSEEEGCG","PC1"]),
                           xend=4.37, yend=median(PC1.a5[PC1.a5$genotype == "TSEEEGCG","PC1"]),
                           day=factor("Session 5", levels=c("Session 1","Session 5")), genotype="TS")
median_line_1 <- data.frame(x=3.63,y= median(PC1.a1[PC1.a1$genotype == "TSEEEGCG","PC1"]),
                           xend=4.37, yend= median(PC1.a1[PC1.a1$genotype == "TSEEEGCG","PC1"]),
                           day=factor("Session 1", levels=c("Session 1","Session 5")), genotype="TS")
# PC1_acq1_5 <- factor (sl_1$genotype)

PC1_acq1_5$genotype <- factor(PC1_acq1_5$genotype, levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG", "TS_fake"), 
                              labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG", "TS_fake"))

boxPlots.PC1.facet <- ggplot(PC1_acq1_5, aes (genotype, PC1, fill = genotype)) + 
  geom_boxplot(show.legend=FALSE) +
#   scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "black", "gray")) +
  scale_fill_manual(name = "Genotype", values=c("darkgreen", "black","magenta", "red", "magenta")) +
  labs(title = "PC1 distribution\n") + xlab ("\nGroups") + ylab("PC1\n") +
  theme (legend.title=element_blank()) + 
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8,10,12), limits=c(-4.5, 13)) +
  geom_segment(data=median_line_5, aes(x = x, y = y, xend=xend, yend=yend), colour="white") +
  geom_segment(data=median_line_1, aes(x = x, y = y, xend=xend, yend=yend), colour="white") +
  facet_grid(. ~ day)

boxPlots.PC1.facet 

# p_cloud_indiv_by_day_facet <- boxPlots.PC1 + facet_wrap(~genotype, ncol = 2)
boxPlots.PC1.line <-boxPlots.PC1.facet  + geom_hline(yintercept = 0, colour="gray") +
  theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis)) +
  panel_border() + theme(panel.border = element_rect(colour = "black"))

boxPlots.PC1.line.stars <- boxPlots.PC1.line + geom_path(data = sl_1, aes(x = x, y = y)) +                  
                                               geom_text(data=stars_plot, aes(x=x_pos, y=y_pos, label=label), 
                                                         size=10, show.legend = FALSE) +
                                               theme(strip.text.x = element_text(size=size_strips, face="bold")) +
                                               theme(plot.title = element_text(size=size_titles)) + 
                                               theme(axis.title.x = element_text(size=size_axis)) +
                                               theme(axis.title.y = element_text(size=size_axis)) +
                                               theme(strip.background = element_blank()) + 
                                               panel_border() + theme(panel.border = element_rect(colour = "black"))                                                                                                      

#PLOT_paper
# ggsave (boxPlots.PC1.a1.line, file=paste(home, "/20151001_ts65_young_MWM/figures/", "boxPlot_PC1_a1", img_format, sep=""), dpi=dpi_q, units="cm", width = 10, height = 7.5)
# ggsave (boxPlots.PC1.a5.line, file=paste(home, "/20151001_ts65_young_MWM/figures/", "boxPlot_PC1_a5", img_format, sep=""), dpi=dpi_q, units="cm", width = 10, height = 7.5)

# ggsave (boxPlots.PC1.a1.line.sig, file=paste(home, "/20151001_ts65_young_MWM/figures/", "boxPlot_PC1_a1", img_format, sep=""), dpi=dpi_q)
# ggsave (boxPlots.PC1.a5.line.sig, file=paste(home, "/20151001_ts65_young_MWM/figures/", "boxPlot_PC1_a5", img_format, sep=""), dpi=dpi_q)

######
#####

img_3plots <- ggdraw() + draw_plot(p_cloud_acq1_5_facet_coord, 0, .5, 1, .5) +
              draw_label("ddd\n") +
              draw_plot(boxPlots.PC1.a1.line.sig, 0, 0, .45, .5) +
              draw_plot(boxPlots.PC1.a5.line.sig, .5, 0, .45, .5) +
              draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = size_titles)

# plot_grid (p_cloud_acq1_5_facet_coord, boxPlots.PC1.line.stars, labels=c("A","B"), nrow = 2)

img_2plots <- ggdraw() + draw_plot(p_cloud_acq1_5_facet_coord, 0, .5, 1, .5) +
              draw_plot(boxPlots.PC1.line.stars, 0, 0, 1, .5) +
              draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = size_titles)
ggsave (img_2plots, file=paste(home, "/20151001_ts65_young_MWM/figures/", "panel_boxPlot", img_format, sep=""), 
        dpi=dpi_q, width=14, height=11)
1100, 700

# Plotting a legend with scuares and colors
l <- ggplot() + geom_point(data=PC1.a5, aes (x=PC1, y=PC2, colour = genotype), shape=15, size=5) +
                scale_colour_manual (values=c("red", "darkgreen", "magenta", "black"))
l <- l + guides(color=guide_legend(title=NULL)) 
l <- l + theme(legend.key = element_blank())
l

# ggsave (l, file=paste(home, "/20151001_ts65_young_MWM/figures/", "legend_squares.jpg", sep=""), dpi=900)

# Plotting a legend with scuares and colors
l <- ggplot() + geom_lines(data=PC1.a5, aes (x=PC1, y=PC2, colour = genotype), shape=15, size=5) +
  scale_colour_manual (values=c("red", "darkgreen", "magenta", "black"))
l <- l + guides(color=guide_legend(title=NULL)) 
l <- l + theme(legend.key = element_blank())
l

## Plotting a legend with lines and colors
l <- ggplot() + geom_line(data=PC1.a5, aes (x=PC1, y=PC2, colour = genotype), shape=15, size=2) +
  scale_colour_manual (values=c("red", "darkgreen", "magenta", "black"))
l <- l + guides(color=guide_legend(title=NULL)) 
l <- l + theme(legend.key = element_blank())
l
## PLOT_paper
# ggsave (l, file=paste(home, "/20151001_ts65_young_MWM/figures/", "lines_legend.jpg", sep=""), 
#          width=14, height=7, dpi=900)
