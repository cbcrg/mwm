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
library("cowplot")

size_titles <- 18
size_axis <- 14
size_strips <- 12
dpi_q <- 300 
dpi_q <- 300
# dpi_q <- 900
# img_format <- ".jpg"
img_format <- ".tiff"

#############
### Latency old
############

jtracks_data = spss.get(paste (home, "/20151001_ts65_young_MWM/data/Jtracks parameters Young TS_SUBCONJ_REV_R_FORMAT.sav", sep=""))
smart_data = spss.get(paste (home, "/20151001_ts65_young_MWM/data/Jtracks parameters Young TS_SUBCONJ.sav", sep=""))
head (smart_data)

# head (data_nxf)
# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))
source (paste (home, "/git/mwm/lib/R/stat_density_2d_function.R", sep=""))

# The last records are empty
tail(jtracks_data, 50)
# tail(smart_data)

head (jtracks_data)
# head(smart_data)

# Percentage in the center is the inverse to time in the periphery, we removed it 
jtracks_data

#####################
# Last 50 rows are empty 
# Parecen que los datos son iguales elimino de la tabla el acquisition
# Me quedo con la tabla jtracks que tiene la misma estructura que la que hemos utilizado en el frontiers
jtracks_data_filt <- head(jtracks_data, -50)
tail(jtracks_data_filt)
dim (jtracks_data_filt)
young_rev <- subset(jtracks_data_filt, grepl("REV", DAY))
dim (young_rev)
head (young_rev)

# I remove the variable percentage center because is exactly equal to 1-percentperi
# We use gall index and not gall distance
young_rev_7var <- subset(young_rev, select=-c(GALL.DIST, PER.CENTER))
head(young_rev_7var)

# I set same labels than in the other script (frontiers)
colnames(young_rev_7var) <- c("id", "gentreat","day", "distance", "gallindex", "latency", "speed", "percentsw", "percentperi", "whishaw")
head (young_rev_7var)

young_rev_7var$gentreat <- gsub("H20", "", young_rev_7var$gentreat)
young_rev_7var$gentreat <- gsub("NE", "", young_rev_7var$gentreat)
young_rev_7var$day <- gsub("REV", "", young_rev_7var$day)

# Removing outlier
young_rev_7var_no_130019287 <- subset (young_rev_7var, !id == "130019287") 
young_rev_7var <- young_rev_7var_no_130019287
tbl4permutation <- young_rev_7var 
tbl4permutation$day <- paste("Day", young_rev_7var$day)

# Saving for permutation test
# write.table(tbl4permutation, "/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young_rev_no_130019287.csv", sep="\t")

# young_rev_7var$gentreat <- factor(young_acq_7var$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
#                                   labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
young_rev_7var$gentreat <- factor(young_rev_7var$gentreat , levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"), 
                                  labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"))

tbl_median <- with (young_rev_7var, aggregate (cbind (distance, gallindex, latency, speed, percentsw, percentperi, whishaw), 
                                               list (gentreat=gentreat, day=day), FUN=median))

rownames (tbl_median) <- paste (tbl_median[,1], gsub ("Day ", "", tbl_median[,2]), sep="")
tbl_ind <- young_rev_7var
tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_median_plus1 <- n_median + 1
pca_medians_rev_aspect_ratio
res = PCA(tbl_med_ind[,(3:9)], scale.unit=TRUE, ind.sup=c(n_median_plus1:length(tbl_med_ind[,1]))) 
pca_medians_rev_aspect_ratio
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

## To have the same direction than the old mice reversal I change the direction on only the X axis
pca_medians_rev <- ggplot(pca2plot, aes(x=-Dim.1, y=Dim.2, colour=gentreat )) + 
  #                           geom_path (size = 1,show.legend = T) + 
                          geom_path (size = 1, show.legend = F) + 
#                           scale_color_manual (values=c("red", "darkgreen", "magenta", "black")) +
                          scale_color_manual (values=c("red", "darkgreen", "magenta", "green")) +                          
                          #                           geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show.legend = T)+
                          geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show.legend = F)+
                          theme (legend.key=element_rect(fill=NA)) +
                          labs (title = "PCA of group medians\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
                                y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
                          #                           guides(colour = guide_legend(override.aes = list(size = 10)))+
                          guides(colour = guide_legend(override.aes = list(size = 1)))+
                          theme(legend.key=element_rect(fill=NA))

#PLOT_paper
pca_medians_rev

# keeping aspect ratio
pca_medians_rev_aspect_ratio <- pca_medians_rev + coord_fixed() + 
                                                  scale_x_continuous (limits=c(-5, 5), breaks=-5:5) + 
                                                  scale_y_continuous (limits=c(-2.5, 2.5), breaks=-2:2)

pca_medians_rev_aspect_ratio <- pca_medians_rev_aspect_ratio + 
                                panel_border() + theme(panel.border = element_rect(colour = "black")) + 
                                theme(plot.title = element_text(size=size_titles)) + 
                                theme(axis.title.x = element_text(size=size_axis)) +
                                theme(axis.title.y = element_text(size=size_axis)) + coord_fixed()
pca_medians_rev_aspect_ratio_leg <- pca_medians_rev_aspect_ratio + geom_path (size = 1, show.legend = T) +
                                    theme(legend.title=element_blank()) +                                     
                                    theme(legend.text = element_text(size = 11)) +
#                                     theme(legend.position = c(0.87, 0.87))
                                    theme(legend.position = c(0.15, 0.2))
pca_medians_rev_aspect_ratio_leg
# ggsave (pca_medians_rev_aspect_ratio, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", 
#         "PCA_medians_legend_rev.jpg", sep=""), width = 10, height = 6, dpi=900)
# ggsave (pca_medians_rev_aspect_ratio, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", 
#         "PCA_medians_NO_legend_rev.jpg", sep=""), width = 9, height = 6, dpi=900)

### Circle Plot
circle_plot <- as.data.frame (res$var$coord)
labels_v <- row.names(res$var$coord)

neg_labels <- labels_v [c(1,2,3,6)]
neg_positions <- circle_plot [c(1,2,3,6), c(1,2)]

# change positions for labels
neg_positions [1,2] <- neg_positions [1,2] + 0.05 

pos_labels <- labels_v [c(4,5,7)]
pos_positions <- circle_plot [c(4,5,7), c(1,2)]

angle <- seq(-pi, pi, length = 50)
df.circle <- data.frame(x = sin(angle), y = cos(angle))

## To have the same direction than the old mice reversal I change the direction on only the X axis  
p_circle_plot <- ggplot(circle_plot) + 
  geom_segment (data=circle_plot, aes(x=0, y=0, xend=-Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, color="red") +
  scale_x_continuous(limits=c(-1.8, 1.8), breaks=(c(-1,0,1))) +
  scale_y_continuous(limits=c(-1.8, 1.8), breaks=(c(-1,0,1))) +
  geom_text (data=neg_positions, aes (x=-Dim.1, y=Dim.2, label=neg_labels, hjust=1.2), show.legend = FALSE, size=5) + 
  geom_text (data=pos_positions, aes (x=-Dim.1, y=Dim.2, label=pos_labels, hjust=-0.3), show.legend = FALSE, size=5) +
  geom_vline (xintercept = 0, linetype="dotted") +
  geom_hline (yintercept=0, linetype="dotted") +
  labs (title = "PCA of the variables\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)

p_circle_plot <- p_circle_plot + panel_border() + theme(panel.border = element_rect(colour = "black")) + 
                 theme(plot.title = element_text(size=size_titles)) + 
                 theme(axis.title.x = element_text(size=size_axis)) +
                 theme(axis.title.y = element_text(size=size_axis)) + coord_fixed()
p_circle_plot

# ggsave (p_circle_plot, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)
# ggsave (p_circle_plot, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)

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
bars_plot <- bars_plot + panel_border() + theme(panel.border = element_rect(colour = "black")) +
             ylim (c(0, 80)) +
             theme(plot.title = element_text(size=size_titles)) + 
             theme(axis.title.x = element_text(size=size_axis)) +
             theme(axis.title.y = element_text(size=size_axis))

#PLOT_paper
# ggsave (bars_plot, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "bar_contribution.jpg", sep=""), dpi=900)

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
bars_plot_PC2 <- bars_plot_PC2 + panel_border() + theme(panel.border = element_rect(colour = "black"))

bars_plot_PC2 <- bars_plot_PC2 + panel_border() + theme(panel.border = element_rect(colour = "black")) +
                 ylim (c(0, 80)) +
                 theme(plot.title = element_text(size=size_titles)) + 
                 theme(axis.title.x = element_text(size=size_axis)) +
                 theme(axis.title.y = element_text(size=size_axis))

bars_plot_PC2

#PLOT_paper
# Final version
# ggsave (bars_plot_PC2, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "bar_contribution_PC2.jpg", sep=""), dpi=900)

#################
#######
## Panel
panel_pca <- ggdraw() + draw_plot (pca_medians_rev_aspect_ratio_leg, 0, .5, 0.5, .5) +
  draw_plot (p_circle_plot, 0.5, 0.5, 0.5, 0.5) +
  draw_plot (bars_plot, 0, 0, 0.5, .5) +
  draw_plot (bars_plot_PC2, 0.5, 0, 0.5, .5) +
#   draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0), c(1, 1, 0.5), size = size_titles)
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = size_titles)
panel_pca

img_format=".tiff"
ggsave (panel_pca, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "panel_PCA_reversal", img_format, sep=""), 
        dpi=dpi_q, width=15, height=11)

###################################
# Plot of supplementary individuals
day <- c()
genotype <- c()

# I only change direction of X axis 
new_coord <- as.data.frame(cbind(-res$ind.sup$coord[,1],res$ind.sup$coord[,2]))
dim(new_coord)

# the info of the genotype and the tt was in this table
new_coord$day <- gsub ("Day ", "", tbl_ind$day)
new_coord$genotype <- tbl_ind$gentreat

new_coord$genotype <- factor(new_coord$genotype , levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"), 
                             labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"))

pca_plot_individuals <- ggplot (data=new_coord, aes (V1, V2)) + 
  geom_text (aes(label=day, colour = genotype), size=5, show.legend = FALSE) +
#   scale_color_manual(values=c("red", "darkgreen", "magenta", "black")) +
  scale_color_manual(values=c("red", "darkgreen", "magenta", "green")) +
  xlim (c(-6, 6)) + ylim (c(-6, 6)) +
  geom_path (data=pca2plot, aes(x=-Dim.1, y=Dim.2, colour=gentreat),size = 1,show.legend = TRUE) +
  labs(title = "Individual as supplementary points, reversal\n", x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
       y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  theme (legend.title=element_blank()) + 
  coord_fixed()

pca_plot_individuals <- pca_plot_individuals + panel_border() + theme(panel.border = element_rect(colour = "black"))

pca_plot_individuals

#PLOT_paper
panel_PCA_reversal
# ggsave (pca_plot_individuals, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "PCA_individuals_rev.tiff", sep=""),
#         height = 5, width = 5, dpi=300)

head(new_coord)
## Individual variation as density plots
p_cloud_indiv_by_day <- ggplot(new_coord, aes(V1, V2, color=genotype, label=day)) + 
  stat_density2d(aes(fill=factor(genotype), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=4, bins=6, show.legend = FALSE) + 
  geom_point(show.legend = FALSE) + 
  scale_color_manual(name='genotype', 
#                      values = c("red", "darkgreen", "magenta", "black"), 
                     values = c("red", "darkgreen", "magenta", "green"), 
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
  #                      values = c("red","gray", "green", "black"),
  #                      labels = c("WT","WTEEEGCG", "TS",  "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     #                      values = c("red","gray", "green", "black"), 
                     #                      labels = c("WT","WTEEEGCG", "TS",  "TSEEEGCG")) + 
#                      values = c("red", "darkgreen", "magenta", "black"),                                            
                     values = c("red", "darkgreen", "magenta", "green"),                   
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
  geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "Density plot of individual variation\n", x = "\nPC1", y="PC2\n") +
  scale_alpha_continuous(range=c(0.3,0.5)) 
# +
#   geom_path (data=pca2plot, aes(x=Dim.1, y=Dim.2, colour=gentreat), size = 0.5, linetype = 2, show.legend = FALSE)

p_cloud_indiv_by_day
p_cloud_indiv_by_day_facet <- p_cloud_indiv_by_day + facet_wrap(~genotype, ncol = 2)

p_cloud_indiv_by_day_facet_lines <- p_cloud_indiv_by_day_facet + geom_vline(xintercept = 0, colour="gray") + 
  geom_hline(yintercept = 0, colour="gray")

#PLOT_presentation
# ggsave (p_cloud_indiv_by_day_facet_lines, file=paste(home, "/20151001_ts65_young_MWM/figures/", "PCA_density_indiv_byDay.jpg", sep=""), 
#          width = 10, height = 10, dpi=900)

# As I have only four groups I can plot all the grous and day 1 and 5 in the same plot
PC1_rev1_3 <- subset(new_coord, day %in% c("1", "3"), c("V1", "V2", "day","genotype"))
colnames (PC1_rev1_3) <- c("PC1","PC2", "day", "genotype_tt") 
PC1_rev1_3$day <- as.factor(PC1_rev1_3$day)

## Labels facet
levels(PC1_rev1_3$day) <- c("Session 1","Session 3")

p_cloud_rev1_3 <- ggplot(PC1_rev1_3, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=5, bins=6, show.legend = F) +
  geom_point(show.legend = F) + 
  scale_y_continuous(limits = c(-5, 8), breaks=seq(-5,8,5)) +
  scale_x_continuous(limits = c(-15, 7.8), breaks=seq(-15,8,5)) +
  scale_color_manual(name='genotype_tt', 
#                      values = c("red", "darkgreen", "magenta", "black"),   
                     values = c("red", "darkgreen", "magenta", "green"),   
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) +
  scale_fill_manual( name='gentreat', 
#                      values = c("red", "darkgreen", "magenta", "black"),
                     values = c("red", "darkgreen", "magenta", "green"),
                     labels = c("WT", "TS", "WTEEEGCG", "TSEEEGCG")) + 
  #   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
#   scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
#   scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
#   coord_cartesian(xlim=c(-8, 11),
#                   ylim=c(-5, 6.5)) +
  labs(title = "PCA coordinates density, reversal session 1 and 3\n", x = "\nPC1", y="PC2\n")

p_cloud_rev1_3_facet <- p_cloud_rev1_3 + facet_grid(day ~ genotype_tt, margins=FALSE) + geom_vline(xintercept = 0, colour="gray") +
  geom_hline(yintercept = 0, colour="gray")
p_cloud_rev1_3_facet <- p_cloud_rev1_3_facet + panel_border() + theme(panel.border = element_rect(colour = "black")) +
                        theme(plot.title = element_text(size=size_titles)) + 
                        theme(axis.title.x = element_text(size=size_axis)) +
                        theme(axis.title.y = element_text(size=size_axis)) 
                        
p_cloud_rev1_3_facet_strips <-  p_cloud_rev1_3_facet + theme(strip.text.x = element_text(size=size_strips, face="bold"), strip.text.y = element_text(size=size_strips, face="bold", angle=90)) +
                                theme(plot.title = element_text(size=size_titles)) + 
                                theme(axis.title.x = element_text(size=size_axis)) +
                                theme(axis.title.y = element_text(size=size_axis)) +
                                theme(strip.background = element_blank()) +
                                panel_border() + theme(panel.border = element_rect(colour = "black"))
p_cloud_rev1_3_facet_strips_coord <- p_cloud_rev1_3_facet_strips + coord_fixed()

#PLOT_presentation
# ggsave (p_cloud_rev1_3_facet, file=paste(home, "/20151001_ts65_young_MWM/figures/", "PCA_density_byGroupAndDay.jpg", sep=""), 
#         width = 10, height = 6, dpi=900)

#### BOXPLOTS of PC1
# I can use the same table
PC1_rev1_3$genotype <- PC1_rev1_3$genotype_tt

# PC1.rev1 <- subset(PC1_rev1_3,  day==1)
# PC1.rev3 <- subset(PC1_rev1_3,  day==3)
PC1.rev1 <- subset(PC1_rev1_3,  day=="Session 1")
PC1.rev3 <- subset(PC1_rev1_3,  day=="Session 3")
box_plot_lab = "\n"

boxPlots.PC1.rev1 <- ggplot(PC1.rev1, aes (genotype, PC1, fill = genotype)) + 
  geom_boxplot(show.legend=FALSE) +
#   scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "black")) +
  scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "green")) +
#   labs(title = "PC1 distribution\n") + xlab ("\nGroups") + ylab("PC1\n") +
#   theme (legend.title=element_blank()) + 
  labs(title = "\n") + xlab (box_plot_lab) + ylab("PC1\n") +
  theme (legend.title=element_blank(), axis.text=element_blank()) +    
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
#   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-5.5, 9.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6), limits=c(-6, 6)) #+
#   geom_segment(aes(x = 3.63, y = median(PC1.rev1[PC1.rev1$genotype == "TSEEEGCG","PC1"]), xend = 4.37, yend = median(PC1.rev1[PC1.rev1$genotype == "TSEEEGCG","PC1"])), colour="white")

# p_cloud_indiv_by_day_facet <- boxPlots.PC1 + facet_wrap(~genotype, ncol = 2)
boxPlots.PC1.rev1.line <- boxPlots.PC1.rev1 + geom_hline(yintercept = 0, colour="gray")
boxPlots.PC1.rev1.line

## Session 3
boxPlots.PC1.rev3 <- ggplot(PC1.rev3, aes (genotype, PC1, fill = genotype)) + 
  geom_boxplot(show.legend=FALSE) +
#   scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "black")) +
  scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "green")) +
#   labs(title = "Session 3 PC1\n") + xlab ("\nGroups") + ylab("PC1\n") +
#   theme (legend.title=element_blank()) + 
  labs(title = "\n") + xlab (box_plot_lab) + ylab("PC1\n") +
  theme (legend.title=element_blank(), axis.text=element_blank()) +   
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6), limits=c(-6, 6)) #+
#   geom_segment(aes(x = 3.63, y = median(PC1.rev3[PC1.rev3$genotype == "TSEEEGCG","PC1"]), xend = 4.37, yend = median(PC1.rev3[PC1.rev3$genotype == "TSEEEGCG","PC1"])), colour="white")

# p_cloud_indiv_by_day_facet <- boxPlots.PC1 + facet_wrap(~genotype, ncol = 2)
boxPlots.PC1.rev3.line <- boxPlots.PC1.rev3 + geom_hline(yintercept = 0, colour="gray")
boxPlots.PC1.rev3.line
max(PC1.rev3$PC1)

#PLOT_paper
# ggsave (boxPlots.PC1.rev1.line, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "boxPlot_PC1_rev1.jpg", sep=""), dpi=900)
# ggsave (boxPlots.PC1.rev3.line, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "boxPlot_PC1_rev3.jpg", sep=""), dpi=900)

################
#### BOXPLOTS facet
PC1_rev1_3$genotype <- PC1_rev1_3$genotype_tt

median_line_r1 <- data.frame(x=3.63,y=median(PC1.rev1[PC1.rev1$genotype == "TSEEEGCG","PC1"]),
                             xend=4.37, yend = median(PC1.rev1[PC1.rev1$genotype == "TSEEEGCG","PC1"]),
                             day=factor("Session 1", levels=c("Session 1","Session 3")), genotype="TS")

median_line_r3 <- data.frame(x=3.63,y= median(PC1.rev3[PC1.rev3$genotype == "TSEEEGCG","PC1"]),
                             xend=4.37, yend= median(PC1.rev3[PC1.rev3$genotype == "TSEEEGCG","PC1"]),
                             day=factor("Session 3", levels=c("Session 1","Session 3")), genotype="TS")

PC1_rev1_3$genotype <- factor(PC1_rev1_3$genotype_tt, levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG", "TS_fake"), 
                              labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG", "TS_fake"))

boxPlots.PC1.rev.facet <- ggplot(PC1_rev1_3, aes (genotype, PC1, fill = genotype)) + 
  geom_boxplot(show.legend=FALSE) +
#   scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "black")) +
  scale_fill_manual(name = "Genotype", values=c("red", "darkgreen", "magenta", "green")) +
#   scale_fill_manual(name = "Genotype", values=c("darkgreen", "black","magenta", "red", "magenta")) +
#   labs(title = "PC1 distribution\n") + xlab ("\nGroups") + ylab("PC1\n") +
#   theme (legend.title=element_blank()) + 
  labs(title = "PC1 distribution\n") + xlab (box_plot_lab) + ylab("PC1\n") +
  theme (legend.title=element_blank(), axis.text=element_blank()) +     
  # Same axis limits in day 1 and day 5
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-6, 0.5)) +
  #   scale_y_continuous(breaks=c(-4,-2,0,2,4,6,8), limits=c(-5.5, 9.5)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6,8,10,12), limits=c(-6, 12)) +
#   geom_segment(data=median_line_r3, aes(x = x, y = y, xend=xend, yend=yend), colour="white") +
#   geom_segment(data=median_line_r1, aes(x = x, y = y, xend=xend, yend=yend), colour="white") +
  facet_grid(. ~ day) +
  geom_hline(yintercept = 0, colour="gray")

# p_cloud_indiv_by_day_facet <- boxPlots.PC1 + facet_wrap(~genotype, ncol = 2)
boxPlots.PC1.rev3.line <- boxPlots.PC1.rev.facet + geom_hline(yintercept = 0, colour="gray")
boxPlots.PC1.rev3.line

sl_1 <- data.frame(x = c(1, 1, 2, 2, 
#                          2, 2, 3, 3, 
                         2, 2, 4, 4,
#                          3, 3, 4, 4,
                         1, 1, 4, 4,
                         # session 3
                         1, 1, 2, 2,
#                          3, 3, 4, 4,
#                          2, 2, 3, 3,
                         2, 2, 4, 4,
                         1, 1, 4, 4
                         ),
                   y = c(3.2, 3.5, 3.5, 3.2, 
                         4.4, 4.7, 4.7, 4.4, 
                         5.6, 5.9, 5.9, 5.6,
#                          3.2, 3.5, 3.5, 3.2,
#                          7.2, 7.5, 7.5, 7.2,
                         # session 3
                         7.2, 7.5, 7.5, 7.2,
#                          7.2, 7.5, 7.5, 7.2,
                         8.4, 8.7, 8.7, 8.4,
                         9.6, 9.9, 9.9, 9.6#,
#                          10.8, 11.1, 11.1, 10.8
                         ), 
#                    genotype=c(rep("WT",4), rep("TS",4), rep("WTEEEGCG",4), rep("TSEEEGCG",4), rep("TS_fake",4),
                     genotype=c(rep("WT",4), rep("TS",4), rep("WTEEEGCG",4),
#                               rep("WT",4), rep("TS",4), rep("TSEEEGCG",4), rep("WTEEEGCG",4), rep("TS_fake",4)),
                                rep("WT",4), rep("TS",4), rep("TSEEEGCG",4)),
                   day=c(rep("Session 1", 12), rep("Session 3", 12)))

sl_1$genotype <- factor(sl_1$genotype, levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"), 
                        labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"))

stars_plot <- data.frame(x_pos = c(1.5, 3,  2.5,
                                  #c(1.5, 2.5, 3, 3.5, 2.5,
                                   1.5, 3, 2.5), 
                         y_pos = c(3.6, 4.8, 6, #3.6, 7.6, 
#                                    7.6, 7.6, 8.8, 10, 11.2),
                                   7.6, 8.8, 10),
                         label=c("***","*", "*",
                               #c("***","***","*","**", "*",
#                                  "***", "***", "***", "*", "*"),
                                "***", "*", "*"),
#                          day=c(rep("Session 1", 5), rep("Session 3", 5)),
                         day=c(rep("Session 1", 3), rep("Session 3", 3)),
#                          genotype=c(rep("WT", 10)))         
                         genotype=c(rep("WT", 6)))   

# p_cloud_indiv_by_day_facet <- boxPlots.PC1 + facet_wrap(~genotype, ncol = 2)
boxPlots.PC1.rev.line <-boxPlots.PC1.rev.facet + geom_hline(yintercept = 0, colour="gray") +
  theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis)) +
  panel_border() + theme(panel.border = element_rect(colour = "black"))

boxPlots.PC1.rev.line.stars <- boxPlots.PC1.rev.line + geom_path(data = sl_1, aes(x = x, y = y)) +                  
  geom_text(data=stars_plot, aes(x=x_pos, y=y_pos, label=label), 
            size=10, show.legend = FALSE) +
  theme(strip.text.x = element_text(size=size_strips, face="bold")) +
  theme(plot.title = element_text(size=size_titles)) + 
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis)) +
  theme(strip.background = element_blank()) + 
  panel_border() + theme(panel.border = element_rect(colour = "black"))                                                                                                      

boxPlots.PC1.rev.line.stars

#######################
#######
## Panel with boxplots
panel_boxPlots <- ggdraw() + draw_plot(p_cloud_rev1_3_facet_strips_coord, 0, .5, 1, .5) +
  draw_plot(boxPlots.PC1.rev.line.stars, 0, 0, 1, .5) +
#   draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = size_titles)
  draw_plot_label(c("E", "F"), c(0, 0), c(1, 0.5), size = size_titles)
panel_boxPlots

ggsave (panel_boxPlots, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "panel_boxPlot_rev", img_format, sep=""), 
        dpi=dpi_q, width=15, height=11)
# size 1100, 700


#####################
####################
###################

# Plotting a legend with scuares and colors
l <- ggplot() + geom_point(data=PC1.rev3, aes (x=PC1, y=PC2, colour = genotype), shape=15, size=5) +
#   scale_colour_manual (values=c("red", "darkgreen", "magenta", "black"))
     scale_colour_manual (values=c("red", "darkgreen", "magenta", "green"))
l <- l + guides(color=guide_legend(title=NULL)) 
l <- l + theme(legend.key = element_blank())
l

# ggsave (l, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_reversal/", "legend_squares.jpg", sep=""), dpi=900)


