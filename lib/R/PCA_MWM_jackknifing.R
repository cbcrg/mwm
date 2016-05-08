#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. February 2016                    ###
#############################################################################
### MWM young PCA jacknifing                                              ###
### Perform a PCA removing each time an animal                            ###
### Assess the stability of the PCA results by computing the angle        ###
### between the original PC and the one when the animal is removed        ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

# To use this script in ant first export this:
# export R_LIBS="/software/R/packages"

##Loading libraries
## Local runs
library(FactoMineR)
library(Hmisc)
library(plyr)
library(cowplot)

# Variables
img_format <- ".tiff"
size_titles <- 18
size_axis <- 14


## Cluster
# library(FactoMineR, lib.loc="/users/cn/jespinosa/R/library")
# library(Hmisc, lib.loc="/users/cn/jespinosa/R/library")

# The data has to be in a general format if I get it in sav format first transform it outside to a csv format
# ma2=spss.get(path2files)
# ma2 <- read.csv("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/ts65_old_3sup_tsegcg_rev.csv", sep="\t")
# ma2 <- read.csv("/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv", sep="\t")
path2files <- "/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv"
ma2 <- read.csv(path2files, sep="\t")

## Get all unique mouse id from the table in bash
## tail -n +2  "/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv" | awk -F '\t' '{print $2}' | sort | uniq

###############
### Functions
# Function to get a PC from the original medians data frame or from the data frame with all rows from an individual deleted
get_pc <- function (df, ori_PC, PC="Dim.1") {
  n_col <- dim (df)[2]
  variables_list <- colnames(df) [4:n_col] 
  
  # Variables interpolated
  tbl_median <- with (df, aggregate (mget(variables_list), list (gentreat=gentreat, day=day), FUN=median))
  
  n_col_tbl_all <- length(tbl_median[1,])
  
  ## PCA of the medians
  res = PCA(tbl_median[,(3:n_col_tbl_all)], scale.unit=TRUE, graph=F) 
  
  var_coord <- as.data.frame (res$var$coord)
  PC_v <- var_coord [, PC]
  
  return (PC_v)
}


# Function to calculate angle between two vectors
# theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
angle_bw_v <- function (a, b){
  theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  return (theta)
}

# PC1 original data
ori_pc1<- get_pc (ma2)

# Calculating PC1 for each one out data frame
tbl_jk <- do.call("cbind", ddply(ma2, c("id"), function(x) { 
  get_pc(subset(ma2, id != x$id))}))

df_jk <- as.data.frame(tbl_jk)

# Calculating angle between all PC1 and the original one
result_jk <- ddply (df_jk, c("id"), function(x) { angle_bw_v (ori_pc1, as.vector(subset(df_jk, id == x$id, select=c(2:8))))})

# Angles are in radians
result_jk$angle_degrees <- result_jk$V1 * 180/3.1416

colnames (result_jk) <- c("dropped_id", "angle_rad", "angle_degree")
setwd("/Users/jespinosa/20151001_ts65_young_MWM/tbl/")
wd <- getwd()
# write.table(result_jk, file = paste(wd, "/jackknife_angles.csv", sep=""), sep="\t", row.names=FALSE, col.names=T)

#################
# PC2 jackknifing -> each time one id is left out
# PC2 original data
ori_pc2<- get_pc (ma2, PC="Dim.2")

# Calculating PC1 for each one out data frame
tbl_jk_PC2 <- do.call("cbind", ddply(ma2, c("id"), function(x) { 
  get_pc(subset(ma2, id != x$id), PC="Dim.2")}))

df_jk_PC2 <- as.data.frame(tbl_jk_PC2)

# Calculating angle between all PC1 and the original one
result_jk_PC2 <- ddply (df_jk_PC2, c("id"), function(x) { angle_bw_v (ori_pc2, as.vector(subset(df_jk_PC2, id == x$id, select=c(2:8))))})

# Angles are in radians
result_jk_PC2$angle_degrees <- result_jk_PC2$V1 * 180/3.1416

colnames (result_jk_PC2) <- c("dropped_id", "angle_rad", "angle_degree")
setwd("/Users/jespinosa/20151001_ts65_young_MWM/tbl/")
wd <- getwd()
# write.table(result_jk_PC2, file = paste(wd, "/jackknife_angles_PC2.csv", sep=""), sep="\t", row.names=FALSE, col.names=T)

###########
### Plots

#####
# PC1
histogram(result_jk$angle_degree)

#######################
# Normal histogram ggplot2
ggplot(result_jk, aes(angle_degree)) + geom_histogram(binwidth=0.1, colour="black", fill="darkgray") +
                                       scale_x_continuous(limits=c(0.10, 1.14), breaks=(seq(0.10, 1.14, by=0.1))) +
  labs(title = "Angle distribution\n", x = paste("\nAngle (degrees)", sep=""), 
       y=paste("Count\n", sep = "")) + panel_border() + theme(panel.border = element_rect(colour = "black"))

PC1_hist<- ggplot(result_jk, aes(angle_degree)) + geom_histogram(binwidth=0.1, colour="black", fill="darkgray") +
  scale_x_continuous(limits=c(0.10, 1.14), breaks=(seq(0.10, 1.14, by=0.1))) +
  labs(title = "Angle distribution to respect original PC1\n", x = paste("\nAngle (degrees)", sep=""), 
       y=paste("Counts\n", sep = "")) + panel_border() + theme(panel.border = element_rect(colour = "black")) +
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis))

# ggsave (PC1_hist, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_jackknifing/", "hist_pc1.tiff", sep=""), 
#         width = 8.5, height = 6, dpi=900)

#####
# PC2
histogram(result_jk_PC2$angle_degree)

#######################
# Normal histogram ggplot2 PC2
PC2_hist <- ggplot(result_jk_PC2, aes(angle_degree)) + geom_histogram(binwidth=1, colour="black", fill="darkgray") +
  scale_x_continuous(limits=c(0, 8), breaks=(seq(0, 8, by=2))) +
  labs(title = "Angle distribution to respect original PC2\n", x = paste("\nAngle (degrees)", sep=""), 
       y=paste("Counts\n", sep = "")) + panel_border() + theme(panel.border = element_rect(colour = "black")) +
  theme(axis.title.x = element_text(size=size_axis)) +
  theme(axis.title.y = element_text(size=size_axis))

# ggsave (PC2_hist, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_jackknifing/", "hist_pc2.jpg", sep=""), 
#         width = 8.5, height = 6, dpi=900)

panel_hist_PC <- ggdraw() + draw_plot (PC1_hist, 0, 0.5, 1, 0.5) +
  draw_plot (PC2_hist, 0, 0, 1, 0.5) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = size_titles)
panel_hist_PC

ggsave (panel_hist_PC, file=paste(home, "/20151001_ts65_young_MWM/figures/fig_jackknifing/", "hist_pc_panel", img_format, sep=""), 
        width = 8.5, height = 12, dpi=300)

##########
# Together
result_jk$PC <- 1
result_jk_PC2$PC <- 2
result_PC1_2 <- rbind(result_jk, result_jk_PC2)
result_PC1_2$PC <- as.factor(result_PC1_2$PC)

# This doesn't work 
# s telling ggplot to construct one histogram using all the values in f0 and 
then color the bars of this single histogram according to the variable utt.
ggplot(result_PC1_2, aes(angle_degree, fill=PC)) + 
geom_histogram(binwidth=0.2,  alpha=0.2, position="dodge") #+ scale_fill_identity()

ggplot(result_PC1_2, aes(angle_degree, fill=PC)) + 
  geom_histogram(binwidth=0.2,  alpha=0.2)

# WORKING other way not dodge
ggplot(result_PC1_2,aes(x=angle_degree)) + 
  geom_histogram(data=subset(result_PC1_2,PC == '1'),binwidth=0.2,colour="black", fill = "red", alpha = 0.2, show.legend =T) +
  geom_histogram(data=subset(result_PC1_2,PC == '2'),binwidth=0.2,colour="black", fill = "blue", alpha = 0.2, show.legend =T) 

ggplot(result_PC1_2, aes(angle_degree)) + 
  geom_histogram(data = result_PC1_2, fill = "blue", alpha = 0.2) +
  geom_histogram(data = result_jk, fill = "red", alpha = 0.2) 

+
  scale_x_continuous(limits=c(0.10, 7), breaks=(seq(0.10, 7, by=0.2))) +
  labs(title = "Angle distribution\n", x = paste("\nAngle to original PC1 (degrees)", sep=""), 
       y=paste("Count\n", sep = "")) + panel_border() + theme(panel.border = element_rect(colour = "black")) 

+
  geom_histogram(data=result_jk_PC2, aes(angle_degree), binwidth=0.2, colour="black", fill="lightgray")

dat <- data.frame(xx = c(runif(100,20,50),runif(100,40,80),runif(100,0,30)),yy = rep(letters[1:3],each = 100))

ggplot(dat,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == 'a'),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(dat,yy == 'b'),fill = "blue", alpha = 0.2) +
  geom_histogram(data=subset(dat,yy == 'c'),fill = "green", alpha = 0.2)

############
# More examples of fancy histograms
ggplot(result_jk, aes(angle_degree)) + geom_density(colour="red", fill="red") +
  scale_x_continuous(limits=c(0, 1.30), breaks=(seq(0, 1.30, by=0.1)))


ggplot(result_jk, aes(x=angle_degree)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.1,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  scale_x_continuous(limits=c(0, 1.30), breaks=(seq(0, 1.30, by=0.1)))

ggplot(result_jk, aes(x=angle_degree)) + 
  geom_density(aes(x=angle_degree, y=..scaled.., fill=grp), alpha=.2, fill="#FF6666") +
  scale_x_continuous(limits=c(0, 1.30), breaks=(seq(0, 1.30, by=0.1)))

ggplot(result_jk, aes(x=angle_degree)) + geom_histogram(aes(y=..count..)) 

ggplot(result_jk, aes(x=angle_degree)) + geom_histogram(aes(y=..count..)) +
geom_density(aes(x=angle_degree, y=..count..), alpha=.2, fill="#FF6666")

### counts frequency
ggplot(result_jk, aes(x=angle_degree)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.1,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  scale_x_continuous(limits=c(0, 1.30), breaks=(seq(0, 1.30, by=0.1)))

# version final
ggplot(result_jk, aes(x=angle_degree)) + 
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=.1,
                 colour="black", fill="white") +
  scale_x_continuous(limits=c(0.15, 1.14), breaks=(seq(0.15, 1.14, by=0.1))) +
  geom_density(aes(x=angle_degree, y=..density..), alpha=.2, fill="#FF6666")


require(scales)
p <- ggplot(result_jk, aes(x = angle_degree)) +  
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  ## version 3.0.9
  # scale_y_continuous(labels = percent_format())
  ## version 3.1.0
  scale_y_continuous(labels=percent) +
  scale_x_continuous(limits=c(0.15, 1.14), breaks=(seq(0.15, 1.14, by=0.1)))
p

min(result_jk$angle_degree)
max(result_jk$angle_degree)
seq(-5, 10, by=5)