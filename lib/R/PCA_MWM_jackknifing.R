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
  PC1_v <- var_coord [, PC]
  
  return (PC1_v)
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
result_jk$angle_degrees <- result_jk$angle * 180/3.1416

colnames (result_jk) <- c("dropped_id", "angle_rad", "angle_degree")
setwd("/Users/jespinosa/20151001_ts65_young_MWM/tbl/")
wd <- getwd()
write.table(result_jk, file = paste(wd, "/jackknife_angles.csv", sep=""), sep="\t", row.names=FALSE, col.names=T)



# t.values[p] <- t_s




