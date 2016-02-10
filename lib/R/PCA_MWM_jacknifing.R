#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Febraury 2016                    ###
#############################################################################
### MWM young PCA jacknifing                                              ###
### Perform a PCA removing each time an animal                            ###
### Assess the stability of the PCA results                               ###
### Cluster runs                                                          ###
### In this case what I can do is the same as in the permutation test     ###
### but in this case what I would have to provide the script is which     ###
### animal has to be drop out from the PCA                                ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

# To use this script in ant first export this:
# export R_LIBS="/software/R/packages"

##Loading libraries
## Local runs
# library(FactoMineR)
# library(Hmisc)
## Cluster
library(FactoMineR, lib.loc="/users/cn/jespinosa/R/library")
library(Hmisc, lib.loc="/users/cn/jespinosa/R/library")

#####################
### VARIABLES
#Reading arguments
args <- commandArgs (TRUE) #if not it doesn't start to count correctly

## Default setting when no arguments passed
if ( length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      PCA_perm_test_nx.R
      
      Arguments:
      --id_drop=someValue       - int
      --path2file=someValue  - character, path to read files
      --help                 - print this text
      
      Example:
      ./PCA_perm_test_nx.R --id_drop=\"id_drop\" \n")
  
  q (save="no")
}

# Use to parse arguments beginning by --
parseArgs <- function(x) 
{
  strsplit (sub ("^--", "", x), "=")
}

#Parsing arguments
argsDF <- as.data.frame (do.call("rbind", parseArgs(args)))
argsL <- as.list (as.character(argsDF$V2))
names (argsL) <- argsDF$V1
print (argsL)

# id_drop is mandatory
{
  if (is.null (argsL$id_drop)) 
  {
    stop ("[FATAL]: id_drop parameter is mandatory")
  }
  else
  {
    id_drop <- as.integer(argsL$id_drop)
  }
}

{
  if (is.null (argsL$path2files)) 
  {
    stop ("[FATAL]: Path to table is mandatory")
    #     path2files <- "/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/Ts65Dn_OLD_ACQ1_ACQ5_SUBCONJ.sav"
  }
  else
  {
    path2files <- argsL$path2files
  }
}

# The data has to be in a general format if I get it in sav format first transform it outside to a csv format
# ma2=spss.get(path2files)
# ma2 <- read.csv("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/ts65_old_3sup_tsegcg_rev.csv", sep="\t")
# ma2 <- read.csv("/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv", sep="\t")
ma2 <- read.csv(path2files, sep="\t")

# New function taking into account the acquisition days in the comparison
# f_t_stat <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=5){
#   group1 <- subset (new_coord, genotype == gen_1 & day==acq_day)
#   group2 <- subset (new_coord, genotype == gen_2 & day==acq_day)
#   t_stat = t.test(group1$V1, group2$V1)$statistic  
#   return <- c(t_stat, gen_1, gen_2, paste (gen_1, gen_2, sep="_"))
# }
ma2
# id_group <- subset(ma2, grepl("Day 1", ma2$day))

## Get all unique mouse id from the table in bash
## tail -n +2  "/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv" | awk -F '\t' '{print $2}' | sort | uniq

# Drop mouse id from the table
id_drop <- "130014264"
id_drop <- as.numeric(id_drop)

ma2_jk <- subset(ma2, id != id_drop)
ma2_jk 

get_pc <- function (df, ori_PC, PC="Dim.1") {
  n_col <- dim (df)[2]
  variables_list <- colnames(ma2_jk) [4:n_col] 
  
  # Variables interpolated
  tbl_median <- with (df, aggregate (mget(variables_list), list (gentreat=gentreat, day=day), FUN=median))
  
  n_col_tbl_all <- length(tbl_median[1,])
  
  ## PCA of the medians
  res = PCA(tbl_median[,(3:n_col_tbl_all)], scale.unit=TRUE, graph=F) 
  
  var_coord <- as.data.frame (res$var$coord)
  PC1_v <- var_coord [, PC]
  
  return (PC1_v)
}


for i in unique (ma2$id){
  
}
ori_pc1<- get_pc (ma2)

library(plyr)
tbl_jk <- do.call("cbind", ddply(ma2, c("id"), function(x) { 
  get_pc(subset(ma2, id != x$id))}))

df_jk <- as.data.frame(tbl_jk)


ddply (df_jk, c("id"), function(x) { angle_bw_v (ori_pc1, as.vector(subset(df_jk, id == x$id, select=c(2:8))))})


subset(df_jk, id == "130014754")
ddply (df_jk, c("id"), function(x) { print(x$id)})

df_jk$id

# theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )

angle_bw_v <- function (a, b){
  theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  return (theta)
}

a <- ori_pc1
b <- as.vector(df_jk [2,c(2:8) ])
angle_bw_v (a,b)


theta

ddply(ma2, .(id), function(x)  subset(ma2, id != x$id))

      data.frame(mean=mean(x[,2])))

n_col <- dim (ma2_jk)[2]
variables_list <- colnames(ma2_jk) [4:n_col] 

# Variables interpolated
tbl_median <- with (ma2_jk, aggregate (mget(variables_list), list (gentreat=gentreat, day=day), FUN=median))
# tbl_median


### I don't need in the jacknifing at least by the moment
# I get all the individuals values from the tbl with the mouse deleted to project them in the pca space
# tbl_ind <- ma2_jk
# 
# # head (tbl_median)
# # head(tbl_ind[,-1])
# # rownames (tbl_ind) <- paste (tbl_ind[,2], gsub ("Day ", "", tbl_ind[,3]), sep="")
# 
# tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
# n_median <- length(tbl_median[,1])
n_col_tbl_all <- length(tbl_median[1,])

## PCA with mice projected in the new space 
# res = PCA(tbl_med_ind[,(3:n_col_tbl_all)], scale.unit=TRUE, ind.sup=c((n_median+1):length(tbl_med_ind[,1])), graph=F) 

## PCA of the medians
res = PCA(tbl_median[,(3:n_col_tbl_all)], scale.unit=TRUE, graph=T) 

var_coord <- as.data.frame (res$var$coord)
PC1_v <- var_coord$Dim.1


labels_v <- row.names(res$var$coord)
which (circle_plot$Dim.1 < 0)


# Assigning again the labels of the groups for t statistic calculation
#res$var$coord[,1]
#new_coord <- cbind(res$ind.sup$coord[,1], res$ind.sup$coord[,2])
day <- c()
genotype <- c()
new_coord <- as.data.frame(cbind(res$ind.sup$coord[,1],res$ind.sup$coord[,2]))
new_coord$day <- gsub ("Day ", "", tbl_ind$day)
new_coord$genotype <- tbl_ind$gentreat

new_coord$genotype <- factor(new_coord$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                             labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))


# myset=which(iddaytreat[,2]==tgt[1] & iddaytreat[,3]==acq[1])
# myset<-c()
# for (j in 1:length(tgt)){
#   for (i in 1:length(acq)){
#     if (i<5){
#       ind=(i-1)*length(tgt)+j
#       ind2=i*length(tgt)+j 
#       
#     }
#     myset=which(iddaytreat[,2]==tgt[j] & iddaytreat[,3]==acq[i])
#     new_coord [myset,c("day")] <- i
#     new_coord [myset,c("genotype")] <- tgt[j]
#     new_coord [myset,c("id")] <- substr(iddaytreat[myset,1], 7, 9)    
#   }
# }

# Get all possible combinations of genontype treatment pairwise comparisons 
gentreat <- sort(unique(ma2$gentreat))
gentreat_pairs <- t (combn (gentreat,2))

result <- c()
result_1 <- c()

for (row in 1:length(gentreat_pairs [,1])) {
  gr1 <- as.character(gentreat_pairs [row,1])
  gr2 <- as.character(gentreat_pairs [row,2])
  
  result_v <- c(f_t_stat (new_coord, gr1, gr2, acq_day=5), seed)
  result_v_1 <- c(f_t_stat (new_coord, gr1, gr2, acq_day=1), seed)
  
  result <- rbind  (result, result_v)
  result_1 <- rbind  (result_1, result_v_1)
  
  colnames(result) <- c ("t", "gr1", "gr2", "comparison", "seed") 
  colnames(result_1) <- c ("t", "gr1", "gr2", "comparison", "seed")               
}

wd <- getwd()
write.table(result, file = paste(wd, "/tbl_t_stat_day5.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)
write.table(result_1, file = paste(wd, "/tbl_t_stat_day1.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)






# t.values[p] <- t_s




