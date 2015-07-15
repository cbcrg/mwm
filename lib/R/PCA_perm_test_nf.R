#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. July 2015                        ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes                 ###
### Cluster runs                                                          ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

# To use this script in ant first export this:
# export R_LIBS="/software/R/packages"

##Loading libraries
#library(FactoMineR)
#library(Hmisc)
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
      --seed=someValue       - int
      --path2file=someValue  - character, path to read files
      --help                 - print this text
      
      Example:
      ./PCA_perm_test_nx.R --seed=\"seed\" \n")
  
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

# seed is mandatory
{
  if (is.null (argsL$seed)) 
  {
    stop ("[FATAL]: seed parameter is mandatory")
  }
  else
  {
    seed <- as.integer(argsL$seed)
  }
}

{
  if (is.null (argsL$path2files)) 
  {
    path2files <- "/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/Ts65Dn_OLD_ACQ1_ACQ5_SUBCONJ.sav"
  }
  else
  {
    path2files <- argsL$path2files
  }
}

ma2=spss.get(path2files)

# Function t statistic calculation
# f_t_stat <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG"){
#   group1 <- subset (new_coord, genotype == gen_1)
#   group2 <- subset (new_coord, genotype == gen_2)
#   t_stat = t.test(group1$V1, group2$V1)$statistic
# 
#   return <- c(t_stat, gen_1, gen_2, paste (gen_1, gen_2, sep="_"))
# }

# New function taking into account the acquisition days in the comparison
f_t_stat <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=5){
  group1 <- subset (new_coord, genotype == gen_1 & day==acq_day)
  group2 <- subset (new_coord, genotype == gen_2 & day==acq_day)
  t_stat = t.test(group1$V1, group2$V1)$statistic  
  return <- c(t_stat, gen_1, gen_2, paste (gen_1, gen_2, sep="_"))
}

id_group <- subset(ma2, grepl("PRE", ma2$day))

id_group <- id_group [,c(1,2)]

# print (paste("ddddddddd", seed))
set.seed (seed)

p_genetreat <- sample(id_group$gentreat)
id_group$gentreat <- p_genetreat

p_ma2 <- merge(id_group, ma2, by="id", all = TRUE)

# I eliminiate the original column with groups labels
p_ma2 <- p_ma2 [,c(-3)]
#head (p_ma2, 50)
#head (ma2, 50)

M=as.matrix(p_ma2[p_ma2$day=="Day 5",c(4,5,7:12)])
treat=as.matrix(p_ma2[p_ma2$day=="Day 5",2])
M4=as.matrix(p_ma2[p_ma2$day=="Day 4",c(4,5,7:12)])
treat4=as.matrix(p_ma2[p_ma2$day=="Day 4",2])
M3=as.matrix(p_ma2[p_ma2$day=="Day 3",c(4,5,7:12)])
treat3=as.matrix(p_ma2[p_ma2$day=="Day 3",2])
M2=as.matrix(p_ma2[p_ma2$day=="Day 2",c(4,5,7:12)])
treat2=as.matrix(p_ma2[p_ma2$day=="Day 2",2])
M1=as.matrix(p_ma2[p_ma2$day=="Day 1",c(4,5,7:12)])
treat1=as.matrix(p_ma2[p_ma2$day=="Day 1",2])

MM5=matrix(0,8,8)

MM5[1,]=apply(M[treat=="WT",],2,median)
MM5[2,]=apply(M[treat=="WTEE",],2,median)
MM5[3,]=apply(M[treat=="WTEGCG",],2,median)
MM5[4,]=apply(M[treat=="WTEEEGCG",],2,median)
MM5[5,]=apply(M[treat=="TS",],2,median)
MM5[6,]=apply(M[treat=="TSEE",],2,median)
MM5[7,]=apply(M[treat=="TSEGCG",],2,median)
MM5[8,]=apply(M[treat=="TSEEEGCG",],2,median)
MM1=matrix(0,8,8)
MM1[1,]=apply(M1[treat1=="WT",],2,median)
MM1[2,]=apply(M1[treat1=="WTEE",],2,median)
MM1[3,]=apply(M1[treat1=="WTEGCG",],2,median)
MM1[4,]=apply(M1[treat1=="WTEEEGCG",],2,median)
MM1[5,]=apply(M1[treat1=="TS",],2,median)
MM1[6,]=apply(M1[treat1=="TSEE",],2,median)
MM1[7,]=apply(M1[treat1=="TSEGCG",],2,median)
MM1[8,]=apply(M1[treat1=="TSEEEGCG",],2,median)
MM2=matrix(0,8,8)
MM2[1,]=apply(M2[treat2=="WT",],2,median)
MM2[2,]=apply(M2[treat2=="WTEE",],2,median)
MM2[3,]=apply(M2[treat2=="WTEGCG",],2,median)
MM2[4,]=apply(M2[treat2=="WTEEEGCG",],2,median)
MM2[5,]=apply(M2[treat2=="TS",],2,median)
MM2[6,]=apply(M2[treat2=="TSEE",],2,median)
MM2[7,]=apply(M2[treat2=="TSEGCG",],2,median)
MM2[8,]=apply(M2[treat2=="TSEEEGCG",],2,median)
MM3=matrix(0,8,8)
MM3[1,]=apply(M3[treat3=="WT",],2,median)
MM3[2,]=apply(M3[treat3=="WTEE",],2,median)
MM3[3,]=apply(M3[treat3=="WTEGCG",],2,median)
MM3[4,]=apply(M3[treat3=="WTEEEGCG",],2,median)
MM3[5,]=apply(M3[treat3=="TS",],2,median)
MM3[6,]=apply(M3[treat3=="TSEE",],2,median)
MM3[7,]=apply(M3[treat3=="TSEGCG",],2,median)
MM3[8,]=apply(M3[treat3=="TSEEEGCG",],2,median)
MM4=matrix(0,8,8)
MM4[1,]=apply(M4[treat4=="WT",],2,median)
MM4[2,]=apply(M4[treat4=="WTEE",],2,median)
MM4[3,]=apply(M4[treat4=="WTEGCG",],2,median)
MM4[4,]=apply(M4[treat4=="WTEEEGCG",],2,median)
MM4[5,]=apply(M4[treat4=="TS",],2,median)
MM4[6,]=apply(M4[treat4=="TSEE",],2,median)
MM4[7,]=apply(M4[treat4=="TSEGCG",],2,median)
MM4[8,]=apply(M4[treat4=="TSEEEGCG",],2,median)
MMALL=rbind(MM1,MM2,MM3,MM4,MM5)
colnames(MMALL)=c("dist","gallindex","latency","speed","percentne","percenter","wishaw","percentperi")
rownames(MMALL)=c("WT1","WTEE1","WTEGCG1","WTEEEGCG1","TS1","TSEE1","TSEGCG1","TSEEEGCG1","WT2","WTEE2","WTEGCG2","WTEEEGCG2","TS2","TSEE2","TSEGCG2","TSEEEGCG2","WT3","WTEE3","WTEGCG3","WTEEEGCG3","TS3","TSEE3","TSEGCG3","TSEEEGCG3","WT4","WTEE4","WTEGCG4","WTEEEGCG4","TS4","TSEE4","TSEGCG4","TSEEEGCG4","WT5","WTEE5","WTEGCG5","WTEEEGCG5","TS5","TSEE5","TSEGCG5","TSEEEGCG5")

# All individuals as supplementary points
acq=c("Day 1","Day 2","Day 3","Day 4","Day 5")
var=c(4,5,7,8,9,11,12)

M.ind=as.matrix(p_ma2[p_ma2$day%in%acq,var])

# Filtering PRE and CUE
iddaytreat=as.matrix(p_ma2[p_ma2$day%in%acq,1:3])
rownames(iddaytreat)=c(1:415)
tgt=sort(unique(iddaytreat[,2]))

M.med=matrix(0,40,7)
rnames=c(1:40)
for (i in 1:length(acq)){
  for (j in 1:length(tgt)){
    ind=(i-1)*length(tgt)+j
    # print (ind)
    M.med[ind,]=apply(M.ind[iddaytreat[,3]==acq[i] & iddaytreat[,2]==tgt[j],],2,median)
    rnames[ind]=paste(tgt[j],acq[i]," ")
  }
}
colnames(M.med)=colnames(M.ind)
rownames(M.med)=rnames

jm=rbind(M.med,M.ind)
# jm[c(1:40),]

res = PCA (jm, scale.unit=TRUE, ind.sup=c(41:455), graph=F) 

# Assigning again the labels of the groups for t statistic calculation
#res$var$coord[,1]
#new_coord <- cbind(res$ind.sup$coord[,1], res$ind.sup$coord[,2])
day <- c()
genotype <- c()
new_coord <- as.data.frame(cbind(res$ind.sup$coord[,1],res$ind.sup$coord[,2]))
new_coord$day <- c(1:415)
new_coord$genotype <- c()

# myset=which(iddaytreat[,2]==tgt[1] & iddaytreat[,3]==acq[1])
myset<-c()
for (j in 1:length(tgt)){
  for (i in 1:length(acq)){
    if (i<5){
      ind=(i-1)*length(tgt)+j
      ind2=i*length(tgt)+j 
      
    }
    myset=which(iddaytreat[,2]==tgt[j] & iddaytreat[,3]==acq[i])
    new_coord [myset,c("day")] <- i
    new_coord [myset,c("genotype")] <- tgt[j]
    new_coord [myset,c("id")] <- substr(iddaytreat[myset,1], 7, 9)    
  }
}

# Get all possible combinations of genontype treatment pairwise comparisons 
gentreat <- unique(ma2$gentreat)
gentreat_pairs <- t (combn (gentreat,2))

result <- c()
result_acq1 <- c()

for (row in 1:length(gentreat_pairs [,1])) {
  gr1 <- as.character(gentreat_pairs [row,1])
  gr2 <- as.character(gentreat_pairs [row,2])
  
  result_v <- c(f_t_stat (new_coord, gr1, gr2, acq_day=5), seed)
  result_v_1 <- c(f_t_stat (new_coord, gr1, gr2, acq_day=1), seed)
  
  result <- rbind  (result, result_v)
  result_1 <- rbind  (result, result_v_1)
  
  colnames(result) <- c ("t", "gr1", "gr2", "comparison", "seed") 
  colnames(result_1) <- c ("t", "gr1", "gr2", "comparison", "seed")               
}


# t_s_ts_tseeegcg <- f_t_stat (new_coord, "TS", "TSEEEGCG")
# 
# t_s_ts_tsee <- f_t_stat (new_coord, "TS", "TSEE")
# 
# t_s_ts_tsegcg <- f_t_stat (new_coord, "TS", "TSEGCG")
# 
# t_s_tseeegcg_tsegcg <- f_t_stat (new_coord, "TSEEEGCG", "TSEGCG")
# 
# t_s_tseeegcg_tsee <- f_t_stat (new_coord, "TSEEEGCG", "TSEE")
# 
# t_s_tsee_tsegcg <- f_t_stat (new_coord, "TSEE", "TSEGCG")
# 
# # t_s_ts_tseeegcg <- f_t_stat (new_coord, "TS", "TSEEEGCG")
# #1
# t_s_ts_wt <- f_t_stat (new_coord, "TS", "WT")
# #3
# t_s_ts_wtee <- f_t_stat (new_coord, "TS", "WTEE")
# #8
# t_s_ts_wtegcg <- f_t_stat (new_coord, "TS", "WTEGCG")
# t_s_ts_wteeegcg <- f_t_stat (new_coord, "TS", "WTEEEGCG")
# #2
# t_s_wt_wtee <- f_t_stat (new_coord, "WT", "WTEE")
# #7
# t_s_wt_wtegcg <- f_t_stat (new_coord, "WT", "WTEGCG")
# t_s_wt_wteeegcg <- f_t_stat (new_coord, "WT", "WTEEEGCG")
# #4
# t_s_wt_tsee <- f_t_stat (new_coord, "WT", "TSEGCG")
# #11
# t_s_wt_tsegcg <- f_t_stat (new_coord, "WT", "TSEE")
# t_s_wt_tseeegcg <- f_t_stat (new_coord, "WT", "TSEEEGCG")
# #9
# t_s_wtee_wtegcg <- f_t_stat (new_coord, "WTEE", "WTEGCG")
# t_s_wtee_wteeegcg <- f_t_stat (new_coord, "WTEE", "WTEEEGCG")
# #6
# t_s_wtee_tsee <- f_t_stat (new_coord, "WTEE", "TSEE")
# t_s_wtee_tsegcg <- f_t_stat (new_coord, "WTEE", "TSEGCG")
# t_s_wtee_tseeegcg <- f_t_stat (new_coord, "WTEE", "TSEEEGCG")
# t_s_wtegcg_wteeegcg <- f_t_stat (new_coord, "WTEGCG", "WTEEEGCG")
# #10
# t_s_wtegcg_tsee <- f_t_stat (new_coord, "WTEGCG", "TSEE")
# t_s_wtegcg_tsegcg <- f_t_stat (new_coord, "WTEE", "TSEGCG")
# t_s_wtegcg_tseeegcg <- f_t_stat (new_coord, "WTEGCG", "TSEEEGCG")
# t_s_wteeegcg_tsee <- f_t_stat (new_coord, "WTEEEGCG", "TSEE")
# t_s_wteeegcg_tsegcg <- f_t_stat (new_coord, "WTEEEGCG", "TSEGCG")
# t_s_wteeegcg_tseeegcg <- f_t_stat (new_coord, "WTEEEGCG", "TSEEEGCG")
# 
# # result <- rbind (t_s_ts_tseeegcg, t_s_ts_tsee, t_s_ts_tsegcg, t_s_ts_wt, t_s_ts_wtee, t_s_ts_wtegcg, t_s_ts_wteeegcg, 
# #                  t_s_wt_wtee,t_s_wt_wtegcg,t_s_wt_wteeegcg ,t_s_wt_tsee,t_s_wt_tsegcg,t_s_wt_tseeegcg,t_s_wtee_wtegcg,
# #                  t_s_wtee_wteeegcg,t_s_wtee_tsee,t_s_wtee_tsegcg,t_s_wtee_tseeegcg,t_s_wtegcg_wteeegcg,t_s_wtegcg_tsee,
# #                  t_s_wtegcg_tsegcg,t_s_wtegcg_tseeegcg,t_s_wteeegcg_tsee,t_s_wteeegcg_tsegcg,t_s_wteeegcg_tseeegcg)
# 
# result <- rbind (t_s_ts_tseeegcg, t_s_ts_tsee, t_s_ts_tsegcg, t_s_tseeegcg_tsegcg, t_s_tseeegcg_tsee, t_s_tsee_tsegcg, 
#                  t_s_ts_wt, t_s_ts_wtee, t_s_ts_wtegcg, t_s_ts_wteeegcg, t_s_wt_wtee, t_s_wt_wtegcg, t_s_wt_wteeegcg, 
#                  t_s_wt_tsee, t_s_wt_tsegcg, t_s_wt_tseeegcg, t_s_wtee_wtegcg, t_s_wtee_wteeegcg, t_s_wtee_tsee, 
#                  t_s_wtee_tsegcg, t_s_wtee_tseeegcg, t_s_wtegcg_wteeegcg, t_s_wtegcg_tsee, t_s_wtegcg_tsegcg,
#                  t_s_wtegcg_tseeegcg, t_s_wteeegcg_tsee, t_s_wteeegcg_tsegcg, t_s_wteeegcg_tseeegcg)
# 
# result <- cbind  (result, seed)
# 
# colnames(result) <- c ("t", "gr1", "gr2", "comparison", "seed")
wd <- getwd()
write.table(result, file = paste(wd, "/tbl_t_stat.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)
write.table(result_1, file = paste(wd, "/tbl_t_stat_day1.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)






# t.values[p] <- t_s



