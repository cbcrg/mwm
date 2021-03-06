#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. July 2015                        ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes                 ###
###                                                                       ###
#############################################################################

#####################
### VARIABLES
library(FactoMineR)
library(Hmisc)
home <- Sys.getenv("HOME")

# Function t statistic calculation real data 
f_t_stat_only_real <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=1){
  group1 <- subset (pca_coord_rem_real, genotype == gen_1 & day==acq_day)
  group2 <- subset (pca_coord_rem_real, genotype == gen_2 & day==acq_day)
  t_stat = t.test (group1$V1, group2$V1)$statistic  
  return (t_stat)
}

# Function t statistic calculation coordinates of permuted labels 
f_t_stat_only <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=1){
  group1 <- subset (df_coord, genotype == gen_1 & day==acq_day)
  group2 <- subset (df_coord, genotype == gen_2 & day==acq_day)
  t_stat = t.test (group1$V1, group2$V1)$statistic  
  return (t_stat)
}

############################
# Reading data from removal
# It is trickiest because data is splitted in 2 data sets
############################
# ma2=spss.get(path2files)

rem_data = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav", sep=""))
ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters except latency.sav", sep=""))
# tail (ma3)
# Last 5 rows are empty 
ma3 <- head(ma3,-5)

ma3_filt_rem_data <- ma3[ , grepl( "REM" , names( ma3 ) ) ]
ma3_filt_rem_data$ID <- ma3 [ , grepl( "ID" , names( ma3 ) ) ]

genotype_ind <- selected_var_rem$GENTREAT
genotype_ind<- gsub("H20", "", genotype_ind)
genotype_ind <- gsub("NE", "", genotype_ind)

# tail (rem_data)
# tail (ma3_filt_rem_data)

# All variables in rem_data are not present in ma3_filt_rem_data
# We need to add them: NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET
# We keep IDs to perform the joining
rem_data_var = rem_data [ , c(1,7:10)]
rem_data_all_var <- merge(rem_data, ma3_filt_rem_data, all=TRUE)

# Guys starting with 1400277xx are missing in the second table, thus when merging they appear as NA, I delete them
rem_data_all_var <- head (rem_data_all_var, -5)

# It is better to unify labels with acquistion scripts here. 
# Here we found TSNEH20 while in acq we found TS
rem_data_all_var$GENTREAT
genotype_ind<- gsub("H20", "", rem_data_all_var$GENTREAT)
genotype_ind <- gsub("NE", "", genotype_ind)
rem_data_all_var$GENTREAT <- genotype_ind

# head (rem_data_all_var)
# tail (rem_data_all_var)

#########################
# I select for the pca the same variables than in the acquisition data of Ionas
#########################

selected_var_rem <- rem_data_all_var [,c(1:8,10,12,14,15,17,18)] 
head (selected_var_rem)

M.ind = rem_data_all_var[,c(7,8,10,12,14,15,17,18)]

genotype_ind <- selected_var_rem$GENTREAT
# genotype_ind<- gsub("H20", "", genotype_ind)
# genotype_ind <- gsub("NE", "", genotype_ind)


############################
# Median calculation
tbl_stat_median <-with (selected_var_rem, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, LATENCY.TARGET, GALLINDEX.REM, SPEED.REM, PERC.NE.REM, PERC.PERI.REM, WISHAW.REM), list (GENTREAT), FUN=function (x) median=median(x)))                                                            
genotype_tt <- as.factor(tbl_stat_median$Group.1)
genotype_tt_ionas <- as.factor(c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
tbl_stat_median <- tbl_stat_median [,-1]
row.names (tbl_stat_median) <- genotype_tt_ionas
M.med <- tbl_stat_median

jm=rbind(M.med,M.ind)

# res = PCA (jm, scale.unit=TRUE, ind.sup=c(41:455), graph=F) 
res = PCA (jm, scale.unit=TRUE, ind.sup=c(9:91), graph=F) 

pca_coord_rem_real <- as.data.frame(cbind(res$ind.sup$coord[,1],res$ind.sup$coord[,2]))
# pca_coord_rem_real$day <- c(1:415)
pca_coord_rem_real$day <- c(1:83)
pca_coord_rem_real$genotype <- genotype_ind

# I fake a session column in order not to change the function f_t_stat_only
pca_coord_rem_real$day <- "rem"

# pca_coord_rem_real 

real_t_stat <- f_t_stat_only (pca_coord_rem_real, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day="rem")

### ADDING SAMPLING OF THE LABELS
# perm <- 5
# set.seed(333)
# t.values = numeric (perm)
# 
# for (p in 1:perm) {
# 
#   #In removal it is much easier, I don't have to control for different sessions to set always the same animal id 
#   #with the same group along the sessions 
#   p_rem_data_all_var <- rem_data_all_var  
# 
#   p_rem_data_all_var$GENTREAT <- sample (p_rem_data_all_var$GENTREAT)
# 
#   p_selected_var_rem <- p_rem_data_all_var [,c(1:8,10,12,14,15,17,18)] 
#   head (p_selected_var_rem)
#   head(selected_var_rem)
#   
#   # Getting gentreat for PCA
#   p.genotype_ind <- p_selected_var_rem$GENTREAT
#   
#   p.M.ind = p_rem_data_all_var [,c(7,8,10,12,14,15,17,18)]  
# 
#   ############################
#   # Median calculation
#   p_tbl_stat_median <- with (p_selected_var_rem, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, LATENCY.TARGET, GALLINDEX.REM, SPEED.REM, PERC.NE.REM, PERC.PERI.REM, WISHAW.REM), list (GENTREAT), FUN=function (x) median=median(x)))                                                            
#   row.names (p_tbl_stat_median) <- as.factor (p_tbl_stat_median$Group.1)  
#   p_tbl_stat_median <- p_tbl_stat_median [,-1]  
#   p.M.med <- p_tbl_stat_median
# #   M.med
# #   p.M.med
#   p.jm = rbind (p.M.med, p.M.ind)
# 
#   #res = PCA (jm, scale.unit=TRUE, ind.sup=c(41:455), graph=F) 
#   p.res = PCA (p.jm, scale.unit=TRUE, ind.sup=c(9:91), graph=F) 
#   p.pca_coord_rem_real <- as.data.frame (cbind (p.res$ind.sup$coord[,1], p.res$ind.sup$coord[,2]))
#   p.pca_coord_rem_real$day <- c(1:83)
#   p.pca_coord_rem_real$genotype <- p.genotype_ind
#   pca_coord_rem_real
#   # I fake a session column in order not to change the function f_t_stat_only
#   p.pca_coord_rem_real$day <- "rem"
#   
#   # pca_coord_rem_real 
#   pca_coord_rem_real
#   p.pca_coord_rem_real
#   f_t_stat_only (p.pca_coord_rem_real, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day="rem") 
#   t_s <- f_t_stat_only (p.pca_coord_rem_real, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day="rem")  
# #   t_s <- f_t_stat_only (new_coord)
#   t.values[p] <- t_s
# }  



###
# Functions for significance calculation

# Function to calculate significance for a given pairwise comparison
significance_perm_tbl (tbl_1111_rem, gr1="TS", gr2="TSEEEGCG", a_day="rem")  

f_t_stat_only (pca_coord_rem_real, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day="rem") #OK
tbl_perm <-tbl_1111_rem

t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1

significance_perm_tbl <- function (tbl_perm, gr1="TS", gr2="TSEEEGCG", a_day="rem")  {
  real_t_stat <- f_t_stat_only (pca_coord_rem_real, gen_1 = gr1, gen_2 = gr2, acq_day=a_day)
  t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
  t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
  sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 <= real_t_stat])) / length (t_perm_gr1_gr2)
  if (sign_thr >= 0.95) { sign_thr <- 1 - sign_thr }
  
  print (sign_thr)
  return (sign_thr)
} 

# Calling the significance calculation function for each pairwise comparison and generating a table
sign_threshold <- function (tbl_perm, day=5){
  df.t_stats <- c()
  
  # In this case I can not set the order
  gentreat <- unique(c(unique(tbl_perm$V2), unique(tbl_perm$V3)))
  #I have to set the order as well in the r nextflow script
  gentreat_order <- c("WT", "TS","WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")
  
  gentreat_pairs <- t (combn (gentreat_order, 2))
  
  for (row in 1:length(gentreat_pairs [,1])) {
    gr1 <- as.character(gentreat_pairs [row,1])
    gr2 <- as.character(gentreat_pairs [row,2])
    
    sign_thr <- significance_perm_tbl (tbl_perm, gr1, gr2, a_day=day)
    
    v <- c (gr1, gr2, paste (gr1, gr2, sep="_vs_"), sign_thr)
    
    df.t_stats <- rbind (df.t_stats, v)
    colnames (df.t_stats) <- c ("gr1", "gr2", "comp", "sign_thresh")
  } 
  
  df.t_stats <- as.data.frame (df.t_stats, row.names=F, stringsAsFactors = F)
  df.t_stats$sign_thresh <- as.numeric (df.t_stats$sign_thresh)
  df.t_stats <- df.t_stats [order(df.t_stats$sign_thresh),]
  row.names(df.t_stats) <- 1:nrow(df.t_stats)
  return (df.t_stats)
}

tbl_1111_rem <- read.table ("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/PCA_t_statistic_REM1111.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
df.sign_threshold.rem <- sign_threshold (tbl_1111_rem, day="rem")


