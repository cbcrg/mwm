#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Oct 2015                         ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes                 ###
### REVERSAL DATA                                                         ###
#############################################################################

library(FactoMineR)
library(Hmisc)
home <- Sys.getenv("HOME")

# Function t statistic calculation
f_t_stat_only <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=1){
  group1 <- subset (new_coord_real_lab, genotype == gen_1 & day==acq_day)
  group2 <- subset (new_coord_real_lab, genotype == gen_2 & day==acq_day)
  t_stat = t.test (group1$V1, group2$V1)$statistic  
  return (t_stat)
}

# young data already adapted from PCA_MWM_youngTs65.R
# we load the table from there
ma3 <- read.csv("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/rev_data_f_6v.csv", sep="\t")

# Setting the order of factor for pairwise comparisons
ma3$gentreat <- factor(ma3$gentreat , levels=c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                       labels=c(c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")))

rev_data_f_6v <- ma3
n_col <- dim (rev_data_f_6v)[2]
variables_list <- colnames(rev_data_f_6v) [4:n_col] 

tbl_median <- with (rev_data_f_6v, aggregate (cbind (distance, gallindex, speed, percentsw, perperi, wishaw, latency), 
                                              list (day=day, gentreat=gentreat), FUN=median))

rownames (tbl_median) <- paste (tbl_median[,2], gsub ("Day ", "", tbl_median[,1]), sep="")

tbl_ind <- rev_data_f_6v
head(tbl_median)
head(tbl_ind[,-1])
tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_median_plus1 <- n_median + 1

res = PCA(tbl_med_ind[,(3:8)], scale.unit=TRUE, ind.sup=c(n_median_plus1:length(tbl_med_ind[,1]))) 

# new_coord_real_lab <- as.data.frame (res$ind$coord[1:n_median,c(1,2)])
new_coord_real_lab <- as.data.frame(cbind(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2]))
new_coord_real_lab$genotype <- tbl_ind$gentreat
new_coord_real_lab$day <- as.numeric (gsub(".*([0-9]+)$", "\\1", tbl_ind$day))

real_t_stat <- f_t_stat_only (new_coord_real_lab, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=3)
real_t_stat
rm(real_t_stat)

###
# Functions for significance calculation

significance_perm_tbl <- function (tbl_perm, gr1="TS", gr2="TSEEEGCG", n_perm=10001, a_day=3)  {
  real_t_stat <- f_t_stat_only (new_coord_real_lab, gen_1 = gr1, gen_2 = gr2, acq_day=a_day)
  t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
  t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
  sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 <= real_t_stat])) / length (t_perm_gr1_gr2)
  if (sign_thr >= 0.95) { sign_thr <- 1 - sign_thr }
  
  print (sign_thr)
  return (sign_thr)
}

# # Calculation of what I though it whas the empirical adjusted p value
# # It is not correct
significance_perm_tbl_emp_adjusted <- function (tbl_perm, gr1="TS", gr2="TSEEEGCG", n_perm=10001, a_day=3)  {
  real_t_stat <- f_t_stat_only (new_coord_real_lab, gen_1 = gr1, gen_2 = gr2, acq_day=a_day)
  t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
  t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
  sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 <= real_t_stat])) / length (t_perm_gr1_gr2)
  
  sign_thr_adj <- 0
  
  if (sign_thr >= 0.95) { 
    sign_thr_adj <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 >= real_t_stat]) + 1) / (length (t_perm_gr1_gr2) + 1) 
    print (paste("tail mod",sign_thr_adj))
  }
  else {
    sign_thr_adj <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 <= real_t_stat]) + 1) / (length (t_perm_gr1_gr2) + 1)
    print (paste("tail",sign_thr_adj))
  }
  
  #   print (sign_thr)
  return (sign_thr_adj)
}

# significance_perm_tbl_emp_adjusted (tbl_1111, a_day = 5)

# Calling the significance calculation function for each 
# pairwise comparison and generating a table
sign_threshold <- function (tbl_perm, day=3){
  df.t_stats <- c()
  
  # In this way the order is given by the input table, easier
  comparisons <- unique(tbl_perm$V4)
  gentreat_pairs <- data.frame(do.call('rbind', strsplit(as.character(comparisons),'_',fixed=TRUE)))
  
  for (row in 1:length(gentreat_pairs [,1])) {
    gr1 <- as.character(gentreat_pairs [row,1])
    gr2 <- as.character(gentreat_pairs [row,2])
    print (paste("===============",gr1, gr2, sep=" "))
    sign_thr <- significance_perm_tbl (tbl_perm, gr1, gr2, a_day=day)
    adj_sign_thr <- significance_perm_tbl_emp_adjusted (tbl_perm, gr1, gr2, a_day=day)
    
    v <- c (gr1, gr2, paste (gr1, gr2, sep="_vs_"), sign_thr, adj_sign_thr)
    
    df.t_stats <- rbind (df.t_stats, v)
    colnames (df.t_stats) <- c ("gr1", "gr2", "comp", "sign_thresh", "adj_sign_thr")
  } 
  
  df.t_stats <- as.data.frame (df.t_stats, row.names=F, stringsAsFactors = F)
  df.t_stats$sign_thresh <- as.numeric (df.t_stats$sign_thresh)
  df.t_stats$adj_sign_thr <- as.numeric (df.t_stats$adj_sign_thr)
  df.t_stats <- df.t_stats [order(df.t_stats$sign_thresh),]
  row.names(df.t_stats) <- 1:nrow(df.t_stats)
  return (df.t_stats)
}

###################################
# Reading results from 10000 permutations 3X for all comparisons
tbl_1111_day1 <- read.table ("/Users/jespinosa/20150515_PCA_old_frotiersPaper/tbl/PCA_t_statistic_reversal_1111_day1.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
tbl_1111_day3 <- read.table ("/Users/jespinosa/20150515_PCA_old_frotiersPaper/tbl/PCA_t_statistic_reversal_1111_day3.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)

df.sign_threshold.day3 <- sign_threshold (tbl_1111_day3, day=3)
df.sign_threshold.day1 <- sign_threshold (tbl_1111_day1, day=1)

write.table(df.sign_threshold.day3, file = paste(home, "/20150515_PCA_old_frotiersPaper/tbl", "/tbl_sign_thresh_perm_3_reversal.csv", sep=""), 
            sep="\t", row.names=FALSE, dec = ",", col.names=TRUE)

write.table(df.sign_threshold.day1, file = paste(home, "/20150515_PCA_old_frotiersPaper/tbl", "/tbl_sign_thresh_perm_1_reversal.csv", sep=""), 
            sep="\t", row.names=FALSE, dec = ",", col.names=TRUE)
