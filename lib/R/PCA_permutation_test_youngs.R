#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Oct 2015                        ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes                 ###
###                                                                       ###
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
young_acq_7var <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv", sep="\t")
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

# new_coord_real_lab <- as.data.frame (res$ind$coord[1:n_median,c(1,2)])
new_coord_real_lab <- as.data.frame(cbind(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2]))
new_coord_real_lab$genotype <- tbl_ind$gentreat
new_coord_real_lab$day <- as.numeric (gsub(".*([0-9]+)$", "\\1", tbl_ind$day))

real_t_stat <- f_t_stat_only (new_coord_real_lab, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=5)
real_t_stat
rm(real_t_stat)

###
# Functions for significance calculation

significance_perm_tbl <- function (tbl_perm, gr1="TS", gr2="TSEEEGCG", n_perm=10001, a_day=5)  {
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
significance_perm_tbl_emp_adjusted <- function (tbl_perm, gr1="TS", gr2="TSEEEGCG", n_perm=10001, a_day=5)  {
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
sign_threshold <- function (tbl_perm, day=5){
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
tbl_1111_day1 <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/tbl/PCA_t_statistic_1111_day1.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
tbl_1111_day5 <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/tbl/PCA_t_statistic_1111_day5.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)

df.sign_threshold.day5 <- sign_threshold (tbl_1111_day5, day=5)
df.sign_threshold.day1 <- sign_threshold (tbl_1111_day1, day=1)

# write.table(df.sign_threshold.day5, file = paste(home, "/20151001_ts65_young_MWM/tbl", "/tbl_sign_thresh_perm_5_revision.csv", sep=""), 
#             sep="\t", row.names=FALSE, dec = ",", col.names=TRUE)

# write.table(df.sign_threshold.day1, file = paste(home, "/20151001_ts65_young_MWM/tbl", "/tbl_sign_thresh_perm_1_revision.csv", sep=""), 
#             sep="\t", row.names=FALSE, dec = ",", col.names=TRUE)


# This part is not used, we do not correct not needed
# How to order make a data frame and order it by significance
# df.sign_threshold.day5 <- as.data.frame(mtx.sign_threshold.day5, row.names=F, stringsAsFactors = F)
# df.sign_threshold.day5$sign_thresh <- as.numeric (df.sign_threshold.day5$sign_thresh)
# # df.sign_threshold.day5$sign_thresh <- as.(df.sign_threshold.day5$sign_thresh)
# df.sign_threshold.day5 <- df.sign_threshold.day5 [order(df.sign_threshold.day5$gr2),]
# df.sign_threshold.day5 <- df.sign_threshold.day5 [order(df.sign_threshold.day5$comp),]
# 
# Ordering for bonferroni holm correction
df.sign_threshold.day5.order <- df.sign_threshold.day5 [order(df.sign_threshold.day5$sign_thresh),]
row.names(df.sign_threshold.day5.order ) <- 1:nrow(df.sign_threshold.day5.order )
df.sign_threshold.day5.order
df.sign_threshold.day5.order$sign_threshBonfHolm <- df.sign_threshold.day5.order$sign_thresh * (length (df.sign_threshold.day5.order$sign_thresh) - as.numeric(row.names(df.sign_threshold.day5.order)) + 1)








# tbl_2222_day1 <- read.table ("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/PCA_t_statistic_2222_day1.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
tbl_2222 <- read.table ("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/PCA_t_statistic_2222.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
head (tbl_2222)
tail (tbl_2222)

################
# Tables session 1

tbl_1111_day1 <- read.table ("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/PCA_t_statistic_1111_day1.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)

mtx.sign_threshold.day1 <- sign_threshold (tbl_1111_day1, day=1)
as.numeric(mtx.sign_threshold.day1 [,4])
as.numeric(mtx.sign_threshold.day1 [,4])/28
0.05/28
tbl_perm <- tbl_1111_day1
head (tbl_1111_day1)
head (tbl_1111)
#df.sign_threshold.session1 <- sign_threshold (tbl_1111_day1, day = 1)
a_day=1
f_t_stat_only (new_coord_real_lab, gen_1 = gr1, gen_2 = gr2, acq_day=a_day)

real_t_stat <- f_t_stat_only (new_coord_real_lab, gen_1 = gr1, gen_2 = gr2, acq_day=a_day)

t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 < real_t_stat])) / length (t_perm_gr1_gr2)
sign_thr

################
gr1<-"WTEEEGCG"
gr2<-"TS"
gr1<-"TS"
gr2<-"TSEEEGCG"
a_day=1
# comparison gives "9.99900009999166e-05" this is like be >10000
tbl_perm <- tbl_1111_day1
real_t_stat <- f_t_stat_only (new_coord_real_lab, gen_1 = gr1, gen_2 = gr2, acq_day=a_day)
t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
head (t_perm_gr1_gr2)
tail (t_perm_gr1_gr2)
sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 < real_t_stat])) / length (t_perm_gr1_gr2)
1-sign_thr

# Function to calculate significance
# t_s_ts_tsee
# tbl_perm <- tbl_3333
gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=1

significance_perm <- function (tbl_perm, gr1="TS", gr2="TSEEEGCG", n_perm=10001, a_day=5)  {
  real_t_stat <- f_t_stat_only (new_coord_real_lab, gen_1 = gr1, gen_2 = gr2, acq_day=a_day)
  t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
  t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
  sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 < real_t_stat])) / length (t_perm_gr1_gr2)
  print (sign_thr)
  
  if (sign_thr >= 0.95) {sign_thr <- 1 - sign_thr}
  
  return (sign_thr)
}

#####
# Calculation of significance using output from cluster
# Output structure
# colnames(result_1) <- c ("t", "gr1", "gr2", "comparison", "seed")  
tbl_1111

# TS vs TSEEEGCG
(28-13+1) * 0.00899910

28*0.00899

0.05/28

write.table(df.sign_threshold.day5.order, file = paste(home, "/20150515_PCA_old_frotiersPaper/tbl", "/tbl_sign_thresh_perm_5.csv", sep=""), 
            sep="\t", row.names=FALSE, dec = ",", col.names=TRUE)

# Get all possible combinations of genontype treatment pairwise comparisons 
gentreat <- unique(ma2$gentreat)
gentreat_pairs <- t (combn (gentreat,2))
gentreat_pairs [1,1]

df.t_stats <- c()

#################
# OLD code in which is based the code above
for (row in 1:length(gentreat_pairs [,1])) {  
  #   v <- c (as.character(gentreat_pairs [row,1]), as.character(gentreat_pairs [row,2]), paste (gentreat_pairs [row,1], gentreat_pairs [row,2], sep="_"), f_t_stat_only (new_coord, gentreat_pairs [row,1], gentreat_pairs [row,2]))
  v <- c (as.character(gentreat_pairs [row,1]), as.character(gentreat_pairs [row,2]), 
          paste (gentreat_pairs [row,1], gentreat_pairs [row,2], sep="_"), 
          significance_perm (tbl_1111, gentreat_pairs [row,1], gentreat_pairs [row,2], a_day=5)) 
  
  df.t_stats <- rbind (df.t_stats, v)
  colnames (df.t_stats) <- c ("gr1", "gr2", "comp", "sign_thresh")
  
}

df.t_stats

f_t_stat_only (new_coord, "TS", "WTEE")
paste(t(combn (gentreat,2))[,1], t(combn (gentreat,2))[,2], sep="_")

f_t_stat_only (new_coord, "TS", "WTEE")

# "TS", "TSEE"
significance_perm (tbl_1111, gr1="TS", gr2="TSEE", a_day=5)
# "TS", "TSEGCG"
significance_perm (tbl_1111, "TS", "TSEGCG", a_day=5)

# "TS", "TSEEEGCG"
significance_perm (tbl_1111, "TS", "TSEEEGCG", a_day=5)

# "TS", "WT"
significance_perm (tbl_1111, "TS", "WT", a_day=5)

# "TSEEEGCG", "TSEGCG"
significance_perm (tbl_1111, "TSEEEGCG", "TSEGCG", a_day=5)

real_t_stat <- f_t_stat_only (new_coord_real_lab, gr1, gr2, acq_day=5)
f_t_stat_only (new_coord_real_lab, "TSEEEGCG", "TSEGCG", acq_day=5)
#"TSEEEGCG", "TSEGCG"
gr1<-"TSEEEGCG"
gr2<-"TSEGCG"

gr1 <- "TS"     
gr2 <- "WTEE"
real_t_stat <- f_t_stat_only (new_coord_real_lab, gr1, gr2, acq_day=5)
tbl_perm<-tbl_1111
t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
# t_perm_gr1_gr2 <- subset (tbl_perm, V3== "TSEGCG")
t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
threshold <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 < real_t_stat])) / length (t_perm_gr1_gr2)
# if (threshold >= 0.95) {threshold <- 1 - threshold} 
threshold
# "TS", "WTEE"
significance_perm (tbl_3333, "TS", "WTEE")

tbl_ts_tsee <- subset (tbl_3333, V4 == "TS_WT")$V1 
t_ts_tsee <-  tbl_ts_tsee$V1
t_ts_tsee <- t_ts_tsee [order(t_ts_tsee)]
t_ts_tsee [9500]
(10000 - length(t_ts_tsee [t_ts_tsee < real_t_s_ts_tsee])) / 10000

t.values [order(t.values)]