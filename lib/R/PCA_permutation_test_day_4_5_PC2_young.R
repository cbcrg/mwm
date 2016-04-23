#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. March 2016                       ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes                 ###
### Comparisons bw genotypes untreated vs treated                         ###
#############################################################################

library(FactoMineR)
library(Hmisc)
home <- Sys.getenv("HOME")

# Compares groups annotated by the interaction bw gen_tt and day
# For PC2!!! V2
f_t_stat <- function (df_coord, gen_tt_day_a = "TS_1", gen_tt_day_b = "TS_2"){
  group1 <- subset (df_coord, gentreat_day == gen_tt_day_a)
  group2 <- subset (df_coord, gentreat_day == gen_tt_day_b)
  t_stat = t.test(group1$V2, group2$V2)$statistic  
  return (t_stat)
}

# young data already adapted from PCA_MWM_youngTs65.R
# we load the table from there
# ma2 <- read.csv("/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv", sep="\t")#del
young_acq_7var <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv", sep="\t")

# young_acq_7var<- ma2
head (young_acq_7var)
# young_acq_7var$gentreat <- factor(young_acq_7var$gentreat , levels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"), 
#                                   labels=c("WT", "TS", "WTEEEGCG", "TSEEEGCG"))

tbl_median <- with (young_acq_7var, aggregate (cbind (distance, gallindex, latency, speed, percentne, percentperi, whishaw), 
                                               list (gentreat=gentreat, day=day), FUN=median))

rownames (tbl_median) <- paste (tbl_median[,1], gsub ("Day ", "", tbl_median[,2]), sep="")
tbl_ind <- young_acq_7var
tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_median_plus1 <- n_median + 1

res = PCA(tbl_med_ind[,(3:9)], scale.unit=TRUE, ind.sup=c(n_median_plus1:length(tbl_med_ind[,1]))) 

new_coord_real_lab <- as.data.frame(cbind(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2]))
new_coord_real_lab$genotype <- tbl_ind$gentreat
new_coord_real_lab$day <- as.numeric (gsub(".*([0-9]+)$", "\\1", tbl_ind$day))

## Generate a variable gentreat_day with gentreat and day together
new_coord_real_lab$gentreat_day <- paste(new_coord_real_lab$genotype, new_coord_real_lab$day, sep="_")

real_t_stat <- f_t_stat(new_coord_real_lab, gen_tt_day_a = "TS_4", gen_tt_day_b = "TS_5")
f_t_stat(new_coord_real_lab, gen_tt_day_a = "TSEEEGCG_4", gen_tt_day_b = "TSEEEGCG_5")
real_t_stat
rm(real_t_stat)

###
# Functions for significance calculation
significance_perm_tbl <- function (tbl_perm, gr1="TS_4", gr2="TS_5", n_perm=10001)  {  

  real_t_stat <- f_t_stat (new_coord_real_lab, gen_tt_day_a = gr1, gen_tt_day_b = gr2)
  t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
  t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
  sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 <= real_t_stat])) / length (t_perm_gr1_gr2)
#   if (sign_thr >= 0.95) { sign_thr <- 1 - sign_thr }
  sign_thr <- 1 - sign_thr
  return (list(sign_thr, real_t_stat))
}

# significance_perm_tbl_emp_adjusted (tbl_1111, a_day = 5)

# Calling the significance calculation function for each 
# pairwise comparison and generating a table
# sign_threshold <- function (tbl_perm, day=5){
sign_threshold <- function (tbl_perm){
  df.t_stats <- c()
  
  # In this way the order is given by the input table, easier
  comparisons <- paste(unique (tbl_perm$V2), unique (tbl_perm$V3), sep="=")

  gentreat_pairs <- data.frame(do.call('rbind', strsplit(as.character(comparisons),'=',fixed=TRUE)))

  for (row in 1:length(gentreat_pairs [,1])) {
    
    gr1 <- as.character(gentreat_pairs [row,1])
    gr2 <- as.character(gentreat_pairs [row,2])
    
    list_sign_thr_real_t <- significance_perm_tbl (tbl_perm, gr1, gr2) 
    sign_thr <- list_sign_thr_real_t[[1]]
    real_t <- list_sign_thr_real_t[[2]]

    v <- c (gr1, gr2, paste (gr1, gr2, sep="_vs_"), sign_thr, real_t)

    df.t_stats <- rbind (df.t_stats, v)
    colnames (df.t_stats) <- c ("gr1", "gr2", "comp", "p_value", "pseudo_t")
  } 
  
  df.t_stats <- as.data.frame (df.t_stats, row.names=F, stringsAsFactors = F)
  df.t_stats$sign_thresh <- as.numeric (df.t_stats$p_value)
  df.t_stats <- df.t_stats [order(df.t_stats$p_value),]
  row.names(df.t_stats) <- 1:nrow(df.t_stats)

  return (df.t_stats)
}

###################################
# Reading results from 10000 permutations 3X for all comparisons

tbl_1111_day4_vs_5 <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/tbl/PCA_t_statistic_PC2_day4_5_1111.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
tbl_2222_day4_vs_5 <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/tbl/PCA_t_statistic_PC2_day4_5_2222.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
tbl_interaction_day4_vs_5 <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/tbl/PCA_t_statistic_interaction_1111.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)
tbl_1111_randomized_days <- read.table ("/Users/jespinosa/20151001_ts65_young_MWM/tbl/PCA_t_statistic_PC2_randomizedDays_1111.csv", sep="\t", dec=".", header=F, stringsAsFactors=F)

df.sign_threshold_day4_vs_5 <- sign_threshold (tbl_1111_day4_vs_5)
df.sign_threshold_day4_vs_5 <- sign_threshold (tbl_2222_day4_vs_5)
write.table(df.sign_threshold_day4_vs_5, file = paste(home, "/20151001_ts65_young_MWM/tbl", "/tbl_sign_thresh_perm_PC2_4vs5.csv", sep=""), 
            sep="\t", row.names=FALSE, dec = ".", col.names=TRUE)

df.sign_threshold_interaction_day4_vs_5 <- sign_threshold (tbl_interaction_day4_vs_5)

df.sign_threshold_randomized_days <- sign_threshold (tbl_1111_randomized_days)
write.table(df.sign_threshold_randomized_days, file = paste(home, "/20151001_ts65_young_MWM/tbl", "/tbl_randomized_days.csv", sep=""), 
            sep="\t", row.names=FALSE, dec = ".", col.names=TRUE)



# tbl_1111_day4_vs_5 [which (tbl_1111_day4_vs_5$V5==1111),]
# gr1="TSEEEGCG_4"
# gr2="TSEEEGCG_5"
# tbl_perm <- tbl_1111_day4_vs_5
# head (tbl_perm)
# t_perm_gr1_gr2 <- subset (tbl_perm, V4 == paste(gr1, gr2, sep="_"))$V1
# min(t_perm_gr1_gr2)
# real_t_stat = -0.3176
# t_perm_gr1_gr2 <- t_perm_gr1_gr2 [order(t_perm_gr1_gr2)]
# head (t_perm_gr1_gr2)
# length(t_perm_gr1_gr2)
# sign_thr <- (length (t_perm_gr1_gr2) - length(t_perm_gr1_gr2 [t_perm_gr1_gr2 <= real_t_stat])) / length (t_perm_gr1_gr2)
# 1-sign_thr
# t_perm_gr1_gr2[469]
# new_coord_real_lab
# group1 <- subset (new_coord_real_lab, gentreat_day == "WT_4")
# group2 <- subset (new_coord_real_lab, gentreat_day == "WT_5")
# summary(group1$V2)
# summary(group2$V2)
# length(group1$V2)
# length(group2$V2)
# 
# t.test(group1$V2, group2$V2)
# summary(new_coord_real_lab)
# 
# group1 <- subset (new_coord_real_lab, gentreat_day == "TS_4")
# group2 <- subset (new_coord_real_lab, gentreat_day == "TS_5")
# summary(group1$V2)
# summary(group2$V2)
# length(group1$V2)
# length(group2$V2)
# 
# t.test(group1$V2, group2$V2)
# summary(new_coord_real_lab)
# f_t_stat <- function (df_coord, gen_tt_day_a = "TS_1", gen_tt_day_b = "TS_2"){
#   group1 <- subset (df_coord, gentreat_day == gen_tt_day_a)
#   group2 <- subset (df_coord, gentreat_day == gen_tt_day_b)
#   t_stat = t.test(group1$V2, group2$V2)$statistic  
#   return (t_stat)
# }
# 
# 
# t.test(y1,y2)
# 
# if (sign_thr >= 0.95) { sign_thr <- 1 - sign_thr }
# 


