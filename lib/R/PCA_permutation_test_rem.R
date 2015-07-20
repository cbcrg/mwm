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

# Function t statistic calculation
f_t_stat_only <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=1){
  group1 <- subset (pca_coord_rem_real, genotype == gen_1 & day==acq_day)
  group2 <- subset (pca_coord_rem_real, genotype == gen_2 & day==acq_day)
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

# tail (rem_data)
# tail (ma3_filt_rem_data)

# All variables in rem_data are not present in ma3_filt_rem_data
# We need to add them: NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET
# We keep IDs to perform the joining
rem_data_var = rem_data [ , c(1,7:10)]
rem_data_all_var <- merge(rem_data, ma3_filt_rem_data, all=TRUE)

# Guys starting with 1400277xx are missing in the second table, thus when merging they appear as NA, I delete them
rem_data_all_var <- head (rem_data_all_var, -5)

# head (rem_data_all_var)
# tail (rem_data_all_var)
#########################
# I select for the pca the same variables than in the acquisition data of Ionas
#########################

selected_var_rem <- rem_data_all_var [,c(1:8,10,12,14,15,17,18)] 
head (selected_var_rem)

M.ind = rem_data_all_var[,c(7,8,10,12,14,15,17,18)]

genotype_ind <- selected_var_rem$GENTREAT
genotype_ind<- gsub("H20", "", genotype_ind)
genotype_ind <- gsub("NE", "", genotype_ind)


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

f_t_stat_only (pca_coord_rem_real, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day="rem")