#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. July 2015                        ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes REMOVAL         ###
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

# {
#   if (is.null (argsL$path2files)) 
#   {
#     path2files <- "/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/Ts65Dn_OLD_ACQ1_ACQ5_SUBCONJ.sav"
#   }
#   else
#   {
#     path2files <- argsL$path2files
#   }
# }

rem_data = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav", sep=""))
ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks_parameters_except_latency.sav", sep=""))
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

# It is better to unify labels with acquistion scripts here. 
# Here we found TSNEH20 while in acq we found TS
rem_data_all_var$GENTREAT
genotype_ind<- gsub("H20", "", rem_data_all_var$GENTREAT)
genotype_ind <- gsub("NE", "", genotype_ind)
rem_data_all_var$GENTREAT <- genotype_ind

# Setting the order of factor for pairwise comparisons
rem_data_all_var$GENTREAT <- factor(rem_data_all_var$GENTREAT , levels=c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                               labels=c(c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")))

# New function taking into account the acquisition days in the comparison
f_t_stat <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=5){
  group1 <- subset (df_coord, genotype == gen_1 & day==acq_day)
  group2 <- subset (df_coord, genotype == gen_2 & day==acq_day)
  t_stat = t.test(group1$V1, group2$V1)$statistic  
  return <- c(t_stat, gen_1, gen_2, paste (gen_1, gen_2, sep="_"))
}

set.seed (seed)

p_rem_data_all_var <- rem_data_all_var
p_selected_var_rem <- p_rem_data_all_var [,c(1:8,10,12,14,15,17,18)] 
# head (p_selected_var_rem)
# head(selected_var_rem)

# Getting gentreat for PCA
p.genotype_ind <- p_selected_var_rem$GENTREAT

p.M.ind = p_rem_data_all_var [,c(7,8,10,12,14,15,17,18)]  

############################
# Median calculation
p_tbl_stat_median <- with (p_selected_var_rem, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, LATENCY.TARGET, GALLINDEX.REM, SPEED.REM, PERC.NE.REM, PERC.PERI.REM, WISHAW.REM), list (GENTREAT), FUN=function (x) median=median(x)))                                                            
row.names (p_tbl_stat_median) <- as.factor (p_tbl_stat_median$Group.1)  
p_tbl_stat_median <- p_tbl_stat_median [,-1]  
p.M.med <- p_tbl_stat_median
#   M.med
#   p.M.med
p.jm = rbind (p.M.med, p.M.ind)

#res = PCA (jm, scale.unit=TRUE, ind.sup=c(41:455), graph=F) 
p.res = PCA (p.jm, scale.unit=TRUE, ind.sup=c(9:91), graph=F) 
p.pca_coord_rem <- as.data.frame (cbind (p.res$ind.sup$coord[,1], p.res$ind.sup$coord[,2]))
p.pca_coord_rem$day <- c(1:83)
p.pca_coord_rem$genotype <- p.genotype_ind

# I fake a session column in order not to change the function f_t_stat_only
p.pca_coord_rem$day <- "rem"

# Get all possible combinations of genontype treatment pairwise comparisons 
gentreat <- sort(unique(rem_data_all_var$GENTREAT))
gentreat_pairs <- t (combn (gentreat,2))

result <- c()


for (row in 1:length(gentreat_pairs [,1])) {
  gr1 <- as.character(gentreat_pairs [row,1])
  gr2 <- as.character(gentreat_pairs [row,2])
  
  result_v <- c(f_t_stat (p.pca_coord_rem, gr1, gr2, acq_day="rem"), seed)
  
  result <- rbind  (result, result_v)
  
  colnames(result) <- c ("t", "gr1", "gr2", "comparison", "seed") 
}

wd <- getwd()
write.table(result, file = paste(wd, "/tbl_t_stat_REM.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)


