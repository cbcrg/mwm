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

path2files
# ma3 <- read.csv("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/rev_data_f_6v.csv", sep="\t")
ma3 <- read.csv(path2files, sep="\t")


# Setting the order of factor for pairwise comparisons
ma3$gentreat <- factor(ma3$gentreat , levels=c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                               labels=c(c("WT","TS","WTEE","TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")))

# New function taking into account the acquisition days in the comparison
f_t_stat <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG", acq_day=5){
  group1 <- subset (new_coord, genotype == gen_1 & day==acq_day)
  group2 <- subset (new_coord, genotype == gen_2 & day==acq_day)
  t_stat = t.test(group1$V1, group2$V1)$statistic  
  return <- c(t_stat, gen_1, gen_2, paste (gen_1, gen_2, sep="_"))
}


id_group <- subset(ma3, grepl("1", ma3$day))

id_group <- id_group [,c(1,2)]

set.seed (seed)

p_genetreat <- sample(id_group$gentreat)
id_group$gentreat <- p_genetreat

p_ma3 <- merge(id_group, ma3, by="id", all = TRUE)

# I eliminiate the original column with groups labels
# p_ma2 <- p_ma2 [,c(-3)]
# Just in case col order changes better this way
p_ma3 <- subset(p_ma3, select=-c(gentreat.y))

colnames(p_ma3)[2] <- "gentreat"

n_col <- dim (p_ma3)[2]
variables_list <- colnames(p_ma3) [4:n_col] 

tbl_median <- with (p_ma3, aggregate (mget(variables_list), list (gentreat=gentreat, day=day), FUN=median))

tbl_ind <- p_ma3

tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_col_tbl_all <- length(tbl_median[1,])
res = PCA(tbl_med_ind[,(3:n_col_tbl_all)], scale.unit=TRUE, ind.sup=c((n_median+1):length(tbl_med_ind[,1])), graph=F) 

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


# Get all possible combinations of genontype treatment pairwise comparisons 
gentreat <- sort(unique(ma2$gentreat))
gentreat_pairs <- t (combn (gentreat,2))

result <- c()
result_1 <- c()

for (row in 1:length(gentreat_pairs [,1])) {
  gr1 <- as.character(gentreat_pairs [row,1])
  gr2 <- as.character(gentreat_pairs [row,2])
  
  result_v <- c(f_t_stat (new_coord, gr1, gr2, acq_day=3), seed)
  result_v_1 <- c(f_t_stat (new_coord, gr1, gr2, acq_day=1), seed)
  
  result <- rbind  (result, result_v)
  result_1 <- rbind  (result_1, result_v_1)
  
  colnames(result) <- c ("t", "gr1", "gr2", "comparison", "seed") 
  colnames(result_1) <- c ("t", "gr1", "gr2", "comparison", "seed")               
}

wd <- getwd()
write.table(result, file = paste(wd, "/tbl_t_stat_rem3.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)
write.table(result_1, file = paste(wd, "/tbl_t_stat_rem1.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)
