#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. March 2016                       ###
#############################################################################
### MWM young ts65dn                                                      ###
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
# ma2 <- read.csv("/Users/jespinosa/20151001_ts65_young_MWM/data/ts65_young.csv", sep="\t")#del

ma2 <- read.csv(path2files, sep="\t")

# Compares groups annotated by the interaction bw gen_tt and day
f_t_stat <- function (df_coord, gen_tt_day_a = "TS_1", gen_tt_day_b = "TS_2"){
  group1 <- subset (df_coord, gentreat_day == gen_tt_day_a)
  group2 <- subset (df_coord, gentreat_day == gen_tt_day_b)
  t_stat = t.test(group1$V1, group2$V1)$statistic  
  return <- c(t_stat, gen_tt_day_a, gen_tt_day_b, paste (gen_tt_day_a, gen_tt_day_b, sep="_"))
}

# gsub(ma2$day, "Day ", "")
ma2$day_n <- gsub ("Day ", "", ma2$day)

## Group labels
id_group <- subset(ma2, grepl("1", ma2$day))

# Lo mas facil es hacer una variable combinada de genotipo treatment
ma2$gentreat_day <- paste(ma2$gentreat, ma2$day_n, sep="_")

# id_group <- id_group [,c(1,2)]
id_gentreat_day <- subset (ma2, select=c(id, gentreat_day))

# seed<-1111
# print (paste("ddddddddd", seed))
set.seed (seed)

# p_genetreat <- sample(id_group$gentreat)
# id_group$gentreat <- p_genetreat
# p_ma2 <- merge(id_group, ma2, by="id", all = TRUE)

p_genetreat_day <- sample(id_gentreat_day$gentreat_day)
id_gentreat_day$gentreat_day <- p_genetreat_day
p_ma2 <- merge(id_gentreat_day, ma2, by="id", all = TRUE)


## gentreat_day.x has the randomize group
## I eliminiate the original column with groups labels
## now I also delete day_n, gentreat and day
p_ma2 <- subset(p_ma2, select=-c(gentreat_day.y, day_n, gentreat, day))

colnames(p_ma2)[2] <- "gentreat_day"

n_col <- dim (p_ma2)[2]
variables_list <- colnames(p_ma2) [3:n_col] 

tbl_median <- with (p_ma2, aggregate (mget(variables_list), list (gentreat_day=gentreat_day), FUN=median))

tbl_ind <- p_ma2

tbl_med_ind <- rbind (tbl_median, tbl_ind[,-1])
n_median <- length(tbl_median[,1])
n_col_tbl_all <- length(tbl_median[1,])
res = PCA(tbl_med_ind[,(3:n_col_tbl_all)], scale.unit=TRUE, ind.sup=c((n_median+1):length(tbl_med_ind[,1])), graph=F) 

# Assigning again the labels of the groups for t statistic calculation
new_coord <- as.data.frame(cbind(res$ind.sup$coord[,1],res$ind.sup$coord[,2]))
new_coord$gentreat_day <- tbl_ind$gentreat_day

new_coord$gentreat_day <- factor(new_coord$gentreat_day , levels=c("WT_1", "WT_2", "WT_3", "WT_4", "WT_5", 
                                                           "TS_1", "TS_2", "TS_3", "TS_4", "TS_5", 
                                                           "WTEEEGCG_1", "WTEEEGCG_2", "WTEEEGCG_3", "WTEEEGCG_4", "WTEEEGCG_5",
                                                           "TSEEEGCG_1", "TSEEEGCG_2", "TSEEEGCG_3", "TSEEEGCG_4", "TSEEEGCG_5"))                     


# Get all possible combinations of genontype treatment pairwise comparisons 
gentreat_day <- sort(unique(ma2$gentreat_day))
gentreat_day_pairs <- t (combn (gentreat_day,2))

result <- c()
result_1 <- c()

for (row in 1:length(gentreat_day_pairs [,1])) {
  gr1 <- as.character(gentreat_day_pairs [row,1])
  gr2 <- as.character(gentreat_day_pairs [row,2])
  
#   result_v <- c(f_t_stat (new_coord, gr1, gr2, acq_day=5), seed)
  result_v <- c(f_t_stat(new_coord, gen_tt_day_a = gr1, gen_tt_day_b = gr2), seed)
  # Before I was comparing only day 1 and day 5 but now I want to compare the interactions
#   result_v_1 <- c(f_t_stat (new_coord, gr1, gr2, acq_day=1), seed)

  result <- rbind  (result, result_v)
#   result_1 <- rbind  (result_1, result_v_1)
  
  colnames(result) <- c ("t", "gr1", "gr2", "comparison", "seed") 
#   colnames(result_1) <- c ("t", "gr1", "gr2", "comparison", "seed")               
}

group1 <- subset (new_coord, gentreat_day == "TS_1")
group2 <- subset (new_coord, gentreat_day == "TS_2")
t_stat = t.test(group1$V1, group2$V1)$statistic

wd <- getwd()
write.table(result, file = paste(wd, "/tbl_t_stat_interaction.csv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE)
