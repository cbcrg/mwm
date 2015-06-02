#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. May 2015                         ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### PCA analysis only acquisition sessions all individuals No MEDIANS     ### 
#############################################################################

# Calling libraries
library(Hmisc)
library(calibrate)
# install.packages("dplyr")
library("dplyr")

##Getting HOME directory
home <- Sys.getenv("HOME") 

# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))

rem_data = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav", sep=""))
ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters except latency.sav", sep=""))

# Last 5 rows are empty 
ma3 <- head(ma3,-5)
tail (ma3)


# Setting ID as the name of the columnn for MICE.ID, otherwise dataframes are not correctly joined
ma3$ID <- ma3 [ , grepl( "ID" , names( ma3 ) ) ]

tail (ma3)
tail (rem_data)

# All variables in rem_data are not present in ma3_filt_rem_data
# We need to add them: NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET
# We keep IDs to perform the joining
rem_data_var = rem_data [ , c(1,7:10)]

# Guys starting with 1400277xx are missing in the second table, thus when merging they appear as NA, I delete them
rem_data_var <- head (rem_data_var, -5)

data_all_sessions <- merge(rem_data_var, ma3, all=TRUE)
tail(data_all_sessions)

# I eliminate all columns but GEN.TREAT and variables
data_all_sessions_var <- data_all_sessions [ , c(9:length(data_all_sessions[1,]))]
tail (data_all_sessions_var)
tbl_median <- aggregate(.~GEN.TREAT,data_all_sessions_var, median)

# another way of doing the same
# median_tbl <- data_all_sessions_var %>% group_by(GEN.TREAT) %>% summarise_each(funs(median))
# tail (median_tbl)

genotype_tt_ionas <- as.factor (tbl_median$GEN.TREAT)
tbl_median <- tbl_median [,-1]
PCA_all <- prcomp (tbl_median, scale=TRUE)
summary(PCA_all)

# Plot PCA color by genotype
pca_all_2plot <- as.data.frame (PCA_all$x)

PCA_all
g_genotype_tt <- ggplot(pca_all_2plot, aes(PC1, PC2)) + geom_point(aes(colour=genotype_tt_ionas), size=4) +                                                           
  labs (title = "PCA") +
  #   scale_color_manual(values=cols, labels=c("1", "2", "3")) +
  #   geom_text (aes (label=genotype_tt), hjust=0, vjust=-0.5)
  geom_text (aes (label=genotype_tt_ionas), hjust=0.5, vjust=-0.5) +
  xlim (c(-10, 10)) + ylim (c(-5,5)) +
  labs (title = "PCA all sessions") +  
  labs (x = "\nPC1", y="PC2\n") +
  theme (legend.key=element_rect(fill=NA), legend.title=element_blank())
#, position=position_jitter(h=0), alpha = 1)

g_genotype_tt







genotype_tt_ionas <- as.factor(tbl_median$GEN.TREAT)
tbl_median <- tbl_median [,-1]
PCA_all <- prcomp(tbl_median, scale=TRUE)
summary(PCA_all)

# Plot PCA color by genotype
pca_all_2plot <- as.data.frame (PCA_all$x)
# pca_all_2plot$PC1_neg <- -pca_rem_2plot$PC1
g_genotype_tt <- ggplot(pca_all_2plot, aes(PC1, PC2)) + geom_point(aes(colour=genotype_tt_ionas), size=4) +                                                           
  labs (title = "PCA") +
  #   scale_color_manual(values=cols, labels=c("1", "2", "3")) +
  #   geom_text (aes (label=genotype_tt), hjust=0, vjust=-0.5)
  geom_text (aes (label=genotype_tt_ionas), hjust=0.5, vjust=-0.5) +
  xlim (c(-10, 10)) + ylim (c(-5,5)) +
  labs (title = "PCA all sessions") +  
  labs (x = "\nPC1", y="PC2\n") +
  theme (legend.key=element_rect(fill=NA), legend.title=element_blank())
#, position=position_jitter(h=0), alpha = 1)

g_genotype_tt

