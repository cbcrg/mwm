#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. June 2015                        ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### Removal session checking for single variables significance:           ###
### N of entries                                                          ###
### Latency to 1st entry to platform area                                 ### 
#############################################################################

# Calling libraries
library(Hmisc)
library(calibrate)
library(multcomp)
library(ggplot2)

##Getting HOME directory
home <- Sys.getenv("HOME") 

# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))

########################
# Tables before revision
# rem_data = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav", sep=""))
# ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters except latency.sav", sep=""))
# 
# Last 5 rows are empty 
# ma3 <- head(ma3,-5)
#######################

# Table after revision
ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters_all_including_extra_3mice.sav", sep=""))
head (ma3)
ma3

ma3_filt_rem_data <- ma3[ , grepl( "REM" , names( ma3 ) ) ]
ma3_filt_rem_data$ID <- ma3 [ , grepl( "ID" , names( ma3 ) ) ]

# tail (rem_data)
tail (ma3_filt_rem_data)

rem_data_all_var <- ma3_filt_rem_data
rem_data_all_var$genotype <- ma3$GEN.TREAT
rem_data_all_var$genotype <- factor(rem_data_all_var$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                                    labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

##################
# Gallagher index
##################
### ANOVA
rem_gallIndex.aov <- aov (GALLAGHER.INDEX.REM ~ genotype  + Error(ID), data = rem_data_all_var)
summary(rem_gallIndex.aov)

# Si elimino los tres ultimos animales a??adidos vuelvo a tener el resultado anterior 
# rem_data_all_var <- head(rem_data_all_var,-3)

######
# Like juanra
# Test like juanra without removing the outlier
lm_rem <- lm (GALLAGHER.INDEX.REM ~ genotype, data = rem_data_all_var)
l2_rem <- glht(lm_rem, linfct = mcp (genotype = "Tukey"))
summary(l2_rem, test = adjusted(type = "BH"))

### Boxplot
bp_gallagher <- ggplot(rem_data_all_var , aes (genotype, GALLAGHER.INDEX.REM, fill = genotype)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Removal gallagher index\n") + xlab ("\nGroups") + ylab("Gallagher index (cm)\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(18, 110), breaks = c(20,40,60,80,100)) +
  
  geom_segment(aes(x = 7.63, y = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLAGHER.INDEX.REM"]), 
                   xend = 8.37, yend = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLAGHER.INDEX.REM"])), 
               colour="white", size =0.8)

bp_gallagher_points <- bp_gallagher + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_gallagher_points

#PLOT_paper
ggsave (bp_gallagher_points, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig3_rem_singleVar/", "rem_gallagher_all_gr.jpg", sep=""), 
         width=14, height=7, dpi=900)


###################################
# Time spent in the target quadrant
###################################
### ANOVA
rem_perc_NE.aov <- aov (PERCENT.NE.REM ~ genotype  + Error(ID), data = rem_data_all_var)
summary(rem_perc_NE.aov)

# Si elimino los tres ultimos animales a??adidos vuelvo a tener el resultado anterior 
# rem_data_all_var <- head(rem_data_all_var,-3)

######
# Like juanra
# Test like juanra without removing the outlier
lm_rem <- lm (PERCENT.NE.REM ~ genotype, data = rem_data_all_var)
l2_rem <- glht(lm_rem, linfct = mcp (genotype = "Tukey"))
summary(l2_rem, test = adjusted(type = "BH"))

### Boxplot
bp_perc_NE <- ggplot(rem_data_all_var , aes (genotype, PERCENT.NE.REM, fill = genotype)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Removal percentage of time in target quadrant\n") + xlab ("\nGroups") + ylab("% target quadrant\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(0, 80), breaks = c(0,20,40,60,80)) +
  
  geom_segment(aes(x = 7.63, y = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","PERCENT.NE.REM"]), 
                   xend = 8.37, yend = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","PERCENT.NE.REM"])), 
               colour="white", size =0.8)

bp_perc_NE_points <- bp_perc_NE + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_perc_NE_points

#PLOT_paper
ggsave (bp_perc_NE_points, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig3_rem_singleVar/", "rem_perc_NE_all_gr.jpg", sep=""), 
        width=14, height=7, dpi=900)















