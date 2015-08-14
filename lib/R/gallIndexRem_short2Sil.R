# Calling libraries
library(Hmisc)
library(calibrate)
library(multcomp)


base_size <- 12

dailyInt_theme <- theme_update (
  #axis.text.x = element_text (angle = 90, hjust = 1, size = base_size * 1.5),
  #axis.text.x = element_text (hjust = 1, size = base_size * 1.5),
  axis.text.x = element_text (size = base_size * 1.5),
  axis.text.y = element_text (size = base_size * 1.5),
  axis.title.x = element_text (size=base_size * 1.5, face="bold"),
  axis.title.y = element_text (size=base_size * 1.5, angle = 90, face="bold"),
  #strip.text.x = element_text (size=base_size * 1.3, face="bold"),#facet titles size 
  strip.text.y = element_text (size=base_size * 1.3, face="bold", angle=90),
  plot.title = element_text (size=base_size * 1.5, face="bold"), 
  legend.text = element_text (size=base_size * 1.2),             
  #   panel.grid.major = theme_line (colour = "grey90"),
  panel.grid.major = element_blank(),
  #                   panel.grid.minor = element_blank(), 
  panel.grid.minor = element_blank(),
  #panel.border = element_blank(),
  panel.border = element_rect(colour = "black",fill=NA),
  panel.background = element_blank())

##Getting HOME directory
home <- Sys.getenv("HOME") 

rem_data = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav", sep=""))
ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters except latency.sav", sep=""))

# Last 5 rows are empty 
ma3 <- head(ma3,-5)

ma3_filt_rem_data <- ma3[ , grepl( "REM" , names( ma3 ) ) ]
ma3_filt_rem_data$ID <- ma3 [ , grepl( "ID" , names( ma3 ) ) ]

tail (rem_data)
tail (ma3_filt_rem_data)

rem_data_var = rem_data [ , c(1,7:10)]
rem_data_all_var <- merge(rem_data, ma3_filt_rem_data, all=TRUE)

# Guys starting with 1400277xx are missing in the second table, thus when merging they appear as NA, I delete them
rem_data_all_var <- head (rem_data_all_var, -5)

head (rem_data_all_var)
rem_data_all_var$GENTREAT
rem_data_all_var$genotype <- gsub("H20", "", rem_data_all_var$GENTREAT)
rem_data_all_var$genotype <- gsub("NE", "", rem_data_all_var$genotype)
rem_data_all_var$genotype <- factor(rem_data_all_var$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                                    labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

# boxplot gallagher index
bp_gallagher <- ggplot(rem_data_all_var , aes (genotype, GALLINDEX.REM, fill = genotype)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Removal gallagher index\n") + xlab ("\nGroups") + ylab("Gallagher index (cm)\n") +
  theme (legend.title=element_blank()) +
  #   scale_y_continuous (breaks=1:10) +
  scale_y_continuous(limits=c(18, 110), breaks = c(20,40,60,80,100)) +
  
  geom_segment(aes(x = 7.63, y = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLINDEX.REM"]), 
                   xend = 8.37, yend = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLINDEX.REM"])), 
               colour="white", size =0.8)

bp_gallagher_points <- bp_gallagher + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_gallagher_points
