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

##Getting HOME directory
home <- Sys.getenv("HOME") 

# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))

rem_data = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav", sep=""))
ma3 = spss.get(paste (home, "/20150515_PCA_old_frotiersPaper/data/Jtracks parameters except latency.sav", sep=""))

# Last 5 rows are empty 
ma3 <- head(ma3,-5)

ma3_filt_rem_data <- ma3[ , grepl( "REM" , names( ma3 ) ) ]
ma3_filt_rem_data$ID <- ma3 [ , grepl( "ID" , names( ma3 ) ) ]

tail (rem_data)
tail (ma3_filt_rem_data)


# All variables in rem_data are not present in ma3_filt_rem_data
# We need to add them: NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET
# We keep IDs to perform the joining
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

########
# ANOVAs
########

# ANOVA with all groups
# data.frame -> rem_data_all_var

rem_nEntries.aov <- aov (NUMBER.ENTRIES ~ genotype  + Error(ID), data = rem_data_all_var)
summary(rem_nEntries.aov)
pairwise.t.test (rem_data_all_var$NUMBER.ENTRIES, rem_data_all_var$genotype, , p.adj="hochberg", paired=F)

rem_latency.aov <- aov (LATENCY.TARGET ~ genotype  + Error(ID), data = rem_data_all_var)
rem_latency.aov <- aov (LATENCY.TARGET ~ genotype,data = rem_data_all_var)
summary (rem_latency.aov) 
rem_nEntries.aov
# Post-hoc
pairwise.t.test (rem_data_all_var$LATENCY.TARGET, rem_data_all_var$genotype, , p.adj="hochberg", paired=F)
TukeyHSD(rem_latency.aov)

# Anova gall index
rem_gallIndex.aov <- aov (GALLINDEX.REM ~ genotype  + Error(ID), data = rem_data_all_var)
summary(rem_gallIndex.aov)

pairwise.t.test (rem_data_all_var$GALLINDEX.REM, rem_data_all_var$genotype, , p.adj="hochberg", paired=F)

# Anova gall index without outlier
rem_data_all_varNoOutlier <- rem_data_all_var [-which (rem_data_all_var$ID == "130054742"),]
pairwise.t.test (rem_data_all_varNoOutlier$GALLINDEX.REM, rem_data_all_varNoOutlier$genotype, , p.adj="hochberg", paired=F)

# Test like juanra without removing the outlier
lm_rem <- lm (GALLINDEX.REM ~ genotype, data = rem_data_all_var)
l2_rem <- glht(lm_rem, linfct = mcp (genotype = "Tukey"))
summary(l2_rem, test = adjusted(type = "BH"))
        
# Like juanra
library(multcomp)
lm_rem_NOoutlier <- lm (GALLINDEX.REM ~ genotype, data = rem_data_all_varNoOutlier)
l2 <- glht(lm_rem_NOoutlier, linfct = mcp (genotype = "Tukey"))
summary(l2, test = adjusted(type = "BH"))

# boxplot number entries
bp_nEntries <- ggplot(rem_data_all_var , aes (genotype, NUMBER.ENTRIES, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "green", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Removal number of entries\n") + xlab ("\ngentreat") + ylab("nEntries\n") +
  theme (legend.title=element_blank())

bp_nEntries + geom_point (position = position_jitter(width = 0.2), colour="red")

# boxplot latency
bp_latency <- ggplot(rem_data_all_var , aes (genotype, LATENCY.TARGET, fill = genotype)) + 
  geom_boxplot() +
  scale_fill_manual(name = "genotype", values = c("red", "green", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Removal latency first entry\n") + xlab ("\ngentreat") + ylab("nEntries\n") +
  theme (legend.title=element_blank())

bp_latency + geom_point (colour="red")
bp_latency + geom_point (position = position_jitter(width = 0.2), colour="red")

rem_data_all_var$NUMBER.ENTRIES
rem_data_all_var_TS <- subset(rem_data_all_var, genotype == "TS" | genotype == "TSEEEGCG" | genotype == "TSEE" | genotype == "TSEGCG")

bp_rem_latency_TS <- ggplot(rem_data_all_var_TS , aes (genotype, LATENCY.TARGET, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("green", "lightblue", "orange", "black")) +
  labs(title = "Removal latency first entry\n") + xlab ("\ngentreat") + ylab("latency\n") +
  theme (legend.title=element_blank())

bp_rem_latency_TS + geom_point (position = position_jitter(width = 0.2), colour="red")

# boxplot gallagher index
bp_gallagher <- ggplot(rem_data_all_var , aes (genotype, GALLINDEX.REM, fill = genotype)) + 
  geom_boxplot (outlier.size=NA, show_guide=FALSE) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Removal gallagher index\n") + xlab ("\ngentreat") + ylab("Gallagher index (cm)\n") +
  theme (legend.title=element_blank()) +
  geom_segment(aes(x = 7.63, y = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLINDEX.REM"]), 
               xend = 8.37, yend = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLINDEX.REM"])), 
               colour="white", size =0.8)

bp_gallagher_points <- bp_gallagher + geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

#PLOT_paper
ggsave (bp_gallagher_points, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "rem_gallagher_all_gr.jpg", sep=""), 
        width=14, height=7, dpi=900)

### Gallagher index without mouse = "130054742" that is considered an outlier according to JR analysis
# Adding mice ID
geom_text (data=neg_positions, aes (x=-Dim.1, y=-Dim.2, label=neg_labels, hjust=1.2), show_guide = FALSE, size=5) + 
outlier <- rem_data_all_var [which (rem_data_all_var$ID == "130054742"),]

# Plot showing the outlier
bp_gallagher_points + geom_text(data=outlier , aes (label=ID))

#Removing the individual from the analysis to make the boxplot
rem_data_all_var_noOutlier <- rem_data_all_var [-which (rem_data_all_var$ID == "130054742"),]
bp_gallagher_noOut <- ggplot(rem_data_all_var_noOutlier, aes (genotype, GALLINDEX.REM, fill = genotype)) + 
  geom_boxplot (outlier.size=NA, show_guide=F) +
  scale_fill_manual(name = "genotype", values = c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black")) +
  labs(title = "Removal gallagher index\n") + xlab ("\ngentreat") + ylab("Gallagher index (cm)\n") +
  theme (legend.title=element_blank()) +
  guides(color=guide_legend(title=NULL)) +
  theme(legend.key = element_blank())+
  geom_segment(aes(x = 7.63, y = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLINDEX.REM"]), 
                   xend = 8.37, yend = median(rem_data_all_var [rem_data_all_var$genotype == "TSEEEGCG","GALLINDEX.REM"])), 
               colour="white", size =0.8) + 
              geom_point (position = position_jitter(width = 0.2), colour="red", show_guide=FALSE)

bp_gallagher_noOut <- bp_gallagher_noOut + geom_point(data=outlier , aes (label=ID), show_guide=FALSE, shape = 2)
bp_gallagher_noOut

#PLOT_paper
ggsave (bp_gallagher_noOut, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "rem_gallagher_all_gr_outlier.jpg", sep=""), 
        width=14, height=7, dpi=900)

# Legend with colours for all group of animals
l <- ggplot() + geom_point(data=rem_data_all_var_noOutlier, aes (x=genotype, y=GALLINDEX.REM, colour = genotype), shape=15, size=5) +
                scale_colour_manual (values=c("red", "darkgreen", "blue", "lightblue", "magenta", "orange", "yellow", "black"))
l <- l + guides(color=guide_legend(title="gentreat")) 
l <- l + theme(legend.key = element_blank())
l

#PLOT_paper
ggsave (l, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "allAnimals_legend.jpg", sep=""), 
        width=14, height=7, dpi=900)

# Only trisomic group
bp_rem_gallagher_TS <- ggplot(rem_data_all_var_TS , aes (genotype, GALLINDEX.REM, fill = genotype)) + 
  geom_boxplot(show_guide=FALSE) +
#   geom_text (aes (label=ID), vjust=-0.5, hjust=1, size=4, show_guide = T)
  geom_text (aes (label=ID))+
  scale_fill_manual(name = "genotype", values = c("green", "lightblue", "orange", "black")) +
  labs(title = "Removal gallagher index\n") + xlab ("\ngentreat") + ylab("Gallagher index (cm)\n") +
  theme (legend.title=element_blank())

bp_rem_gallagher_TS + geom_point (position = position_jitter(width = 0.2), colour="red")







#### 
# Anova with only TS group
rem_data_all_var_TS <- subset(rem_data_all_var, grepl("TS", rem_data_all_var$genotype))
rem_TS_nEntries.aov <- aov(NUMBER.ENTRIES ~ genotype  + Error(ID), data = rem_data_all_var_TS)
summary(rem_TS_nEntries.aov)

bp_ts_rem_nEntries <- ggplot(rem_data_all_var_TS , aes (genotype, NUMBER.ENTRIES, fill = genotype)) + 
  geom_boxplot() +
  scale_fill_manual(name = "genotype", values = c("green","lightblue","orange",  "black")) +
  labs(title = "Removal number of entries target\n") + xlab ("\ngentreat") + ylab("latency\n") +
  theme (legend.title=element_blank())

bp_ts_rem_nEntries + geom_point (colour="red")
bp_ts_rem_nEntries + geom_point (position = position_jitter(width = 0.2), colour="red")

rem_TS_latency.aov <- aov(LATENCY.TARGET ~ genotype  + Error(ID), data = rem_data_all_var_TS)
summary(rem_TS_latency.aov)

# Post-hoc
pairwise.t.test (rem_data_all_var_TS$LATENCY.TARGET, rem_data_all_var_TS$genotype, , p.adj="hochberg", paired=F)

bp_ts_rem_latency <- ggplot(rem_data_all_var_TS , aes (genotype, LATENCY.TARGET, fill = genotype)) + 
  geom_boxplot() +
  scale_fill_manual(name = "genotype", values =c("green","lightblue","orange",  "black")) +
  labs(title = "Removal latency first entry\n") + xlab ("\ngentreat") + ylab("latency\n") +
  theme (legend.title=element_blank())

bp_ts_rem_latency + geom_point (colour="red")
bp_ts_rem_latency + geom_point (position = position_jitter(width = 0.2), colour="red")


#### 
# Anova with only TS group
rem_data_all_var_TS <- subset(rem_data_all_var, grepl("TS", rem_data_all_var$genotype))
rem_TS_nEntries.aov <- aov(NUMBER.ENTRIES ~ genotype  + Error(ID), data = rem_data_all_var_TS)
summary(rem_TS_nEntries.aov)

###########
### Gallagher Index
head(rem_data_all_var_TS)

rem_TS_gallIndex.aov <- aov(GALLINDEX.REM ~ genotype  + Error(ID), data = rem_data_all_var_TS)
summary(rem_TS_gallIndex.aov)

# boxplot gallagher
bp_ts_rem_gallIndex <- ggplot(rem_data_all_var_TS, aes (genotype, GALLINDEX.REM, fill = genotype)) + 
  geom_boxplot() +
  scale_fill_manual(name = "genotype", values =c("green","lightblue","orange",  "black")) +
  labs(title = "Removal latency first entry\n") + xlab ("\ngentreat") + ylab("latency\n") +
  theme (legend.title=element_blank())

bp_ts_rem_gallIndex

pairwise.t.test(rem_data_all_var_TS$GALLINDEX.REM, rem_data_all_var_TS$genotype, , p.adj="hochberg", paired=F)

# Post-hoc
# df.anova$interaction <- paste(df.anova$group, df.anova$time, sep="_")
pairwise.t.test(new_coord_rem$V1, new_coord_rem$genotype, , p.adj="hochberg", paired=F)

