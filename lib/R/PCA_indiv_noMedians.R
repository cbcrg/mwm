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

# Reading latency, the structure of the table is not the correct, days are in rows
ma2=spss.get("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/Ts65Dn OLD ACQ1_ACQ5_SUBCONJ.sav")
tail (ma2)

# filter by day
lat_1 <- subset(ma2, day=="Day 1", c("id", "latency"))
lat_2 <- subset(ma2, day=="Day 2", c("id", "latency"))
lat_3 <- subset(ma2, day=="Day 3", c("id", "latency"))
lat_4 <- subset(ma2, day=="Day 4", c("id", "latency"))
lat_5 <- subset(ma2, day=="Day 5", c("id", "latency"))

acq_latency <- c()
acq_latency <- merge(lat_1, lat_2, by=c("id"))
colnames(acq_latency) <- c ("id", "lat.ACQ1", "lat.ACQ2")
acq_latency <- merge(acq_latency, lat_3, by=c("id"))
acq_latency <- merge(acq_latency, lat_4, by=c("id"))
acq_latency <- merge(acq_latency, lat_5, by=c("id"))
colnames(acq_latency) <- c ("id", "lat.ACQ1", "lat.ACQ2", "lat.ACQ3", "lat.ACQ4", "lat.ACQ5")
head (data_all_sessions_var)

# I eliminate all columns but GEN.TREAT and variables
df.ACQ_data_except_latency <- data_all_sessions [ , grepl( "ACQ" , names( data_all_sessions ) ) ]
head(df.ACQ_data_except_latency)
# I add the ids to the table
id_and_groups<- data_all_sessions[,c(6:9)]
colnames(id_and_groups) [1] <- "id"
df.ACQ_data_except_latency <- cbind(id_and_groups,df.ACQ_data_except_latency)
head (df.ACQ_data_except_latency)

df.ACQ_data <- merge (df.ACQ_data_except_latency, acq_latency, by=c("id"))

# Perform a PCA of all the variables and inviduals during acquisition sessions
# filtering group labels
df.ACQ_data_no_groups <- df.ACQ_data [ , c(5:length(df.ACQ_data[1,]))]

# PCA individuals
PCA_individuals_acq <- prcomp (df.ACQ_data_no_groups, scale=TRUE)
summary(PCA_individuals_acq)

# Plot
genotype_tt <- as.factor (df.ACQ_data$GEN.TREAT)
genotype_tt <- factor(genotype_tt , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                      labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_indiv_2plot <- as.data.frame (PCA_individuals_acq$x)

p_genotype_tt <- ggplot(pca_indiv_2plot, aes(PC1, PC2)) + geom_point(aes(colour=genotype_tt), size=4) +                                                           
  #   scale_color_manual(values=cols, labels=c("1", "2", "3")) +
  #   geom_text (aes (label=genotype_tt), hjust=0, vjust=-0.5)
#   geom_text (aes (label=genotype_tt), hjust=0.5, vjust=-0.5) +
  scale_color_manual(values=c("red", "green", "blue", "lightblue", 
                              "magenta", "orange", "gray", "black")) +
  xlim (c(-10, 10)) + ylim (c(-5,5)) +
  labs (title = "PCA individuals") +  
  labs (x = "\nPC1", y="PC2\n") +
  theme (legend.key=element_rect(fill=NA), legend.title=element_blank())
#, position=position_jitter(h=0), alpha = 1)

p_genotype_tt

p_genotype_tt <- ggplot(pca_indiv_2plot, aes(PC1, PC2, color=genotype_tt, label=rownames(pca_indiv_2plot))) + 
                        stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                                       geom="polygon", color=NA, n=200, h=4, bins=6) + 
                        geom_smooth(se=F, method='lm') + 
                        geom_point() + 
                        scale_color_manual(name='genotype_tt', 
                                           values = c("red", "green", "blue", "lightblue", 
                                                      "magenta", "orange", "yellow", "black"), 
                                           labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
                        scale_fill_manual( name='mutation', 
                                           values = c("red", "green", "blue", "lightblue", 
                                                      "magenta", "orange", "yellow", "black"), 
                                           labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
                        geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
                        scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
                        scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
                        coord_cartesian(xlim=c(-14, 10),
                                        ylim=c(-8, 8))                    
#     coord_cartesian(xlim=c(-10, 10)) + ylim=c(-5,5))
                        
# 
p_genotype_tt

# Statistical comparison of the groups
# PCA results --> pca_indiv_2plot
pca_indiv_2plot

# I can just add to this table the genotype and treatment
id <- df.ACQ_data$id
df.pca2anova <- cbind (id, genotype_tt, pca_indiv_2plot)
head (df.pca2anova)
pca2anova.aov <- aov(PC1 ~ genotype_tt + Error(id), data = df.pca2anova)
summary(pca2anova.aov)

pairwise.t.test(df.pca2anova$PC1, df.pca2anova$genotype_tt, , p.adj="hochberg", paired=F)

