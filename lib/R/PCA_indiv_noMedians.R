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
df.ACQ_data <- df.ACQ_data []

data_all_sessions_var <- data_all_sessions [ , c(9:length(data_all_sessions[1,]))]




data_all_sessions [ , c(9:length(data_all_sessions[1,]))]
cbind (data_all_sessions_var$GEN.TREAT, data_all_sessions_var$)
merge (df.ACQ_data_except_latency, acq_latency, by=c("id"))







acq_latency <- merge(lat_1, lat_2, by=c("id"))



acq_latency <- merge(acq_latency, lat_4, by=c("id"))
acq_latency <- merge(acq_latency, lat_5, by=c("id"))



                     lat_3, lat_4, lat_5, by=c("id"))

""


day1_latency <- subset(dat, Subject==2, latency)


ma2 []
studentdata[studentdata$Drink == 'water',]



session <- colnames (ma2)
gentreat <- ma2$gentreat


# transpose all but the first column (name)
t.ma2 <- as.data.frame(t(ma2))
head(t.ma2)
row.names(t.ma2)

# Get only the rowname of the variable that we want and add to the first table as column
# lantency is row 7
latency <- t.ma2[7,]
id <- t.ma2[1,]

day1[studentdata$Drink == 'water',]
t(rbind (id,latency))

rbind (row.names (t.ma2)[4], t.ma2$day)



colnames(df.aree) <- n
df.aree$myfactor <- factor(row.names(df.aree))

str(df.aree)

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

