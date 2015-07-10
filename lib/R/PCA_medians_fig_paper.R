#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. May 2015                         ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### PCA analysis figures for frontiers paper                              ### 
#############################################################################

library("ggplot2")
library("Hmisc")

##Getting HOME directory
home <- Sys.getenv("HOME")

#path of R session is /users/cn/ierb/work/MaraDierssen/Silvina/

# Loading functions:
source (paste (home, "/git/mwm/lib/R/plot_param_public.R", sep=""))

library(Hmisc)
#install.packages("calibrate")
library(calibrate)
ma2=spss.get("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/Ts65Dn OLD ACQ1_ACQ5_SUBCONJ.sav")

M=as.matrix(ma2[ma2$day=="Day 5",c(4,5,7,8,9,10,11,12)])
treat=as.matrix(ma2[ma2$day=="Day 5",2])
M4=as.matrix(ma2[ma2$day=="Day 4",c(4,5,7,8,9,10,11,12)])
treat4=as.matrix(ma2[ma2$day=="Day 4",2])
M3=as.matrix(ma2[ma2$day=="Day 3",c(4,5,7,8,9,10,11,12)])
treat3=as.matrix(ma2[ma2$day=="Day 3",2])
M2=as.matrix(ma2[ma2$day=="Day 2",c(4,5,7,8,9,10,11,12)])
treat2=as.matrix(ma2[ma2$day=="Day 2",2])
M1=as.matrix(ma2[ma2$day=="Day 1",c(4,5,7,8,9,10,11,12)])
treat1=as.matrix(ma2[ma2$day=="Day 1",2])

MM5=matrix(0,8,8)
MM5[1,]=apply(M[treat=="WT",],2,median)
MM5[2,]=apply(M[treat=="WTEE",],2,median)
MM5[3,]=apply(M[treat=="WTEGCG",],2,median)
MM5[4,]=apply(M[treat=="WTEEEGCG",],2,median)
MM5[5,]=apply(M[treat=="TS",],2,median)
MM5[6,]=apply(M[treat=="TSEE",],2,median)
MM5[7,]=apply(M[treat=="TSEGCG",],2,median)
MM5[8,]=apply(M[treat=="TSEEEGCG",],2,median)
MM1=matrix(0,8,8)
MM1[1,]=apply(M1[treat1=="WT",],2,median)
MM1[2,]=apply(M1[treat1=="WTEE",],2,median)
MM1[3,]=apply(M1[treat1=="WTEGCG",],2,median)
MM1[4,]=apply(M1[treat1=="WTEEEGCG",],2,median)
MM1[5,]=apply(M1[treat1=="TS",],2,median)
MM1[6,]=apply(M1[treat1=="TSEE",],2,median)
MM1[7,]=apply(M1[treat1=="TSEGCG",],2,median)
MM1[8,]=apply(M1[treat1=="TSEEEGCG",],2,median)
MM2=matrix(0,8,8)
MM2[1,]=apply(M2[treat2=="WT",],2,median)
MM2[2,]=apply(M2[treat2=="WTEE",],2,median)
MM2[3,]=apply(M2[treat2=="WTEGCG",],2,median)
MM2[4,]=apply(M2[treat2=="WTEEEGCG",],2,median)
MM2[5,]=apply(M2[treat2=="TS",],2,median)
MM2[6,]=apply(M2[treat2=="TSEE",],2,median)
MM2[7,]=apply(M2[treat2=="TSEGCG",],2,median)
MM2[8,]=apply(M2[treat2=="TSEEEGCG",],2,median)
MM3=matrix(0,8,8)
MM3[1,]=apply(M3[treat3=="WT",],2,median)
MM3[2,]=apply(M3[treat3=="WTEE",],2,median)
MM3[3,]=apply(M3[treat3=="WTEGCG",],2,median)
MM3[4,]=apply(M3[treat3=="WTEEEGCG",],2,median)
MM3[5,]=apply(M3[treat3=="TS",],2,median)
MM3[6,]=apply(M3[treat3=="TSEE",],2,median)
MM3[7,]=apply(M3[treat3=="TSEGCG",],2,median)
MM3[8,]=apply(M3[treat3=="TSEEEGCG",],2,median)
MM4=matrix(0,8,8)
MM4[1,]=apply(M4[treat4=="WT",],2,median)
MM4[2,]=apply(M4[treat4=="WTEE",],2,median)
MM4[3,]=apply(M4[treat4=="WTEGCG",],2,median)
MM4[4,]=apply(M4[treat4=="WTEEEGCG",],2,median)
MM4[5,]=apply(M4[treat4=="TS",],2,median)
MM4[6,]=apply(M4[treat4=="TSEE",],2,median)
MM4[7,]=apply(M4[treat4=="TSEGCG",],2,median)
MM4[8,]=apply(M4[treat4=="TSEEEGCG",],2,median)
MMALL=rbind(MM1,MM2,MM3,MM4,MM5)
colnames(MMALL)=c("dist","gallindex","latency","speed","percentne","percenter","wishaw","percentperi")
rownames(MMALL)=c("WT1","WTEE1","WTEGCG1","WTEEEGCG1","TS1","TSEE1","TSEGCG1","TSEEEGCG1","WT2","WTEE2","WTEGCG2","WTEEEGCG2","TS2","TSEE2","TSEGCG2","TSEEEGCG2","WT3","WTEE3","WTEGCG3","WTEEEGCG3","TS3","TSEE3","TSEGCG3","TSEEEGCG3","WT4","WTEE4","WTEGCG4","WTEEEGCG4","TS4","TSEE4","TSEGCG4","TSEEEGCG4","WT5","WTEE5","WTEGCG5","WTEEEGCG5","TS5","TSEE5","TSEGCG5","TSEEEGCG5")
ResMMALL=prcomp(MMALL,scale=TRUE)

pca2plot <- as.data.frame(ResMMALL$x)
row.names(pca2plot)
pca2plot$gen_day <- row.names(pca2plot)

pca2plot$days <-  as.factor(as.numeric (gsub(".*([0-9]+)$", "\\1", pca2plot$gen_day)))
pca2plot$gentreat <-  as.factor(gsub("([A-Z]+).*$", "\\1", pca2plot$gen_day))
pca2plot$PC1

pca2plot$gentreat <- factor(pca2plot$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                            labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_medians_acq <- ggplot(pca2plot, aes(x=PC1, y=PC2, colour=gentreat )) + 
                          geom_path (size = 1,show_guide = T) + 
#                           geom_path (size = 1,show_guide = F) + 
                          scale_color_manual(values=c("red", "darkgreen", "blue", "lightblue", 
                                                      "magenta", "orange", "gray", "black")) +
                          geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show_guide = T)+
#                           geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show_guide = F)+
                          theme(legend.key=element_rect(fill=NA)) +
                          labs(title = "PCA of group medians\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
#                           guides(colour = guide_legend(override.aes = list(size = 10)))+
                          guides(colour = guide_legend(override.aes = list(size = 1)))+
                          theme(legend.key=element_rect(fill=NA))

#PLOT_paper
pca_medians_acq 
ggsave (pca_medians_acq, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig1_PCA/", "PCA_medians_legend.jpg", sep=""), dpi=900)
# ggsave (pca_medians_acq, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig1_PCA/", "PCA_medians_NO_legend.jpg", sep=""), width = 10, height = 10, dpi=900)

# Plot for presentation
# setwd("/Users/jespinosa/Dropbox (Personal)/presentations_2015/20150630_GM_Cedric/figures")
pca_medians_acq <- ggplot(pca2plot, aes(x=PC1, y=PC2, colour=gentreat )) + 
  geom_path (size = 2,show_guide = T) + 
  scale_color_manual(values=c("red", "green", "blue", "lightblue", 
                              "magenta", "orange", "gray", "black")) +
  geom_text (aes (label=days), vjust=-0.5, hjust=1, size=5, show_guide = T)+
  theme(legend.key=element_rect(fill=NA)) +
  labs(title = "PCA of group medians\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
  #                           guides(colour = guide_legend(override.aes = list(size = 10)))+
  guides(colour = guide_legend(override.aes = list(size = 2)))+
  theme(legend.key=element_rect(fill=NA))

#PLOT_presentation
pca_medians_acq 
# ggsave (pca_medians_acq, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_medians_legend.jpg", sep=""), dpi=900)
# ggsave (pca_medians_acq, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_medians_NO_legend.jpg", sep=""), width = 10, height = 10, dpi=900)







############################

mWT=rbind(M1[treat1=="WT",],M2[treat2=="WT",],M3[treat3=="WT",],M4[treat4=="WT",],M[treat=="WT",])

sMMALL=matrix(0,dim(MMALL)[1],dim(MMALL)[2])
for (i in 1:dim(MMALL)[1]){
  sMMALL[i,]=(MMALL[i,]-apply(MMALL,2,mean))/apply(MMALL,2,sd)
}

#this is the same scaling that is employed when doing prcomp

#Now: we have to use these means and variances for the original data
smWT=matrix(0,dim(mWT)[1],dim(mWT)[2])
for (i in 1:dim(mWT)[1]){
  smWT[i,]=(mWT[i,]-apply(MMALL,2,mean))/apply(MMALL,2,sd)
}
#check:
apply(smWT[1:10,],2,median)%*%ResMMALL$rotation[,1]==ResMMALL$x[1,1]
apply(smWT[11:20,],2,median)%*%ResMMALL$rotation[,1]==ResMMALL$x[9,1]

############################################################################################################

#plot all data as supplementary points

acq=c("Day 1","Day 2","Day 3","Day 4","Day 5")
var=c(4,5,7,8,9,11,12)

M.ind=as.matrix(ma2[ma2$day%in%acq,var])
iddaytreat=as.matrix(ma2[ma2$day%in%acq,1:3])
rownames(iddaytreat)=c(1:415)
tgt=sort(unique(iddaytreat[,2]))



M.med=matrix(0,40,7)
rnames=c(1:40)
for (i in 1:length(acq)){
  for (j in 1:length(tgt)){
    ind=(i-1)*length(tgt)+j
    print (ind)
    M.med[ind,]=apply(M.ind[iddaytreat[,3]==acq[i] & iddaytreat[,2]==tgt[j],],2,median)
    rnames[ind]=paste(tgt[j],acq[i]," ")
  }
}
colnames(M.med)=colnames(M.ind)
rownames(M.med)=rnames

jm=rbind(M.med,M.ind)


library(FactoMineR)
res = PCA(jm, scale.unit=TRUE, ind.sup=c(41:455)) 
plot(res,choix="var")
cols=c("green","lightblue","black","orange","red","blue","darkgrey","magenta")


# png("figures/PCAmed_variables.png",res=300,width=15,height=15,unit="cm")
plot(res$var$coord[,1],res$var$coord[,2],asp=1,ylim=c(-1.1,1.1),type="n",main="Mapped variables",xlab="PC1 (80% of variance)",ylab="PC2 (12% of variance)")
arrows(0,0,-res$var$coord[,1],-res$var$coord[,2],length=0.1,angle=10)
textxy(-res$var$coord[,1],-res$var$coord[,2],colnames(M.med),cex=0.9,offset=0.5,col="red")
circle(1,c(0,0))
origin(lty=2)
# dev.off()

row.names(res$var$coord)
circle_plot <- as.data.frame (res$var$coord)
labels_v <- row.names(res$var$coord)
labels_v[6] <- "whishaw"
circle_plot$labels <- lab_names

neg_labels <- labels_v [c(1,2,3,7)]
neg_positions <- circle_plot [c(1,2,3,7), c(1,2)]
# change positions for labels
neg_positions [2,2] <-neg_positions [2,2] + 0.04 
neg_positions [3,2] <-neg_positions [3,2] - 0.04
pos_labels <- labels_v [c(4,5,6)]
pos_positions <- circle_plot [c(4,5,6), c(1,2)]

angle <- seq(-pi, pi, length = 50)
df.circle <- data.frame(x = sin(angle), y = cos(angle))

#aes(x=PC1, y=PC2, colour=gentreat )) 
p_circle_plot <- ggplot(circle_plot) + 
                       geom_segment (data=circle_plot, aes(x=0, y=0, xend=-Dim.1, yend=-Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, color="red") +
                       xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
                       geom_text (data=neg_positions, aes (x=-Dim.1, y=-Dim.2, label=neg_labels, hjust=1.2), show_guide = FALSE, size=5) + 
                       geom_text (data=pos_positions, aes (x=-Dim.1, y=-Dim.2, label=pos_labels, hjust=-0.3), show_guide = FALSE, size=5) +
                       geom_vline (xintercept = 0, linetype="dotted") +
                       geom_hline (yintercept=0, linetype="dotted") +
                       labs (title = "PCA of the variables\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
                #        geom_polygon(aes(x, y), data = df, inherit.aes = F, Fill=NA)
#                         scale_x_continuous(breaks=1:10)  
                       geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)
p_circle_plot
# ggsave (p_circle_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)


############
## BARPLOT
df.bars <- cbind (as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE)), names(res$var$coord[,1])[order(res$var$coord[,1]^2,decreasing=TRUE)])
df.bars[6,2] <- "whishaw"
df.bars_to_plot <- as.data.frame(df.bars)
df.bars_to_plot$index <- as.factor (df.bars_to_plot$V2)
class (df.bars_to_plot$V1)
df.bars_to_plot$value <- as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE))

df.bars_to_plot$index <- factor(df.bars_to_plot$index , levels=c("gallindex", "latency", "percentne", "dist", "percentperi", "whishaw", "speed"), 
                            labels=c("gallindex", "latency", "percentne", "dist", "percentperi", "whishaw", "speed"))


bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
                    geom_bar (stat="identity", fill="gray", width=0.8) + 
                    labs (title = "Variable contribution to PC1\n", x = "", y="Contribution in %\n") +
                     theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot

#PLOT_paper
# ggsave (bars_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "bar_contribution.jpg", sep=""), dpi=900, height=5, width=10)

df.bars_PC2 <- cbind (as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE)), names(res$var$coord[,2])[order(res$var$coord[,2]^2,decreasing=TRUE)])
df.bars_PC2[2,2] <- "whishaw"
df.bars_to_plot_PC2 <- as.data.frame(df.bars_PC2)
df.bars_to_plot_PC2$index <- as.factor (df.bars_to_plot_PC2$V2)
# class (df.bars_to_plot_PC2$V1)
# df.bars_to_plot_PC2$value <- as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE))
df.bars_to_plot_PC2$value <- as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE))

df.bars_to_plot_PC2$index
df.bars_to_plot_PC2$index <- factor(df.bars_to_plot$index , levels=c("speed", "whishaw", "percentperi", "dist", "percentne", "latency", "gallindex"), 
                                labels=c("speed", "whishaw", "percentperi", "dist", "percentne", "latency", "gallindex"))


df.bars_to_plot_PC2$value <- rev(df.bars_to_plot_PC2$value)
bars_plot_PC2 <- ggplot (data=df.bars_to_plot_PC2, aes(x=index, y=value)) + 
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = "Variable contribution to PC2\n", x = "", y="Contribution in %\n") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1) )
bars_plot_PC2

#PLOT_paper
# ggsave (bars_plot_PC2, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "bar_contribution_PC2.jpg", sep=""), dpi=900, height=5, width=10)

names(-res.pca$ind$coord[,1])

#############################
# Plot of all gays
# png("figures/PCAmed_individuals.png",res=300,width=15,height=15,unit="cm")
# dataframe creation

new_coord <- cbind(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2])
new_coord

res

plot(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2],type="n",main="Individual variation as supplementary points",xlab="PC1 (80% of variance)",ylab="PC2 (12% of variance)")

day <- c()
genotype <- c()
new_coord <- as.data.frame( cbind(-res$ind.sup$coord[,1],-res$ind.sup$coord[,2]))
new_coord$day <- c(1:415)
new_coord$genotype <- c()

# myset=which(iddaytreat[,2]==tgt[1] & iddaytreat[,3]==acq[1])
myset<-c()
for (j in 1:length(tgt)){
  for (i in 1:length(acq)){
    if (i<5){
      ind=(i-1)*length(tgt)+j
      ind2=i*length(tgt)+j 
      #       lines(-res.pca$ind$coord[c(ind,ind2),1],-res.pca$ind$coord[c(ind,ind2),2],col=cols[j])
#       day <- c(day, i)
#       genotype <- c(genotype, tgt[j])
    }
    myset=which(iddaytreat[,2]==tgt[j] & iddaytreat[,3]==acq[i])
    new_coord [myset,c("day")] <- i
    new_coord [myset,c("genotype")] <- tgt[j]
    new_coord [myset,c("id")] <- substr(iddaytreat[myset,1], 7, 9)
#     text(-res.pca$ind.sup$coord[myset,1],-res.pca$ind.sup$coord[myset,2],i,col=cols[j],cex=0.5)
    
  }
}

# substr(new_coord$id, 7,9)

new_coord$genotype <- factor(new_coord$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                            labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_plot_individuals <- ggplot (data=new_coord, aes (V1, V2)) + 
      geom_text (aes(label=day, colour = genotype), size=5, show_guide = FALSE) +
      scale_color_manual(values=c("red", "darkgreen", "blue", "lightblue", 
                        "magenta", "orange", "gray", "black")) +
      xlim (c(-6, 6.5)) + ylim (c(-8, 6.5)) +
      geom_path (data=pca2plot, aes(x=PC1, y=PC2, colour=gentreat),size = 1,show_guide = FALSE) +
      labs(title = "Individual as supplementary points\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
      coord_fixed()
pca_plot_individuals

#PLOT_paper
ggsave (pca_plot_individuals, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig_1_PCA_sup/", "PCA_individuals.jpg", sep=""), height = 10, width = 10, dpi=900)

####
# Adding cloud to the individuals plot
# By group, instead of the line let's see how it looks like
pca_plot_individuals <- ggplot (data=new_coord, aes (V1, V2)) + 
  geom_text (aes(label=day, colour = genotype), size=3, show_guide = FALSE) +
  scale_color_manual(values=c("red", "green", "blue", "lightblue", 
                              "magenta", "orange", "gray", "black")) +
  xlim (c(-6, 6.5)) + ylim (c(-8, 6.5)) +
  geom_path (data=pca2plot, aes(x=PC1, y=PC2, colour=gentreat),size = 0.5,show_guide = FALSE) +
  labs(title = "Individual variation as supplementary points\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") 
pca_plot_individuals
head(new_coord)

pca2plot$genotype <- pca2plot$gentreat

p_cloud_indiv_by_day <- ggplot(new_coord, aes(V1, V2, color=genotype, label=day)) + 
  stat_density2d(aes(fill=factor(genotype), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype', 
                     values = c("red", "green", "blue", "lightblue", 
                                "magenta", "orange", "yellow", "black"), 
                     labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("red", "green", "blue", "lightblue", 
                                "magenta", "orange", "yellow", "black"), 
                     labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
  geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "Density plot of individual variation\n", x = "\nPC1", y="PC2\n") +
  scale_alpha_continuous(range=c(0.3,0.5)) +
  geom_path (data=pca2plot, aes(x=PC1, y=PC2), colour="black", size = 0.5, linetype = 2, show_guide = FALSE)
#  geom_path (data=pca2plot, aes(x=PC1, y=PC2, colour=gentreat), size = 0.6, show_guide = FALSE) +
#             scale_color_manual (name='genotype', values =c (rep ("black",7), "white")) 

p_cloud_indiv_by_day_facet <- p_cloud_indiv_by_day + facet_wrap(~genotype, ncol = 2)
p_cloud_indiv_by_day_facet

#PLOT_presentation
# ggsave (p_cloud_indiv_by_day_facet, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_cloud_indiv_byDay.jpg", sep=""), width = 10, height = 10, dpi=900)


##############
# Adding the plot of PCA cloud for day A1 and day A5
new_coord

PC1_acq1 <- subset(new_coord, day=="1", c("V1", "V2", "day","genotype", "id"))
colnames (PC1_acq1) <- c("PC1","PC2", "day", "genotype_tt", "id")
p_cloud_acq1 <- ggplot(PC1_acq1, aes(PC1, PC2, color=genotype_tt, label=id)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("red", "green", "blue", "lightblue", 
                                "magenta", "orange", "yellow", "black"), 
                     labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("red", "green", "blue", "lightblue", 
                                "magenta", "orange", "yellow", "black"), 
                     labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
  geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "PCA coordinates density, day 1\n", x = "\nPC1", y="PC2\n")

p_cloud_acq1

PC1_acq5 <- subset(new_coord, day=="5", c("V1", "V2", "day","genotype", "id"))
colnames (PC1_acq5) <- c("PC1","PC2", "day", "genotype_tt", "id")
p_cloud_acq5 <- ggplot(PC1_acq5, aes(PC1, PC2, color=genotype_tt, label=id)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) + 
#   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("red", "green", "blue", "lightblue", 
                                "magenta", "orange", "yellow", "black"), 
                     labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("red", "green", "blue", "lightblue", 
                                "magenta", "orange", "yellow", "black"), 
                     labels = c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG")) + 
  geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "PCA coordinates density, trisomic group, day 5\n", x = "\nPC1", y="PC2\n")

p_cloud_acq1
p_cloud_acq5

#PLOT_presentation
# ggsave (p_cloud_acq1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_acq1_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)
# ggsave (p_cloud_acq5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_acq5_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)

##############
# Adding the plot of PCA cloud for day A1 and day A5 only for trisomics
new_coord
indiv <- rownames(pca_indiv_2plot)
new_coord_indiv <- rbind(indiv, new_coord)
PC1_TS  <- subset(new_coord, grepl("TS", genotype))
PC1_TS_acq1 <- subset(PC1_TS, day=="1", c("V1", "V2", "day","genotype", "id"))

colnames (PC1_TS_acq1) <- c("PC1","PC2", "day", "genotype_tt", "id")

# Link to the explanation of how to perform density plots with ggplot2
# http://stackoverflow.com/questions/19791181/density-shadow-around-the-data-with-ggplot2-r
# Quick facts
# n controls the smoothness of the density polygon.
# h is the bandwidth of the density estimation.
# bins controls the number of density levels.
# alpha transparency allows to make the cloud more or less transparent depending on the level
# level, Computed density, is the amount you have to increase to change from one level to the next

# p_cloud_ts_acq1 <- ggplot(PC1_TS_acq1, aes(PC1, PC2, color=genotype_tt, label=id)) + 
p_cloud_ts_acq1 <- ggplot(PC1_TS_acq1, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
#                  geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) +
                 geom="polygon", color=NA, h=5, n=100, bins=6, show_guide = FALSE) +
#   scale_alpha(range = c(0.00, 1)) +
#   scale_size(range = c(0, 0.01), guide = "none")  
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values=c("darkgreen", "lightblue", "darkorange", "black"), 
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values=c("darkgreen", "lightblue", "orange", "black"),
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
#   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "PCA coordinates density, trisomic group, day 1\n", x = "\nPC1", y="PC2\n")

p_cloud_ts_acq1 <- p_cloud_ts_acq1 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_ts_acq1

PC1_TS_acq5 <- subset(PC1_TS, day=="5", c("V1", "V2", "day","genotype", "id"))

colnames (PC1_TS_acq5) <- c("PC1","PC2", "day", "genotype_tt", "id")
# p_cloud_ts_acq5 <- ggplot(PC1_TS_acq5, aes(PC1, PC2, color=genotype_tt, label=id)) + 
p_cloud_ts_acq5 <- ggplot(PC1_TS_acq5, aes(PC1, PC2, color=genotype_tt)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=5, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("darkgreen", "lightblue", "darkorange", "black"), 
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("darkgreen", "lightblue","orange", "black"),
                     labels = c("TS", "TSEE", "TSEGCG", "TSEEEGCG")) + 
#   geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-6, 8.5),
                  ylim=c(-9, 6)) +
  labs(title = "PCA coordinates density, trisomic group, day 5\n", x = "\nPC1", y="PC2\n")

p_cloud_ts_acq1 + scale_alpha_continuous(range=c(0.3,0.5))
# p_cloud_ts_acq5 + scale_alpha_continuous(range=c(0.3,0.5))
p_cloud_ts_acq5 <- p_cloud_ts_acq5 + facet_wrap(~genotype_tt, ncol = 2)  + geom_vline(xintercept = 0, colour="gray") + geom_hline(yintercept = 0, colour="gray")
p_cloud_ts_acq5

# Plot for paper
setwd("/Users/jespinosa/20150515_PCA_old_frotiersPaper/figures/fig2_PCA")

#PLOT_paper
ggsave (p_cloud_ts_acq1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "PCA_acq1_ts_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)
ggsave (p_cloud_ts_acq5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/fig2_PCA/", "PCA_acq5_ts_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)

## Only comparison between TS and TSEEEGCG day 1 vs day 5
# Day 1
PC1_TS_acq1 <- subset(PC1_TS, day=="1", c("V1", "V2", "day","genotype", "id"))

PC1_TS_TSEEEGCG_acq1  <- subset(PC1_TS_acq1, genotype=="TS" | genotype=="TSEEEGCG")

colnames (PC1_TS_TSEEEGCG_acq1) <- c("PC1","PC2", "day", "genotype_tt", "id")

p_cloud_ts_tseeegcg_acq1 <- ggplot(PC1_TS_TSEEEGCG_acq1, aes(PC1, PC2, color=genotype_tt, label=id)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("green", "black"), 
                     labels = c("TS", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("green", "black"),
                     labels = c("TS",  "TSEEEGCG")) + 
  geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "PCA coordinates density, trisomic not treated\nversus trisomic EE and EGCG, day 1\n", x = "\nPC1", y="PC2\n")

p_cloud_ts_tseeegcg_acq1

#PLOT_presentation
# ggsave (p_cloud_ts_tseeegcg_acq1, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_acq1_ts_tseeegcg_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)

# Day 5
PC1_TS_acq5 <- subset(PC1_TS, day=="5", c("V1", "V2", "day", "genotype", "id"))

PC1_TS_TSEEEGCG_acq5  <- subset(PC1_TS_acq5, genotype=="TS" | genotype=="TSEEEGCG")

colnames (PC1_TS_TSEEEGCG_acq5) <- c("PC1","PC2", "day", "genotype_tt", "id")
p_cloud_ts_tseeegcg_acq5 <- ggplot(PC1_TS_TSEEEGCG_acq5, aes(PC1, PC2, color=genotype_tt, label=id)) + 
  stat_density2d(aes(fill=factor(genotype_tt), alpha = ..level..), 
                 geom="polygon", color=NA, n=100, h=4, bins=6, show_guide = FALSE) + 
  #   geom_smooth(se=F, method='lm', show_guide = FALSE) + 
  geom_point(show_guide = FALSE) + 
  scale_color_manual(name='genotype_tt', 
                     values = c("green", "black"), 
                     labels = c("TS", "TSEEEGCG")) + 
  scale_fill_manual( name='gentreat', 
                     values = c("green", "black"),
                     labels = c("TS",  "TSEEEGCG")) + 
  geom_text(hjust=0.5, vjust=-1 ,size=3, color="black") + 
  scale_x_continuous(expand=c(0.3, 0)) + # Zooms out so that density polygons
  scale_y_continuous(expand=c(0.3, 0)) + # don't reach edges of plot.
  coord_cartesian(xlim=c(-7, 9),
                  ylim=c(-10, 10)) +
  labs(title = "PCA coordinates density, trisomic not treated\nversus trisomic EE and EGCG, day 1\n", x = "\nPC1", y="PC2\n")

p_cloud_ts_tseeegcg_acq5

#PLOT_presentation
# ggsave (p_cloud_ts_tseeegcg_acq5, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_acq5_ts_tseeegcg_cloud.jpg", sep=""), width = 10, height = 10, dpi=900)












lat_2 <- subset(ma2, day=="Day 2", c("id", "latency"))




new_coord
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

