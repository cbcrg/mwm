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
                          geom_path (size = 1.5) + 
                          scale_color_manual(values=c("red", "green", "blue", "lightgreen", 
                                                      "pink", "orange", "gray", "black")) +
                          geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4)+
                          theme(legend.key=element_rect(fill=NA)) +
                          labs(title = "PCA of group medians\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
#                           guides(colour = guide_legend(override.aes = list(size = 10)))+
                          guides(colour = guide_legend(override.aes = list(size = 1)))+
                          theme(legend.key=element_rect(fill=NA))
                         
pca_medians_acq 
ggsave (pca_medians_acq, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_medians.jpg", sep=""), dpi=900)
