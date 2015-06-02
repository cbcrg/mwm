library("ggplot2")

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
                          scale_color_manual(values=c("red", "green", "blue", "lightblue", 
                                                      "magenta", "orange", "gray", "black")) +
                          geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show_guide = T)+
                          theme(legend.key=element_rect(fill=NA)) +
                          labs(title = "PCA of group medians\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
#                           guides(colour = guide_legend(override.aes = list(size = 10)))+
                          guides(colour = guide_legend(override.aes = list(size = 1)))+
                          theme(legend.key=element_rect(fill=NA))
 
#PLOT_paper
pca_medians_acq 
ggsave (pca_medians_acq, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_medians_legend.jpg", sep=""), dpi=900)
ggsave (pca_medians_acq, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_medians_NO_legend.jpg", sep=""), width = 10, height = 10, dpi=900)







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


png("figures/PCAmed_variables.png",res=300,width=15,height=15,unit="cm")
plot(res$var$coord[,1],res$var$coord[,2],asp=1,ylim=c(-1.1,1.1),type="n",main="Mapped variables",xlab="PC1 (80% of variance)",ylab="PC2 (12% of variance)")
arrows(0,0,-res$var$coord[,1],-res$var$coord[,2],length=0.1,angle=10)
textxy(-res$var$coord[,1],-res$var$coord[,2],colnames(M.med),cex=0.9,offset=0.5,col="red")
circle(1,c(0,0))
origin(lty=2)
dev.off()

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
                       labs (title = "Variable contributions\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
                #        geom_polygon(aes(x, y), data = df, inherit.aes = F, Fill=NA)
#                         scale_x_continuous(breaks=1:10)  
                       geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1)
p_circle_plot
ggsave (circle_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "circle_plot.jpg", sep=""), width = 10, height = 10, dpi=900)


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
ggsave (bars_plot, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "bar_contribution.jpg", sep=""), dpi=900, height=5, width=10)

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
ggsave (bars_plot_PC2, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "bar_contribution_PC2.jpg", sep=""), dpi=900, height=5, width=10)

names(-res.pca$ind$coord[,1])
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
#     text(-res.pca$ind.sup$coord[myset,1],-res.pca$ind.sup$coord[myset,2],i,col=cols[j],cex=0.5)
    
  }
}

new_coord$genotype <- factor(new_coord$genotype , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
                            labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))

pca_plot_individuals <- ggplot (data=new_coord, aes (V1, V2)) + 
      geom_text (aes(label=day, colour = genotype), size=3, show_guide = FALSE) +
      scale_color_manual(values=c("red", "green", "blue", "lightblue", 
                        "magenta", "orange", "gray", "black")) +
      xlim (c(-6, 6.5)) + ylim (c(-8, 6.5)) +
      geom_path (data=pca2plot, aes(x=PC1, y=PC2, colour=gentreat),size = 0.5,show_guide = FALSE) +
      labs(title = "Individual variation as supplementary points\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") 
pca_plot_individuals
#PLOT_paper
ggsave (pca_plot_individuals, file=paste(home, "/20150515_PCA_old_frotiersPaper/figures/", "PCA_individuals.jpg", sep=""), dpi=900)



      
      




#permutation test:
di=2
for (a in 1:5){
  set1=which(iddaytreat[,2]=="TS" & iddaytreat[,3]==acq[a])
  set2=which(iddaytreat[,2]=="TSEGCG" & iddaytreat[,3]==acq[a])
  M1=mean(res.pca$ind.sup$coord[set1,di])
  M2=mean(res.pca$ind.sup$coord[set2,di])
  m=length(set1)
  n=length(set2)
  S=sqrt((sum((res.pca$ind.sup$coord[set1,di]-M1)^2)+sum((res.pca$ind.sup$coord[set1,di]-M2)^2))/(m+n-2))
  truestat=((M1-M2)*sqrt(m*n))/(S*sqrt(m+n))
  
  set=c(set1,set2)
  set.seed(317)
  st=c(1:10000)
  st[1]=truestat
  for (i in 2:10000){
    #print(i)
    sam=sample(set)
    set1=sam[1:m]
    set2=sam[(m+1):(m+n)]
    M1=mean(res.pca$ind.sup$coord[set1,di])
    M2=mean(res.pca$ind.sup$coord[set2,di])
    S=sqrt((sum((res.pca$ind.sup$coord[set1,di]-M1)^2)+sum((res.pca$ind.sup$coord[set1,di]-M2)^2))/(m+n-2))
    stat=((M1-M2)*sqrt(m*n))/(S*sqrt(m+n))
    st[i]=stat
  }
  p=max(which (sort(st,decreasing=TRUE)==st[1])/10000)
  print (acq[a])
  print (p)
}

#R objects for Jose:

set=which(iddaytreat[,2]=="TS" & iddaytreat[,3]==acq[a])
m=length(set)
s1=matrix(0,m,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="TS" & iddaytreat[,3]==acq[a])
  s1[,a]=-res.pca$ind.sup$coord[set,1]
}

set=which(iddaytreat[,2]=="TSEEEGCG" & iddaytreat[,3]==acq[a])
n=length(set)
s2=matrix(0,n,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="TSEEEGCG" & iddaytreat[,3]==acq[a])
  s2[,a]=-res.pca$ind.sup$coord[set,1]
}

set=which(iddaytreat[,2]=="TSEE" & iddaytreat[,3]==acq[a])
n=length(set)
s3=matrix(0,n,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="TSEE" & iddaytreat[,3]==acq[a])
  s3[,a]=-res.pca$ind.sup$coord[set,1]
}

set=which(iddaytreat[,2]=="TSEGCG" & iddaytreat[,3]==acq[a])
n=length(set)
s4=matrix(0,n,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="TSEGCG" & iddaytreat[,3]==acq[a])
  s4[,a]=-res.pca$ind.sup$coord[set,1]
}

set=which(iddaytreat[,2]=="WT" & iddaytreat[,3]==acq[a])
n=length(set)
s5=matrix(0,n,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="WT" & iddaytreat[,3]==acq[a])
  s5[,a]=-res.pca$ind.sup$coord[set,1]
}

TS=s1
TSEEEGCG=s2
TSEE=s3
TSEGCG=s4
WT=s5
save(TS,TSEEEGCG,TSEE,TSEGCG,WT,file="5setsPC1.R")

set=which(iddaytreat[,2]=="WTEE" & iddaytreat[,3]==acq[a])
n=length(set)
s6=matrix(0,n,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="WTEE" & iddaytreat[,3]==acq[a])
  s6[,a]=-res.pca$ind.sup$coord[set,1]
}

WTEE=s6

set=which(iddaytreat[,2]=="WTEGCG" & iddaytreat[,3]==acq[a])
n=length(set)
s7=matrix(0,n,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="WTEGCG" & iddaytreat[,3]==acq[a])
  s7[,a]=-res.pca$ind.sup$coord[set,1]
}

WTEGCG=s7

set=which(iddaytreat[,2]=="WTEEEGCG" & iddaytreat[,3]==acq[a])
n=length(set)
s8=matrix(0,n,5)
for (a in 1:5){
  set=which(iddaytreat[,2]=="WTEEEGCG" & iddaytreat[,3]==acq[a])
  s8[,a]=-res.pca$ind.sup$coord[set,1]
}

WTEEEGCG=s8

save(WTEEEGCG,WTEE,WTEGCG,file="3setsPC1.R")




boxplot(s1,ylim=c(-6.2,7.2),ylab="PC1",xlab="Acquisition day")
boxplot(s2,add=TRUE,col=rgb(0,1,0,0.5))
legend("topleft",pch=0,col=c("green","white"),c("TS","TSEEEGCG"))





#analyse removal session:

ma3=spss.get("Jtracks parameters except latency.sav")
MR=ma3[1:83,c("DIST.REM","GALLINDEX.REM","SPEED.REM","PERC.NE.REM","PERC.PERI.REM","WISHAW.REM")]
gt=sort(unique(idgentreat[,2]))
idgentreat=ma3[1:83,c("Mouse.ID","GEN.TREAT")]

MR.med=matrix(0,8,6)
for (i in 1:length(gt)){
  print (i)
  MR.med[i,]=apply(MR[idgentreat[1:83,2]==gt[i],],2,median)
  
}
colnames(MR.med)=colnames(MR)
rownames(MR.med)=gt

joint=rbind(MR.med,MR)
resR = PCA(joint, scale.unit=TRUE, ind.sup=c(9:91))

#permutation test:
di=3
set1=which(idgentreat[,2]=="TS")
set2=which(idgentreat[,2]=="TSEEEGCG")
M1=mean(resR$ind.sup$coord[set1,di])
M2=mean(resR$ind.sup$coord[set2,di])
m=length(set1)
n=length(set2)
S=sqrt((sum((resR$ind.sup$coord[set1,di]-M1)^2)+sum((resR$ind.sup$coord[set1,di]-M2)^2))/(m+n-2))
truestat=((M1-M2)*sqrt(m*n))/(S*sqrt(m+n))

set=c(set1,set2)
set.seed(317)
st=c(1:10000)
st[1]=truestat
for (i in 2:10000){
  #print(i)
  sam=sample(set)
  set1=sam[1:m]
  set2=sam[(m+1):(m+n)]
  M1=mean(resR$ind.sup$coord[set1,di])
  M2=mean(resR$ind.sup$coord[set2,di])
  S=sqrt((sum((resR$ind.sup$coord[set1,di]-M1)^2)+sum((resR$ind.sup$coord[set1,di]-M2)^2))/(m+n-2))
  stat=((M1-M2)*sqrt(m*n))/(S*sqrt(m+n))
  st[i]=stat
}
p=max(which (sort(st)==st[1])/10000)
print (p)


#check original variables: Gallagher index: 

a=4
set1=which(iddaytreat[,2]=="TS" & iddaytreat[,3]==acq[a])
set2=which(iddaytreat[,2]=="TSEEEGCG" & iddaytreat[,3]==acq[a])
GA1=as.numeric(M.ind[set1,"gallindex"])
GA2=as.numeric(M.ind[set2,"gallindex"])

a=4
set1=which(iddaytreat[,2]=="TS" & iddaytreat[,3]==acq[a])
set2=which(iddaytreat[,2]=="TSEEEGCG" & iddaytreat[,3]==acq[a])
LA1=as.numeric(M.ind[set1,"latency"])
LA2=as.numeric(M.ind[set2,"latency"])


#discriminant analysis trisomic:

M.ind=as.matrix(ma2[ma2$day%in%acq,var])
iddaytreat=as.matrix(ma2[ma2$day%in%acq,1:3])
rownames(iddaytreat)=c(1:415)

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



