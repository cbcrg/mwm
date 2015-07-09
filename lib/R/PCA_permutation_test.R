#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. July 2015                        ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes                 ###
###                                                                       ###
#############################################################################

ma2=spss.get("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/Ts65Dn OLD ACQ1_ACQ5_SUBCONJ.sav")

head(ma2)

# permutations
p<-1000

set.seed(111)
# id_group<-data.frame()

id_group <- subset(ma2, grepl("PRE", ma2$day))

id_group <- id_group [,c(1,2)]

p_genetreat <- sample(id_group$gentreat)
id_group$gentreat <- p_genetreat

p_ma2 <- merge(id_group, ma2, by="id", all = TRUE)
p_ma2 <- p_ma2 [,c(-1)]
head (p_ma2)
M=as.matrix(p_ma2[p_ma2$day=="Day 5",c(4,5,7:12)])
treat=as.matrix(p_ma2[p_ma2$day=="Day 5",2])
M4=as.matrix(p_ma2[p_ma2$day=="Day 4",c(4,5,7:12)])
treat4=as.matrix(p_ma2[p_ma2$day=="Day 4",2])
M3=as.matrix(p_ma2[p_ma2$day=="Day 3",c(4,5,7:12)])
treat3=as.matrix(p_ma2[p_ma2$day=="Day 3",2])
M2=as.matrix(p_ma2[p_ma2$day=="Day 2",c(4,5,7:12)])
treat2=as.matrix(p_ma2[p_ma2$day=="Day 2",2])
M1=as.matrix(p_ma2[p_ma2$day=="Day 1",c(4,5,7:12)])
treat1=as.matrix(p_ma2[p_ma2$day=="Day 1",2])

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


# pca2plot <- as.data.frame(ResMMALL$x)
# row.names(pca2plot)
# pca2plot$gen_day <- row.names(pca2plot)

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






for (i in 1:p) {
  id_group <- subset(ma2, grepl("PRE", ma2$day))
  id_group <- id_group [,c(1,2)]  
  p_genetreat <- sample(id_group$gentreat)
  id_group$gentreat <- p_genetreat
  
  p_ma2 <- merge(id_group, ma2, by="id", all = TRUE)
  
  # gentreat.y hold new labels
  M=as.matrix(p_ma2[p_ma2$day=="Day 5",c(5,6,8:13)])
  treat=as.matrix(p_ma2[p_ma2$day=="Day 5",3])
  M4=as.matrix(p_ma2[p_ma2$day=="Day 4",c(5,6,8:13)])
  treat4=as.matrix(p_ma2[p_ma2$day=="Day 4",3])
  M3=as.matrix(p_ma2[p_ma2$day=="Day 3",c(5,6,8:13)])
  treat3=as.matrix(p_ma2[p_ma2$day=="Day 3",3])
  M2=as.matrix(p_ma2[p_ma2$day=="Day 2",c(5,6,8:13)])
  treat2=as.matrix(p_ma2[p_ma2$day=="Day 2",3])
  M1=as.matrix(p_ma2[p_ma2$day=="Day 1",c(5,6,8:13)])
  treat1=as.matrix(p_ma2[p_ma2$day=="Day 1",3])
  
  MM5=matrix(0,8,8)
  ?apply
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
  
  
  
  
  
} 

# t statistic calculation
PCA
set.seed(111)

sample(ma2$gentreat)

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
