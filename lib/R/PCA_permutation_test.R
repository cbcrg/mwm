#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. July 2015                        ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### t-test statistic calculation for two PCA group comparison by doing    ### 
### several permutations of the original animal genotypes                 ###
###                                                                       ###
#############################################################################

library(FactoMineR)

ma2=spss.get("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/Ts65Dn OLD ACQ1_ACQ5_SUBCONJ.sav")

# Function t statistic calculation
f_t_stat <- function (df_coord, gen_1 = "TS", gen_2 = "TSEEEGCG"){
  group1 <- subset (new_coord, genotype == gen_1)
  group2 <- subset (new_coord, genotype == gen_2)
  t_stat = t.test(group1$V1, group2$V1)$statistic  
  return (t_stat)
}

# head(ma2)

#################
# Original set 
#################

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

acq=c("Day 1","Day 2","Day 3","Day 4","Day 5")
var=c(4,5,7,8,9,11,12)

M.ind=as.matrix(ma2[ma2$day%in%acq,var])
# Filtering PRE and CUE
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
# jm[c(1:40),]

res = PCA (jm, scale.unit=TRUE, ind.sup=c(41:455), graph=F) 

# Assigning again the labels of the groups for t statistic calculation
#res$var$coord[,1]
#new_coord <- cbind(res$ind.sup$coord[,1], res$ind.sup$coord[,2])
day <- c()
genotype <- c()
new_coord <- as.data.frame(cbind(res$ind.sup$coord[,1],res$ind.sup$coord[,2]))
new_coord$day <- c(1:415)
new_coord$genotype <- c()

# myset=which(iddaytreat[,2]==tgt[1] & iddaytreat[,3]==acq[1])
myset<-c()
for (j in 1:length(tgt)){
  for (i in 1:length(acq)){
    if (i<5){
      ind=(i-1)*length(tgt)+j
      ind2=i*length(tgt)+j 
      
    }
    myset=which(iddaytreat[,2]==tgt[j] & iddaytreat[,3]==acq[i])
    new_coord [myset,c("day")] <- i
    new_coord [myset,c("genotype")] <- tgt[j]
    new_coord [myset,c("id")] <- substr(iddaytreat[myset,1], 7, 9)    
  }
}

real_t_stat <- f_t_stat (new_coord)

group1 <- subset (new_coord, genotype == "TS")
group2 <- subset (new_coord, genotype == "TSEEEGCG")

#################
# Original set end
#################

#################
# Permutations
#################
# permutations
perm <- 1000

set.seed(333)

t.values = numeric (perm)

for (p in 1:perm) {

  id_group <- subset(ma2, grepl("PRE", ma2$day))
  
  id_group <- id_group [,c(1,2)]
  
  p_genetreat <- sample(id_group$gentreat)
  id_group$gentreat <- p_genetreat
  
  p_ma2 <- merge(id_group, ma2, by="id", all = TRUE)
  
  # I eliminiate the original column with groups labels
  p_ma2 <- p_ma2 [,c(-3)]
  #head (p_ma2, 50)
  #head (ma2, 50)
  
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
  
  
  #####
  # PLOT
  #####
#   ResMMALL=prcomp(MMALL,scale=TRUE)
#   pca2plot <- as.data.frame(ResMMALL$x)
#   row.names(pca2plot)
#   pca2plot$gen_day <- row.names(pca2plot)
#   
#   pca2plot$days <-  as.factor(as.numeric (gsub(".*([0-9]+)$", "\\1", pca2plot$gen_day)))
#   pca2plot$gentreat <-  as.factor(gsub("([A-Z]+).*$", "\\1", pca2plot$gen_day))
#   pca2plot$PC1
#   pca2plot [pca2plot$gentreat == "WT",]
#   
#   pca2plot$gentreat <- factor(pca2plot$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
#                               labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
#   
#   pca2plot <- as.data.frame(ResMMALL$x)
#   row.names(pca2plot)
#   pca2plot$gen_day <- row.names(pca2plot)
#   
#   pca2plot$days <-  as.factor(as.numeric (gsub(".*([0-9]+)$", "\\1", pca2plot$gen_day)))
#   pca2plot$gentreat <-  as.factor(gsub("([A-Z]+).*$", "\\1", pca2plot$gen_day))
#   pca2plot$PC1
#   
#   pca2plot$gentreat <- factor(pca2plot$gentreat , levels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"), 
#                               labels=c("WT", "TS", "WTEE", "TSEE", "WTEGCG", "TSEGCG", "WTEEEGCG", "TSEEEGCG"))
#   
#   pca_medians_acq <- ggplot(pca2plot, aes(x=PC1, y=PC2, colour=gentreat )) + 
#     geom_path (size = 1,show_guide = T) + 
#     #                           geom_path (size = 1,show_guide = F) + 
#     scale_color_manual(values=c("red", "darkgreen", "blue", "lightblue", 
#                                 "magenta", "orange", "gray", "black")) +
#     geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show_guide = T)+
#     #                           geom_text (aes (label=days), vjust=-0.5, hjust=1, size=4, show_guide = F)+
#     theme(legend.key=element_rect(fill=NA)) +
#     labs(title = "PCA of group medians\n", x = "\nPC1 (80% of variance)", y="PC2 (12% of variance)\n") +
#     #                           guides(colour = guide_legend(override.aes = list(size = 10)))+
#     guides(colour = guide_legend(override.aes = list(size = 1)))+
#     theme(legend.key=element_rect(fill=NA))
#   
#   #PLOT_paper
#   pca_medians_acq 
  #####
  # PLOT END
  #####
  
  # All individuals as supplementary points
  acq=c("Day 1","Day 2","Day 3","Day 4","Day 5")
  var=c(4,5,7,8,9,11,12)
  
  M.ind=as.matrix(p_ma2[p_ma2$day%in%acq,var])
  
  # Filtering PRE and CUE
  iddaytreat=as.matrix(p_ma2[p_ma2$day%in%acq,1:3])
  rownames(iddaytreat)=c(1:415)
  tgt=sort(unique(iddaytreat[,2]))
  
  M.med=matrix(0,40,7)
  rnames=c(1:40)
  for (i in 1:length(acq)){
    for (j in 1:length(tgt)){
      ind=(i-1)*length(tgt)+j
      # print (ind)
      M.med[ind,]=apply(M.ind[iddaytreat[,3]==acq[i] & iddaytreat[,2]==tgt[j],],2,median)
      rnames[ind]=paste(tgt[j],acq[i]," ")
    }
  }
  colnames(M.med)=colnames(M.ind)
  rownames(M.med)=rnames
  
  jm=rbind(M.med,M.ind)
  # jm[c(1:40),]
  
  res = PCA (jm, scale.unit=TRUE, ind.sup=c(41:455), graph=F) 

  # Assigning again the labels of the groups for t statistic calculation
  #res$var$coord[,1]
  #new_coord <- cbind(res$ind.sup$coord[,1], res$ind.sup$coord[,2])
  day <- c()
  genotype <- c()
  new_coord <- as.data.frame(cbind(res$ind.sup$coord[,1],res$ind.sup$coord[,2]))
  new_coord$day <- c(1:415)
  new_coord$genotype <- c()
  
  # myset=which(iddaytreat[,2]==tgt[1] & iddaytreat[,3]==acq[1])
  myset<-c()
  for (j in 1:length(tgt)){
    for (i in 1:length(acq)){
      if (i<5){
        ind=(i-1)*length(tgt)+j
        ind2=i*length(tgt)+j 
        
      }
      myset=which(iddaytreat[,2]==tgt[j] & iddaytreat[,3]==acq[i])
      new_coord [myset,c("day")] <- i
      new_coord [myset,c("genotype")] <- tgt[j]
      new_coord [myset,c("id")] <- substr(iddaytreat[myset,1], 7, 9)    
      }
  }
  
  t_s <- f_t_stat (new_coord)
#   t_s_ts_tseeegcg <- f_t_stat (new_coord)
#   t_s_ts_tsee <- f_t_stat (new_coord, "TS", "TSEE")
#   t_s_ts_tsegcg <- f_t_stat (new_coord, "TS", "TSEGCG")
#   t_s_ts_tseeegcg <- f_t_stat (new_coord, "TS", "TSEGCG")
   
  t.values[p] <- t_s
}

# Hacer esta funcion que contenga ya dentro las comparaciones y genere una tabla que tenga
# dentro todas las comparaciones para luego poder hacer subset

t.values <- t.values [order(t.values)]
hist (t.values)
t.values [950]
real_t_stat
(perm - length(t.values [t.values < real_t_stat])) / perm

#t.values_111 <- t.values
#t.values_222 <- t.values
t.values_333 <- t.values

hist (t.values)
t.values_111 [950]
t.values_222 [950]
t.values_333 [950]

(perm - length(t.values_111 [t.values_111 < real_t_stat])) / perm
(perm - length(t.values_222 [t.values_222 < real_t_stat])) / perm
(perm - length(t.values_333 [t.values_333 < real_t_stat])) / perm



real_t_stat
(perm - length(t.values [t.values < real_t_stat])) / perm



