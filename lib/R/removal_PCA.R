#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. May 2015                         ###
#############################################################################
### MWM Paper Silvina frontiers                                           ###
### Modified from Ionas script called PCA_clean for PCA analysis          ### 
#############################################################################

# Calling libraries
library(Hmisc)
library(calibrate)

##Getting HOME directory
home <- Sys.getenv("HOME") 
  
# Loading functions:
source (paste (home, "/git/phecomp/lib/R/plotParamPublication.R", sep=""))

rem_data=spss.get("/Users/jespinosa/20150515_PCA_old_frotiersPaper/data/TS_old_removal.sav")

# rem_data_var = rem_data [ , c(7:10)]

# tbl_stat_mean <-with (df.act_sum, aggregate (cbind (V6), list (index=index, V4=V4), FUN=function (x) c (mean=mean(x), std.error=std.error(x))))
# tbl_stat_mean$mean <- tbl_stat_mean$V18 [,1]
# tbl_stat_mean$std.error <- tbl_stat_mean$V18 [,2]
tbl_stat_median <-with (rem_data, aggregate (cbind (NUMBER.ENTRIES, PERM.TIME, PERCENT.PERM.TIME, LATENCY.TARGET), list (GENTREAT), FUN=function (x) median=median(x)))

# First column are labels, transform to row.names
rownames(tbl_stat_median) <- tbl_stat_median$Group.1
tbl_stat_median <- tbl_stat_median [,-1]
res_rem <- prcomp(tbl_stat_median)
summary(res_rem)
res_rem$x
# png("figures/PCAmedians.png",res=300,width=15,height=15,unit="cm")
dev.off()
plot(res_rem$x[,1],res_rem$x[,2],type="n",main="PCA of group medians",xlab="PA1 (97% of variance)",ylab="PA2 (2% of variance)")
points(res_rem$x[1,1],res_rem$x[1,2],col="red", pch=21, bg="red")
textxy(res_rem$x[1,1],res_rem$x[1,2], row.names(res_rem$x)[1], cex=0.7,col="black")
points(res_rem$x[2,1],res_rem$x[2,2],col="green", pch=21, bg="green")
textxy(res_rem$x[2,1]+2,res_rem$x[2,2]+0.1, row.names(res_rem$x)[2], cex=0.7,col="black")
points(res_rem$x[3,1],res_rem$x[3,2],col="yellow", pch=21, bg="yellow")
textxy(res_rem$x[3,1]-0.7,res_rem$x[3,2]+0.1, row.names(res_rem$x)[3], cex=0.7,col="black")
points(res_rem$x[4,1],res_rem$x[4,2],col="yellow", pch=21, bg="yellow")
textxy(res_rem$x[4,1],res_rem$x[4,2], row.names(res_rem$x)[4], cex=0.7,col="black")
points(res_rem$x[5,1],res_rem$x[5,2],col="magenta", pch=21, bg="magenta")
textxy(res_rem$x[5,1]-0.4,res_rem$x[5,2], row.names(res_rem$x)[5], cex=0.7,col="black")
points(res_rem$x[6,1],res_rem$x[6,2],col="magenta", pch=21, bg="lightblue")
textxy(res_rem$x[6,1]+0.5,res_rem$x[6,2], row.names(res_rem$x)[6], cex=0.7,col="black")
points(res_rem$x[7,1],res_rem$x[7,2],col="blue", pch=21, bg="blue")
textxy(res_rem$x[7,1]-3.5,res_rem$x[7,2]+0.05, row.names(res_rem$x)[7], cex=0.7,col="black")
points(res_rem$x[8,1],res_rem$x[8,2],col="gray", pch=21, bg="gray")
textxy(res_rem$x[8,1]+1,res_rem$x[8,2]-0.05, row.names(res_rem$x)[8], cex=0.7,col="black")

# Plot PCA color by genotype
pca_rem_2plot <- as.data.frame (res_rem$x)
genotype_tt <- as.factor(row.names(res_rem$x))
length (genotype_tt)

g_genotype_tt <- ggplot(pca_rem_2plot, aes(PC1, PC2)) + geom_point(aes(colour=genotype_tt), size=4) +                                                           
  labs (title = "PCA") +
#   scale_color_manual(values=cols, labels=c("1", "2", "3")) +
#   geom_text (aes (label=genotype_tt), hjust=0, vjust=-0.5)
  geom_text (aes (label=genotype_tt), hjust=0.3, vjust=-0.5)
#, position=position_jitter(h=0), alpha = 1)
g_genotype_tt

textxy(res_rem$x[2,1],res_rem$x[2,2], row.names(res_rem$x)[1], cex=0.7,col="green")
lines(res_rem$x[2,1], res_rem$x[2,2],col="black")

textxy(res_rem$x[c(5,13,21,29,37),1],res_rem$x[c(5,13,21,29,37),2],c(1:5),cex=0.7,col="green")
lines(res_rem$x[c(5,13,21,29,37),1],res_rem$x[c(5,13,21,29,37),2],col="green")
textxy(res_rem$x[c(7,15,23,31,39),1],res_rem$x[c(7,15,23,31,39),2],c(1:5),cex=0.7,col="yellow")
lines(res_rem$x[c(7,15,23,31,39),1],res_rem$x[c(7,15,23,31,39),2],col="yellow")
textxy(res_rem$x[c(3,11,19,27,35),1],res_rem$x[c(3,11,19,27,35),2],c(1:5),cex=0.7,col="magenta")
lines(res_rem$x[c(3,11,19,27,35),1],res_rem$x[c(3,11,19,27,35),2],col="magenta")
textxy(res_rem$x[c(6,14,22,30,38),1],res_rem$x[c(6,14,22,30,38),2],c(1:5),cex=0.7,col="lightblue")
lines(res_rem$x[c(6,14,22,30,38),1],res_rem$x[c(6,14,22,30,38),2],col="lightblue")
textxy(res_rem$x[c(2,10,18,26,34),1],res_rem$x[c(2,10,18,26,34),2],c(1:5),cex=0.7,col="blue")
lines(res_rem$x[c(2,10,18,26,34),1],res_rem$x[c(2,10,18,26,34),2],col="blue")
textxy(res_rem$x[c(4,12,20,28,36),1],res_rem$x[c(4,12,20,28,36),2],c(1:5),cex=0.7,col="blue")
lines(res_rem$x[c(4,12,20,28,36),1],res_rem$x[c(4,12,20,28,36),2],col="grey")
legend(x="bottomleft",c("WT","TS","WTEE","TSEE","WTEGCG","TSEGCG","WTEEEGCG","TSEEEGCG"),col=c("red","green","blue","lightblue","magenta","yellow","grey","black"),lty=1,cex=0.6)
altnamesM=c("dist","","latency\ngallindex","speed","percentne","percenter","wishaw","percentperi")
dev.off()
png("figures/PCAmedians_variables.png",res=300,width=15,height=15,unit="cm")
plot(res_rem$rotation[,1],res_rem$rotation[,2],type="n",xlim=c(-1,1),asp=1,main="Variable contributions to principal axes",xlab="PA1 (80% of variance)",ylab="PA2 (11% of variance)")
arrows(0,0,res_rem$rotation[,1],res_rem$rotation[,2],angle=5,length=0.1)
textxy(res_rem$rotation[,1],res_rem$rotation[,2],altnamesM,cex=1)
dev.off()





ResMMALL=prcomp(MMall,scale=TRUE)

MM5[1,]=apply(M[treat=="WT",],2,median)

M=as.matrix(ma2[rem_data$day=="Day 5",c(4,5,7,8,9,10,11,12)])
treat=as.matrix(ma2[ma2$day=="Day 5",2])
M4=as.matrix(ma2[ma2$day=="Day 4",c(4,5,7,8,9,10,11,12)])
treat4=as.matrix(ma2[ma2$day=="Day 4",2])
M3=as.matrix(ma2[ma2$day=="Day 3",c(4,5,7,8,9,10,11,12)])
treat3=as.matrix(ma2[ma2$day=="Day 3",2])
M2=as.matrix(ma2[ma2$day=="Day 2",c(4,5,7,8,9,10,11,12)])
treat2=as.matrix(ma2[ma2$day=="Day 2",2])
M1=as.matrix(ma2[ma2$day=="Day 1",c(4,5,7,8,9,10,11,12)])
treat1=as.matrix(ma2[ma2$day=="Day 1",2])



