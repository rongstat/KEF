#######################################################################
###############  data analysis: germline data
###############     July 13, 2023
#######################################################################


library(Seurat)
library(VGAM)
library(TSCAN)
library(SCORPIUS)
library(ggplot2)

#load RNA seq data
load("data/germline-human-female_guo.RData")
#get cell info
info=colnames(data)
label = which(substr(info,1,5)=="M_PGC")

info = data.frame(id=info[label], label=gsub("_","",substr(info[label],7,9)) )
dim(info)
table(info[,2])

data = data[,label]
#normalization and QC
data <- CreateSeuratObject(counts = data, project = "TI", min.cells = 3, min.features = 500)
data <- NormalizeData(data)
n=dim(scale.data.var)[2]

data.names = factor(info[,2])
levels(data.names)=c(3,4,5,1,2)



out  = matrix(ncol=6, nrow=10)
feature.num = seq(2500, 4000, length.out = 10)
for(num.f in 1:10){
  data.s <- FindVariableFeatures(data, selection.method = "vst", nfeatures = feature.num[num.f])
  all.genes <- rownames(data.s)
  data.s <- ScaleData(data.s, features = all.genes)
  scale.data = data.s@assays[["RNA"]]@scale.data
  scale.data.var = scale.data[match(data.s@assays[["RNA"]]@var.features,rownames(scale.data)),]
  
  
  #KEF
  dist.mat = dist(t(scale.data.var))
  K.mat = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
  eigen.K.mat = eigen(K.mat)
  u=order(eigen.K.mat$vectors[,2],decreasing =F)
  ure = order(eigen.K.mat$vectors[,2],decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,1]=max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                     cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  

  #PCA
  u=order(svd(scale.data.var)$v[,1],decreasing =F)
  ure = order(svd(scale.data.var)$v[,1],decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,2]=max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                     cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall") ))
  
  ###Laplacian
  D=rowSums(K.mat)
  L = diag(D)-K.mat
  u.L=order(eigen(L)$vectors[,n-1],decreasing = F)
  u.Lre=order(eigen(L)$vectors[,n-1],decreasing = T)
  u2.L= match(1:n,u.L)
  u2.Lre= match(1:n,u.Lre)
  out[num.f,3]=max(c(cor(as.numeric(data.names), as.numeric(data.names)[u2.L], method="kendall"), 
                     cor(as.numeric(data.names), as.numeric(data.names)[u2.Lre], method="kendall")))
  
  #tscan
  lpsmclust <- exprmclust(scale.data.var)
  out.tscan=TSCANorder(lpsmclust)
  u.tscan=match((out.tscan), colnames(scale.data.var))
  u.tscanre=match(rev(out.tscan), colnames(scale.data.var))
  u2.tscan=match(colnames(scale.data.var),out.tscan)
  out[num.f,4]=max(c(cor(as.numeric(data.names)[!is.na(u2.tscan)], as.numeric(data.names)[u.tscan], method="kendall"), 
                     cor(as.numeric(data.names)[!is.na(u2.tscan)], as.numeric(data.names)[u.tscanre], method="kendall")))
  
 
  #scorpius
  space <- reduce_dimensionality(t(scale.data.var), "spearman")
  space <- reduce_dimensionality(t(scale.data.var), "spearman")
  traj <- infer_trajectory(space)
  u.sc = order(traj$time, decreasing =T)
  u.scre = order(traj$time, decreasing =F)
  u2.sc= match(1:n,u.sc)
  u2.scre= match(1:n,u.scre)
  out[num.f,5]=max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2.sc], method="kendall"),
                     cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2.scre], method="kendall")))
}


out.data.f = data.frame(tau = c(out[,1],out[,2],out[,3],out[,4], out[,5]), 
                        method = rep(c("prop","PCA","SerialRank","TSCAN","SCORPIUS"), each=10),
                        num.feature = rep(seq(2500,4000, length.out = 10), 5))

out.data.f$method = factor(out.data.f$method, levels=c("prop","PCA","SerialRank","TSCAN","SCORPIUS"))
p <- ggplot(out.data.f, aes(num.feature, tau, group = method,colour = method))

p + geom_line(aes(linetype=method))+labs(linetype = "method", shape = "method")+
  geom_point(aes(shape=method),size=5) + ylab("Kendall.tau") + 
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500
############# figures

plot(svd(scale.data.var)$v[,1:2], pch="", col=factor(info[,2]),xlab="PC1", ylab="PC2",font.lab=2, cex.lab = 1.2)
text(svd(scale.data.var)$v[,1:2],  labels=as.character(info[,2]),cex=1, col = as.character(data.names))


plot(eigen(K.mat)$vectors[,c(2,1)], pch="", col=factor(info[,2]),xlab="KEF2", ylab="KEF1", 
     xlim=c(-0.12,0.22), ylim=c(-0.12,0.0), font.lab=2, cex.lab = 1.2)
text(eigen(K.mat)$vectors[,c(2,1)],  labels=as.character(info[,2]),cex=1,col = as.character(data.names))
