#######################################################################
###############  data analysis: cell cycle data
###############     July 13, 2023
#######################################################################


library(lle)
library(dimRed)
library(ggplot2)
library(Seurat)

setwd("data/E-MTAB-2805")
#read data
g2m=read.table("G2M_singlecells_counts.txt", header =T)
s=read.table("S_singlecells_counts.txt", header =T)
g1=read.table("G1_singlecells_counts.txt", header =T)
dim(g2m)
dim(s)
dim(g1)

gene = g1[,1:4]
g1=g1[,-(1:4)]
s=s[,-(1:4)]
g2m=g2m[,-(1:4)]
combined = cbind(g1, s, g2m)

#remove ERCCs
combined = combined[-which(substr(gene[,1], 1,4)=="ERCC"),]
gene= gene[-which(substr(gene[,1], 1,4)=="ERCC"),]
dim(combined)

data=combined
colnames(data)

info = c(rep(c("G1","S","G2M"), each=96))
table(info)

#normalization and QC
data <- CreateSeuratObject(counts = data, project = "cycle", min.cells = 3, min.features = 20)
data <- NormalizeData(data)

n=288
data.names = factor(info)
levels(data.names)=c(1,3,2)

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

out = matrix(ncol=6, nrow=10)


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
  u1=eigen.K.mat$vectors[,2] * eigen.K.mat$values[2]
  u2=eigen.K.mat$vectors[,3] * eigen.K.mat$values[3]
  
  u=order(atan(u1/u2),decreasing =F)
  ure = order(atan(u1/u2),decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,1]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  for(i in 1:(n-1)){
    u2= match(shifter(1:n, n =i),u)
    u2re = match(shifter(1:n, n =2),ure)
    temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
    if(temp >  out[num.f,1]){
      out[num.f,1]=temp
    }
  }
  
  #PCA
  u1=svd(scale.data.var)$v[,1] * svd(scale.data.var)$d[1]
  u2=svd(scale.data.var)$v[,2] * svd(scale.data.var)$d[2]
  u=order(atan(u1/u2),decreasing =F)
  ure = order(atan(u1/u2),decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,2]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                      cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  for(i in 1:n-1){
    u2= match(shifter(1:n, n =i),u)
    u2re = match(shifter(1:n, n =2),ure)
    temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
    if(temp >  out[num.f,2]){
      out[num.f,2]=temp
    }
  }
  
  ###Laplacian
  D=rowSums(K.mat)
  L = diag(D)-K.mat
  u1=eigen(L)$vectors[,n-1] * eigen(L)$values[n-1]
  u2=eigen(L)$vectors[,n-2] * eigen(L)$values[n-2]
  u=order(atan(u1/u2),decreasing =F)
  ure = order(atan(u1/u2),decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,3]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                      cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  for(i in 1:n-1){
    u2= match(shifter(1:n, n =i),u)
    u2re = match(shifter(1:n, n =2),ure)
    temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
    if(temp >  out[num.f,3]){
      out[num.f,3]=temp
    }
  }
  
  #MDS
  dim.red.data = cmdscale(dist(t(scale.data.var)), k=2)
  u1=dim.red.data[,1]
  u2=dim.red.data[,2]
  u=order(atan(u1/u2),decreasing =F)
  ure = order(atan(u1/u2),decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,4]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                      cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  for(i in 1:n-1){
    u2= match(shifter(1:n, n =i),u)
    u2re = match(shifter(1:n, n =2),ure)
    temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
    if(temp >  out[num.f,4]){
      out[num.f,4]=temp
    }
  }
  
  
  #LLE
  lle.data = lle(t(scale.data.var), m=2, k=40, reg=2)
  dim.red.data =  lle.data$Y
  u1=dim.red.data[,1]
  u2=dim.red.data[,2]
  u=order(atan(u1/u2),decreasing =F)
  ure = order(atan(u1/u2),decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,5]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                      cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  for(i in 1:n-1){
    u2= match(shifter(1:n, n =i),u)
    u2re = match(shifter(1:n, n =2),ure)
    temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
    if(temp >  out[num.f,5]){
      out[num.f,5]=temp
    }
  }
  
  #DM
  dim.red.data = embed(t(scale.data.var), "DiffusionMaps")@data@data
  u1=dim.red.data[,1]
  u2=dim.red.data[,2]
  u=order(atan(u1/u2),decreasing =F)
  ure = order(atan(u1/u2),decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  out[num.f,6]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                      cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  for(i in 1:n-1){
    u2= match(shifter(1:n, n =i),u)
    u2re = match(shifter(1:n, n =2),ure)
    temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
    if(temp >  out[num.f,6]){
      out[num.f,6]=temp
    }
  }
  
  
  
}


out.data.f = data.frame(tau = c(out[,1],out[,2],out[,3],out[,4], out[,5], out[,6]), 
                        method = rep(c("prop","PCA","Laplacian","MDS","LLE","DM"), each=10),
                        num.feature = rep(seq(2500,4000, length.out = 10), 6))

out.data.f$method = factor(out.data.f$method, levels=c("prop","PCA","Laplacian","MDS","LLE","DM"))

p <- ggplot(out.data.f, aes(num.feature, tau, group = method,colour = method))

p + geom_line(aes(linetype=method))+labs(linetype = "method", shape = "method")+
  geom_point(aes(shape=method),size=5) + ylab("Kendall.tau") + 
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500




