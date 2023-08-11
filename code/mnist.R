#######################################################################
###############  data analysis: mnist data
###############     July 13, 2023
#######################################################################


library(dimRed)
library(Rfast)
library(rARPACK)
library(cluster)
library(ggplot2)


data=read.csv("data/mnist_test.csv", header =F)
table(data[,1])
data=data[c(which(data[,1]==8), which(data[,1]==2),which(data[,1]==6),which(data[,1]==4)),]
data.nolabel = data[,-1]
data.nolabel = data.nolabel[,which(colSums(data.nolabel)!=0)]
info = data[,1]

(n=dim(data.nolabel)[1])
(p=dim(data.nolabel)[2])

dist.mat = dist(data.nolabel)
out = matrix(ncol=7, nrow=7)

feature.num = c(5, 10, 15, 20, 25, 30, 35)
for(num.f in 1:7){
  
  #prop-0.25
  K.mat = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.25)^2)
  eigen.K.mat = eigs(K.mat, k=40)
  embed.data=eigen.K.mat$vectors[,2:(1+feature.num[num.f])] %*% diag(eigen.K.mat$values[2:(1+feature.num[num.f])] )
  out[num.f,1]= summary(silhouette(as.numeric(factor(info)), dmatrix = as.matrix(dist(embed.data)), FUN=mean))$avg.width
  
  #prop-0.75
  K.mat = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.75)^2)
  eigen.K.mat = eigs(K.mat, k=40)
  embed.data=eigen.K.mat$vectors[,2:(1+feature.num[num.f])] %*% diag(eigen.K.mat$values[2:(1+feature.num[num.f])] )
  out[num.f,7]= summary(silhouette(as.numeric(factor(info)), dmatrix = as.matrix(dist(embed.data)), FUN=mean))$avg.width
  
  
  #KPCA
  J = matrix(rep(1/n,n^2), ncol=n)
  K.mat.mean = K.mat - J %*% K.mat - K.mat %*% J + J %*% K.mat %*% J
  eigen.K.mat = eigs(K.mat.mean, k=40)
  embed.data=eigen.K.mat$vectors[,2:(1+feature.num[num.f])] %*% diag(eigen.K.mat$values[2:(1+feature.num[num.f])] )
  out[num.f,2]= summary(silhouette(as.numeric(factor(info)), dmatrix = as.matrix(dist(embed.data)), FUN=mean))$avg.width
  
  
  ##PCA
  eigen.K.mat=svds(as.matrix(data.nolabel),k=40)
  embed.data=eigen.K.mat$u[,1:feature.num[num.f]] %*% diag(eigen.K.mat$d[1:feature.num[num.f]])
  out[num.f,3]= summary(silhouette(as.numeric(factor(info)), dmatrix = as.matrix(dist(embed.data)), FUN=mean))$avg.width
  
  
  ####Laplacian
  D=rowSums(K.mat)
  L = diag(D)-K.mat
  eigen.K.mat=eigen(L)
  embed.data=eigen.K.mat$vectors[,n-(1:feature.num[num.f])] %*% diag(eigen.K.mat$values[n-(1:feature.num[num.f])])
  out[num.f,4]= summary(silhouette(as.numeric(factor(info)), dmatrix = as.matrix(dist(embed.data)), FUN=mean))$avg.width
  
  ##MDS
  dim.red.data = cmdscale(dist(data.nolabel), k=feature.num[num.f])
  out[num.f,5]= summary(silhouette(as.numeric(factor(info)), dmatrix = as.matrix(dist(dim.red.data)), FUN=mean))$avg.width
  
  
  ##DM
  dim.red.data = embed(data.nolabel, "DiffusionMaps", ndim=feature.num[num.f])@data@data
  out[num.f,6]= summary(silhouette(as.numeric(factor(info)), dmatrix = as.matrix(dist(dim.red.data)), FUN=mean))$avg.width
  
  
}


out.data.f = data.frame(tau = c(out[,1],out[,2],out[,3],  out[,4], out[,5], out[,6], out[,7]), 
                        method = rep(c("Prop-0.25","KPCA","PCA","Laplacian","MDS","DM","Prop-0.75"), each=7),
                        n.dim = rep(c(5, 10, 15, 20, 25, 30, 35), 7))

p <- ggplot(out.data.f, aes(n.dim, tau, group = method,colour = method))+
  geom_point(aes(shape = method),size= 3)+
  geom_line(aes(linetype = method),size=1)+
  labs(linetype = "method", shape = "method")+ylab("Average Silhouette Index")
p + theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15,face="bold"),
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14))


