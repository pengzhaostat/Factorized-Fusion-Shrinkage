rm(list=ls())
set.seed(2021)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_adaptive.R')
setwd(file_location)



trend_generate = function(T,d,sigma_0){
  X = matrix(rep(0,T*d),nrow = T)
  initial_components = sample(1:2,prob=c(0.5,0.5),size=1, replace=TRUE)
  initial_mean <- matrix(c(1,0,-1,0), nrow=2)
  X[1,] = rnorm(n=1, mean = initial_mean[,initial_components], sd = rep(sigma_0,d) )
  rho = 0.99
  for (t in 2:T) {
    increment = rep(sample(c(0,0.5),prob = c(rho, 1-rho),size=1),d)
    X[t,] <- X[t-1,]+increment*sign(runif(d,-1,1))
  }
  return(X)
}


positions_to_edges_binary = function(X,Z,beta,sigma,t){
  n = length(X)
  p = length(Z)
  edge_mat = matrix(rep(0,n*p),nrow = n)
  for (i in 1:n){
    for(j in 1:p){
      prob = 1/(1+exp(-beta-t(X[[i]][t,])%*%Z[[j]][t,]))
        edge_mat[i,j] = rbinom(n=1, size=1, prob=prob)
    }
  }
  return(edge_mat)
}


evaulate_error = function(model,n,p,response){
  T0 = dim(model$Mean_X[[1]])[1]
  Inner_true = rep(0,n*p)
  Inner_prod = rep(0,n*p)
  k = 1
  for ( i in 1:n){
    for (j in 1:p){
      Inner_prod[k] = 1/(1+exp(model$mean_beta- t(model$Mean_X[[i]][T0,])%*%(model$Mean_Z[[j]][T0,])))
      Inner_true [k] = response[i,j]
      k = k+1  
    }
  }
  error = auc(roc(Inner_true,Inner_prod))
  return(error)
}
#--------------------------------Data Generation---------------------

n = 20 

p = 10

T = 100

rho = 1

beta = 0 

d = 2   

sigma = 0.3


X <- vector("list", n)

Z = vector("list", p)

for ( i in 1:n){
  X[[i]] = trend_generate(T,d,sigma)

}



for ( i in 1:p){
  Z[[i]] = trend_generate(T,d,sigma)
  
}

Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_binary(X,Z,beta,sigma,t)
}







MF_list = MP_binary_weighted_adaptive (Y, gap =0.01, max_iter=200,global_prior='Cauthy',gap_per_iter=1e-3)


svd_pred_list = vector("list", T)

for (t in 1:T){
  svd_list = svd(Y[[t]],nu=d,nv=d)
  svd_pred_list[[t]] = svd_list$u %*% diag(svd_list$d[1:d]) %*% t(svd_list$v)
  svd_pred_list[[t]][svd_pred_list[[t]]<0] =0
  svd_pred_list[[t]][svd_pred_list[[t]]>1] =1
}

svd_pred = rep(0,T*p*n)


pred_mean_MF = rep(0,T*p*n)


res = rep(0,T*p*n)
r=1

for (i in 1:n){
  for (j in 1:p){
      for (t in 1:T){
        pred_mean_MF[r] = 1/(1+exp(MF_list$mean_beta-t(MF_list$Mean_X[[i]][t,])%*% MF_list$Mean_Z[[j]][t,]))
        
        res[r] = 1/(1+exp(-t(X[[i]][t,])%*% Z[[j]][t,]))
        svd_pred[r] = svd_pred_list[[t]][i,j]
        r=r+1
      }
  }
}



cat('PCC Our', cor(res, pred_mean_MF,method = "pearson"))
  
cat('PCC SVD', cor(res, svd_pred,method = "pearson"))
  
  

 

 
 par(mfrow=c(1,2))
 plot(pred_mean_MF,res,xlim = c(0,1),xlab = 'Estimated Prob',ylab = 'True Prob',cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
 title(main = 'FFS',cex.main=1.5)
 lines(seq(0.01,1,by=0.01),seq(0.01,1,by=0.01),col='red',lwd = 2)
 
 
 plot(svd_pred,res,xlab = 'Estimated Prob',ylab = 'True Prob',cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
 title(main = 'SVD',cex.main=1.5)
 lines(seq(0.01,1,by=0.01),seq(0.01,1,by=0.01),col='red',lwd = 2)
 
 

 
