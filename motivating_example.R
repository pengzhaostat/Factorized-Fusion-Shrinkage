rm(list=ls())
set.seed(1234)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)
library(animation)


file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_sym_adaptive.R')
source('mix_DN_adaptive_invgamma.R')
source('helper.R')
setwd(file_location)



trend_generate = function(T,d,sigma_0){
  X = matrix(rep(0,T*d),nrow = T)
  initial_components = sample(1:2,prob=c(0.5,0.5),size=1, replace=TRUE)
  initial_mean <- matrix(c(-1,0,1,0), nrow=2)
  X[1,] = rnorm(n=1, mean = initial_mean[,initial_components], sd = rep(sigma_0,d) )
  for (t in 2:T) {
    increment = rep(sample(c(0,0.5),prob = c(1,0),size=1),d)
    X[t,] <- X[t-1,]+increment*sign(runif(d,-1,1))
  }
  return(X)
}



positions_to_edges_binary = function(X,beta,sigma,t){
  n = length(X)
  edge_mat = matrix(rep(0,n*n),nrow = n)
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {
        edge_mat[j,i] = rbinom(n=1, size=1, prob=1/(1+exp(-beta-t(X[[i]][t,])%*%X[[j]][t,])))
      }
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}
#--------------------------------Data Generation---------------------

n = 10  

T = 100

rho = 1

beta = 0 

d = 2  

sigma = 0.2


X <- vector("list", n)

for ( i in 1:n){
  X[[i]] = trend_generate(T,d,sigma)

}

Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_binary(X,beta,sigma,t)
}





MP_list = mix_DN_adaptive_invgamma (Y,mean_beta_prior=0, sigma_beta_prior=sqrt(10), gap = 0.01 )



MF_list = MP_binary_weighted_adaptive (Y, gap =0.01, max_iter=200,global_prior='Cauthy')




pred_mean_MF = rep(T*(n-1)*n,0)
pred_mean_MP = rep(T*(n-1)*n,0)
res = rep(T*(n-1)*n,0)
r=1

for (i in 1:n){
  for (j in 1:n){
    if ( j != i){
      for (t in 1:T){
        pred_mean_MF[r] = 1/(1+exp(MF_list$mean_beta-t(MF_list$Mean_X[[i]][t,])%*% MF_list$Mean_X[[j]][t,]))
        pred_mean_MP[r] = 1/(1+exp(MP_list$mean_beta-t(MP_list$Mean_X[[t]][i,])%*% MP_list$Mean_X[[t]][j,]))
        res[r] = 1/(1+exp(-t(X[[i]][t,])%*% X[[j]][t,]))
        r=r+1
      }
    }
  }
}


  cor(res, pred_mean_MP,method = "pearson")
  cor(res, pred_mean_MF,method = "pearson")

 Xm = vector("list", T)
 for(t in 1:T){
   Xm [[t]] =matrix(rep(0,n*d),nrow=n)
   for (i in 1:n){
     Xm [[t]][i,] = MP_list$Mean_X[[t]][i,]
   }
 }
 for ( t in 2:T){
   Xm[[t]] = t(procrustes_r(t(Xm[[t-1]]),t(Xm[[t]]))$B.transformed)
 }
 
 

 Xm2 = vector("list", T)
 for(t in 1:T){
   Xm2 [[t]] =matrix(rep(0,n*d),nrow=n)
   for (i in 1:n){
     Xm2 [[t]][i,] = MF_list$Mean_X[[i]][t,]
   }
 }
 for ( t in 2:T){
   Xm2 [[t]] = t(procrustes_r(t(Xm2[[t-1]]),t(Xm2[[t]]))$B.transformed)
 }
 
 par(mfrow=c(2,3))
 par(mar=rep(2,4))
 
 plot_clus_igraph(Xm[[1]],Y[[1]],rep(1,n),1,1:n)
 plot_clus_igraph(Xm[[50]],Y[[50]],rep(1,n),50,1:n)
 plot_clus_igraph(Xm[[100]],Y[[100]],rep(1,n),100,1:n)
 plot_clus_igraph(Xm2[[1]],Y[[1]],rep(1,n),1,1:n)
 plot_clus_igraph(Xm2[[50]],Y[[50]],rep(1,n),50,1:n)
 plot_clus_igraph(Xm2[[100]],Y[[100]],rep(1,n),100,1:n)
