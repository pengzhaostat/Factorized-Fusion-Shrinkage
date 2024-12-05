rm(list=ls())
set.seed(2024)
library(Rcpp)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)
library(animation)
library(fossil)
library(ggplot2)
library(gridExtra)
library(grid)
library(dynsbm)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
# source('MP_gibbs_multi_Sigma_adaptive.R')
# source('MP_gibbs_network_sym_adaptive.R')
# source('mix_DN_adaptive.R')
#source('mix_DN_adaptive_invgamma.R')

source('helper.R')
#setwd(file_location)
sourceCpp('mix_DN_adaptive_inv_gamma_new.cpp')
sourceCpp('MP_gibbs_network_sym_adaptive_cpp.cpp')

n = 50
prob =0.9
d=2
k0 = 5
beta = 0.5

        
binary_generate = function(T,d,prob){
  X = matrix(rep(0,T*d),nrow = T)
  X[1,] <- rbinom(n=d,size=1,prob=0.5)
  for (t in 2:T) {
    increment = rbinom(n=1,size=1,prob=(prob-1/(2^d))/(1-1/(2^d)))
      if (increment == 1){
        X[t,] <- X[t-1,]
      } else{
        X[t,] <- rbinom(n=d,size=1,prob=0.5)
      }
  }
  X[X==0] = -1
  return(X)
}

positions_to_edges_binary = function(X,beta,t){
  n = length(X)
  edge_mat = matrix(rep(0,n*n),nrow = n)
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {
        edge_mat[j,i] = rbinom(n=1, size=1, prob=1/(1+exp(2-beta*t(X[[i]][t,])%*%X[[j]][t,])))
      }
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}


# n = 20  

T = 100


X <- vector("list", n)

for ( i in 1:k0){
  X[[i]] = binary_generate(T,d,prob=prob)
}

for ( i in (k0+1):n){
  X[[i]] = binary_generate(T,d,prob=1)
}
Ym = vector("list", T)

Y = array(0,dim=c(n,n,T))

for (t in 1:T){
  Y[,,t] = positions_to_edges_binary(X,beta,t)
  
  Ym[[t]] = Y[,,t]
  

}



Xt = vector("list", T)
for(t in 1:T){
  Xt [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xt [[t]][i,] = X[[i]][t,]
  }
}

Zt = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
  km = kmeans(Xt[[t]],centers=dim(unique(Xt[[t]]))[1])
  Zt[t,] = km$cluster
}


MP_list = mix_DN_adaptive_invgamma (Ym, gap = 1e-3, max_iter=200 )


FFS_list = MP_binary_weighted_adaptive (Y, gap =1e-3, max_iter=200)

Y_swapped <- aperm(Y, c(3, 2, 1))

list.dynsbm <- select.dynsbm(Y_swapped,Qmin=2^d,Qmax=2^d,edge.type="binary")





Zm = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
  Xm  =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm [i,] = MP_list$Mean_X[[t]][i,]
  }
  Xm  = Xm /sqrt(apply(Xm^2,1,sum))
  km = kmeans(Xm,centers=2^d,nstart = 100)
  Zm[t,] = km$cluster
}




Zm_FFS = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
  Xm_FFS =  FFS_list$Mean_X[,,t]/sqrt(apply(FFS_list$Mean_X[,,t]^2,1,sum))
  km = kmeans(Xm_FFS,centers=2^d,nstart = 100)
  Zm_FFS [t,] = km$cluster
}




RI_IG = rep(0,T)
RI_dynsbm = rep(0,T)
RI_FFS = rep(0,T)

for(t in 1:T){
RI_IG[t] = rand.index( Zm[t,], Zt[t,] )
RI_dynsbm[t]= rand.index( list.dynsbm[[1]]$membership[,t], Zt[t,] )
RI_FFS[t] =  rand.index( Zm_FFS[t,], Zt[t,] )
}


c(mean(RI_IG),mean(RI_dynsbm),mean(RI_FFS))


