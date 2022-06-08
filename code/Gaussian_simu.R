rm(list=ls())
set.seed(2021)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(genlasso)
library(MSGLasso)


file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_Gaussnetwork_adaptive.R')
setwd(file_location)




trend_generate = function(T,d){
  X = matrix(rep(0,T*d),nrow = T)
  initial_mean <- rnorm(d)
  X[1,] = initial_mean
  for (t in 2:T) {
    increment = rep(sample(c(0,1),prob = c(0.99,0.01),size=1),d)
      X[t,] <- X[t-1,]+increment*sign(runif(d,-1,1))
  }
  return(X)
}




positions_to_edges_Gaussian = function(X,Z,beta,sigma,t){
  n = length(X)
  edge_mat = matrix(rep(0,n*p),nrow = n)
  for (i in 1:n){
    for(j in 1:p){
      { 
        edge_mat[i,j] = rnorm(n=1,mean=beta+t(X[[i]][t,])%*%Z[[j]][t,],sd=sigma)}
    }
  }
  return(edge_mat)
}
#--------------------------------Data Generation---------------------

evaulate_error = function(model,n,p,response){
  T0 = dim(model$Mean_X[[1]])[1]
  Inner_true = rep(0,n*p)
  Inner_prod = rep(0,n*p)
  k = 1
  for ( i in 1:n){
    for (j in 1:p){
        Inner_prod[k] = (model$mean_beta+ t(model$Mean_X[[i]][T0,])%*%(model$Mean_Z[[j]][T0,]))
        Inner_true [k] = response[i,j]
      k = k+1  
    }
  }
  error = sum((Inner_prod-Inner_true)^2/n/p/T0)
  return(error)
}


n = 20  

p = 10

T = 100

rho = 1

beta = 0 

d = 2  

sigma = 0.3



initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)


X <- vector("list", n)

for ( i in 1:n){
  X[[i]] = trend_generate(T,d)
}
Z = vector("list", p)


for ( i in 1:p){
  Z[[i]] = trend_generate(T,d)
}

Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_Gaussian(X,Z,beta,sigma,t)
}


 MF_list = MP_Gaussain_weighted_adaptive (Y=Y, gap =1e-3, max_iter=100,global_prior='Cauthy')
 

svd_pred_list = vector("list", T)

for (t in 1:T){
  svd_list = svd(Y[[t]],nu=d,nv=d)
  svd_pred_list[[t]] = svd_list$u %*% diag(svd_list$d[1:d]) %*% t(svd_list$v)
}


Inner_true = matrix(rep(0,n*p*T),nrow= n*p)
Inner_prod = matrix(rep(0,n*p*T),nrow= n*p)
svd_pred = matrix(rep(0,n*p*T),nrow= n*p)


k = 1
for ( i in 1:n){
  for (j in 1:p){

     for (t in 1:T){
  Inner_prod[k,t] = (MF_list$mean_beta+ t(MF_list$Mean_X[[i]][t,])%*%(MF_list$Mean_Z[[j]][t,]))
  Inner_true [k,t] = (t(X[[i]][t,])%*%(Z[[j]][t,]))
  svd_pred[k,t] =  svd_pred_list[[t]][i,j]
     }
  k = k+1  
 }
}

Y_reshape = matrix(rep(0,n*p*T), ncol = T)
k=1
for (i in 1:n){
  for(j in 1:p){
    for(t in 1:T){
      Y_reshape[k,t] = Y[[t]][i,j]
    }
    k = k+1
  }
}


M_reshape = matrix(rep(0,n*p*T), ncol = T)
for(k in 1:(n*p)){
  trend_f = trendfilter(Y_reshape[k,])
  lambda_min = cv.trendfilter(trend_f,verbose = FALSE)$lambda.min
  M_reshape[k,] = trend_f$beta[,which(trend_f$lambda == lambda_min)[1]]
}


Mfused_svd_list = vector("list", T)

for (t in 1:T){
  M_t = matrix(M_reshape[,t],nrow =n,byrow = TRUE)
  svd_list = svd(M_t,nu=d,nv=d)
  Mfused_svd_list[[t]] = svd_list$u %*% diag(svd_list$d[1:d]) %*% t(svd_list$v)
}

Mfused_svd_pred = matrix(rep(0,n*p*T),nrow= n*p)

k = 1
for ( i in 1:n){
  for (j in 1:p){
    for (t in 1:T){
      Mfused_svd_pred[k,t] = Mfused_svd_list[[t]][i,j]
    }
    k = k+1
  }
}


cat('RMSE for FFS',sqrt(sum((Inner_prod-Inner_true)^2/n/p/T)),'\n')


cat('RMSE for SVD',sqrt(sum((svd_pred-Inner_true)^2/n/p/T)),'\n')


cat('RMSE for Flasso1',sqrt(sum((M_reshape-Inner_true)^2/n/p/T)),'\n')


cat('RMSE for Flasso2',sqrt(sum((Mfused_svd_pred-Inner_true)^2/n/p/T)),'\n')
