rm(list=ls())
set.seed(2021)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(rTensor)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_tensor_adaptive.R')
setwd(file_location)



T = 100

m = 3

n = 10

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma = 0.3


trend_generate_tensor = function(T,n,d,m){
  prob = 0.95
  X = rand_tensor(c(n,d,m,T))
  for (t in 2:T) {
    increment = rep(sample(c(0,0.25),prob = c(prob, 1-prob),size=1),d)
    for(k in 1:m){
    X[,,k,t] <- X[,,k,t-1]+increment*sign(runif(d,-1,1))
    }
  }
  return(X)
}



positions_to_edges_tensor = function(X,sigma){
  n = dim(X)[1]
  m = dim(X)[3]
  T = dim(X)[4]
  Y = array(rep(0,),dim=c(rep(n,m),T))
  tensor_mean = array(rep(0,),dim=c(rep(n,m),T))
  for ( t in 1:T){
    for (d0 in 1:d){
      a = X[,d0,1,t]@data
  for (i in 2:m){
    a = outer(a,X[,d0,i,t]@data)
  }
      tensor_mean[,,,t] = tensor_mean[,,,t] +a
    }
  Y[,,,t] = rnorm(n=n^m, mean=tensor_mean[,,,t],sd=sigma)
  }
  return(list(Y=Y,tensor_mean=tensor_mean))
}


X = trend_generate_tensor(T,n,d,m)
  

Y_list = positions_to_edges_tensor(X,sigma)

Y = Y_list$Y

tensor_mean = Y_list$tensor_mean

MF_list = MP_Gaussain_tensor_adaptive (Y, gap =1e-3, max_iter=100,sigma=sigma)


MF_error = sqrt(sum((as.vector(MF_list$pred_mean)-as.vector(tensor_mean))^2/n^3/T))

MF_error


CP_error = 0
for (t in 1:T){ 
  CPD <- cp(as.tensor(Y[,,,t]),num_components =d,max_iter = 500, tol = 1e-3)
  CP_error = CP_error+ sum((as.vector(CPD$est@data-tensor_mean[,,,t]))^2)
}
CP_error = sqrt(CP_error/(n^3)/T)

CP_error