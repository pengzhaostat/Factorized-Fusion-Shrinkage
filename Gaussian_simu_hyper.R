rm(list=ls())
set.seed(2023)
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
source('MP_gibbs_Gaussnetwork_adaptive.R')
source('helper.R')
setwd(file_location)
# 

SMF_errors_all <- list()

num_simulations <- 25
for (sim in 1:num_simulations) {


alpha_values <- seq(0.5, 1, by = 0.05)

#alpha_values = c(alpha_values,1)


SMF_errors <- list()

#--------------------------------Data Generation---------------------

for (alpha in alpha_values){
  

   

  
  n = 20  
  
  p = 10
  
  T = 100
  
  rho = 1
  
  beta = 0 
  
  d = 2  
  
  sigma = 0.3

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = 0.1  # sd for initial state

initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

initial_mean <- matrix(c(0,0,0,0), nrow=2)

initial_Sigma <- sigma_0^2*diag(d)

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




#--------------------------------Model Fitting---------------------


mean_beta_prior = 0    #prior mean of beta, 

sigma_beta_prior = sqrt(2)   #prior sd of beta, 

sigma_X1 = sigma_0;   # initial sd of X[[1]], 

gap =1e-3

min_iter = 1

MF_list =  MP_Gaussain_weighted_adaptive (Y=Y, alpha=alpha, gap =1e-4, max_iter=200,global_prior='Cauthy')



#--------------------------------Computation Efficiency and Estimation Error Comparison---------------------

pred_mean_MF = matrix(rep(0,T*n*p),nrow = T)
res =  matrix(rep(0,T*n*p),nrow = T)

Mix_rmse_sm = rep(0,T)

for (t in 1:T){
  r=1
  for (i in 1:n){
    for (j in 1:p){
        pred_mean_MF[t,r] = MF_list$mean_beta+t(MF_list$Mean_X[[i]][t,])%*% MF_list$Mean_Z[[j]][t,]
        res[t,r] = beta+t(X[[i]][t,])%*% Z[[j]][t,]
        r=r+1
    }
  }
  Mix_rmse_sm[t] = sum((pred_mean_MF[t,]-res[t,])^2/n/T/p)
}

SMF_errors[[as.character(alpha)]] <- sqrt(sum(Mix_rmse_sm))

}

SMF_errors_all[[sim]] <-  SMF_errors

}

SMF_errors_df <- data.frame(matrix(ncol = length(alpha_values), nrow = 0))
colnames(SMF_errors_df) <- as.character(alpha_values)

for (sim in 1:num_simulations) {
  SMF_errors_numeric <- as.numeric(unlist(SMF_errors_all[[sim]]))
  SMF_errors_df[sim, ] <- SMF_errors_numeric
}

boxplot(SMF_errors_df,  main = "Boxplot of RMSE Errors vs fractional power",
        xlab = expression(Fractional~power~alpha), ylab = "RMSE Errors",
        col = "lightblue", border = "black", outline = FALSE)

