rm(list=ls())
set.seed(2023)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)
library(BayesLogit)
library(animation)
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('mix_DN_adaptive_invgamma.R')
source('helper.R')
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_sym_adaptive.R')
source('mcmc_DN_adaptive.R')
source('helper.R')
setwd(file_location)

auc_node =NULL

auc_non_node = NULL

#for (tau in c(0.01,0.05,0.1,0.5,1)){

for (tau in 1-c(0.01,0.05,0.1,0.5,1)){
  
for (simulation_id in 1:25){

#--------------------------------Data Generation---------------------

n = 20   # n = 100; 

T = 100   # T = 100;

beta = 0 # 

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = sqrt(0.1)  # sd for initial state

initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

initial_mean <- matrix(c(1,0,-1,0), nrow=2)

initial_Sigma <- sigma_0^2*diag(d)



trend_generate = function(T,d,sigma_0,prob){
  X = matrix(rep(0,T*d),nrow = T)
  initial_components = sample(1:2,prob=c(0.5,0.5),size=1, replace=TRUE)
  initial_mean <- matrix(c(1,0,-1,0), nrow=2)
  X[1,] = rnorm(n=1, mean = initial_mean[,initial_components], sd = rep(sigma_0,d) )
  rho = prob
  for (t in 2:T) {
    increment = rep(sample(c(0,0.5),prob = c(rho, 1-rho),size=1),d)
    X[t,] <- X[t-1,]+increment*sign(runif(d,-1,1))
  }
  return(X)
}

X <- vector("list", n)

for ( i in 1:n){
  X[[i]] = trend_generate(T,d,sigma_0,prob=tau)
}


# for ( i in 1:n){
#   X[[1]] = rbind(X[[1]], rmvnorm(n=1,mean=initial_mean[,initial_components[i]], sigma=initial_Sigma))
# }
# 
# for(t in 1:(T-1)){X[[t+1]] = X[[1]]}
# 
# 
# sig_eps = diag(rep(1,T))
# 
# 
# tau_Sigma = tau^2*sig_eps
# 
# 
# for(i in 1:n){
#   eps =  rmvnorm(n=d,mean=rep(0,T), sigma= tau_Sigma)
#   for (t in 2:T) {
#     X[[t]][i,] <- X[[t-1]][i,] +eps[,t]
#   }
# }



Y = vector("list", T)

# for (t in 1:T){
#   Y[[t]] = positions_to_edges(X[[t]],beta)
# }

positions_to_edges_binary = function(X,beta,t){
  n = length(X)
  edge_mat = matrix(rep(0,n*n),nrow = n)
  for (i in 1:n){
    for(j in 1:n){
      prob = 1/(1+exp(-beta-t(X[[i]][t,])%*%X[[j]][t,]))
      edge_mat[i,j] = rbinom(n=1, size=1, prob=prob)
    }
  }
  return(edge_mat)
}


Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_binary(X,beta,t)
}


#### Priors 

mean_beta_prior = 0    #prior mean of beta_0 

sigma_beta_prior = sqrt(10)   #prior sd of beta_0, 


gap = 1e-3

start_time_MF <- Sys.time()


MF_list = mix_DN_adaptive_invgamma (Y,mean_beta_prior=0, sigma_beta_prior=sqrt(10), gap =1e-3 )


end_time_MF <- Sys.time()




start_time_Mix <- Sys.time()


Mix_list = MP_binary_weighted_adaptive (Y, gap =1e-3, max_iter=200,global_prior='Cauthy')


end_time_Mix <- Sys.time()





# pearson correlation coefficient
auc_mean_mf = rep((n-1)*n/2,0)
auc_mean_mix = rep((n-1)*n/2,0)
auc_res = rep((n-1)*n/2,0)

auc_mf = NULL


auc_mix = NULL

for(t in 1:T){
  r=1
  for (i in 1:n){
    for (j in 1:n){
      if (j<i){
        auc_mean_mf[r] = 1/(1+exp(-MF_list$mean_beta-t( MF_list$Mean_X[[t]][i,])%*% MF_list$Mean_X[[t]][j,]))
        auc_mean_mix[r] = 1/(1+exp(-Mix_list$mean_beta-t( Mix_list$Mean_X[[i]][t,])%*% Mix_list$Mean_X[[j]][t,]))
        auc_res[r] = 1/(1+exp(-beta-t(X[[i]][t,])%*% X[[j]][t,]))
        r=r+1
      }
    }
  }
  auc_mf <- c(auc_mf,cor(auc_res, auc_mean_mf,method = "pearson"))

  
  auc_mix <- c(auc_mix,cor(auc_res, auc_mean_mix,method = "pearson"))
}


cat('Cycles for non-nodewise SMF:',MF_list$iter,'\n')
cat('Running time for non-nodewise SMF:',end_time_MF-start_time_MF,'\n')
cat('Pearson correlation for non-nodewise SMF:',mean(auc_mf),'\n')


cat('Cycles for nodewise SMF:',Mix_list$iter,'\n')
cat('Running time for nodewise SMF:',end_time_Mix-start_time_Mix,'\n')
cat('Pearson correlation for nodewise SMF:',mean(auc_mix),'\n')


auc_node = c(auc_node,mean(auc_mix))
auc_non_node = c(auc_non_node,mean(auc_mf))
}
}

auc_node= auc_node[1:100]
auc_non_node = auc_non_node[1:100]


tau_list = as.character(1-rep(c(rep(0.01,5),rep(0.05,5),rep(0.1,5),rep(0.5,5)),2))

 df = data.frame(list(tau = factor(tau_list),
                      method=c(rep('FFS',100),rep('IGLSM',100)),PCC = c(auc_node,auc_non_node)  ))
 

 library(ggplot2)
 ggplot(df, aes(x = tau, y = PCC, fill = method)) +
   geom_boxplot(outlier.shape = NA) +
   labs(x = expression(paste("Stopping probability ",rho)), y = "PCC scores", fill = "method",
        title = "FFS VS IGLSM under sparse transition") +
   theme_bw()+
   theme(plot.title = element_text(hjust = 0.5))
