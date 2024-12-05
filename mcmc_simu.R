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
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_sym_adaptive.R')
source('mcmc_DN_adaptive.R')
source('helper.R')
setwd(file_location)


#--------------------------------Data Generation---------------------

n = 20   # n = 100; 

T = 100   # T = 100;

tau = 0.1 # tau=0.01,0.02,0.05,0.1,0.2,0.3

rho = 0.5

beta = -2 # 

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = sqrt(0.1)  # sd for initial state

initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

initial_mean <- matrix(c(0.5,0,-0.5,0), nrow=2)

initial_Sigma <- sigma_0^2*diag(d)

X <- vector("list", T)


for ( i in 1:n){
  X[[1]] = rbind(X[[1]], rmvnorm(n=1,mean=initial_mean[,initial_components[i]], sigma=initial_Sigma))
}

for(t in 1:(T-1)){X[[t+1]] = X[[1]]}


sig_eps = diag(rep(1,T))

for ( t1 in 1:(T)){
  for (t2 in 1:(T)){
    if (abs(t1-t2)!=0){sig_eps[t1,t2] = rho }
  }
}

tau_Sigma = tau^2*sig_eps


for(i in 1:n){
  eps =  rmvnorm(n=d,mean=rep(0,T), sigma= tau_Sigma)
  for (t in 2:T) {
    X[[t]][i,] <- X[[t-1]][i,] +  eps[,t-1]
  }
}



Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges(X[[t]],beta)
}



#### Priors 

mean_beta_prior = 0    #prior mean of beta_0 

sigma_beta_prior = sqrt(10)   #prior sd of beta_0, 


gap = 1e-4

start_time_MF <- Sys.time()


MF_list =mcmc_DN_adaptive (Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                           gap = gap,max_iter=100)

end_time_MF <- Sys.time()


# New code for max_iter=50
start_time_MF_50 <- Sys.time()

MF_list_50 = mcmc_DN_adaptive(Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                              gap = gap,max_iter=50)
end_time_MF_50 <- Sys.time()

# New code for max_iter=200
start_time_MF_200 <- Sys.time()

MF_list_200 = mcmc_DN_adaptive(Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                               gap = gap,max_iter=200)
end_time_MF_200 <- Sys.time()



start_time_Mix <- Sys.time()


Mix_list = MP_binary_weighted_adaptive(Y = Y,  mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                           gap = gap,first_index = TRUE)

end_time_Mix <- Sys.time()





# pearson correlation coeffcient
auc_mean_mf = rep((n-1)*n/2,0)
auc_mean_mf_50 = rep((n-1)*n/2,0)
auc_mean_mf_200 = rep((n-1)*n/2,0)
auc_mean_mix = rep((n-1)*n/2,0)
auc_res = rep((n-1)*n/2,0)

auc_mf = NULL
auc_mf_50 = NULL
auc_mf_200 = NULL

auc_mix = NULL

for(t in 1:T){
  r=1
  for (i in 1:n){
    for (j in 1:n){
      if (j<i){
        auc_mean_mf[r] = 1/(1+exp(-MF_list$mean_beta-t( MF_list$Mean_X[[t]][i,])%*% MF_list$Mean_X[[t]][j,]))
        auc_mean_mf_50[r] = 1/(1+exp(-MF_list_50$mean_beta-t( MF_list_50$Mean_X[[t]][i,])%*% MF_list_50$Mean_X[[t]][j,]))
        auc_mean_mf_200[r] = 1/(1+exp(-MF_list_200$mean_beta-t( MF_list_200$Mean_X[[t]][i,])%*% MF_list_200$Mean_X[[t]][j,]))  
        auc_mean_mix[r] = 1/(1+exp(-Mix_list$mean_beta-t( Mix_list$Mean_X[[t]][i,])%*% Mix_list$Mean_X[[t]][j,]))
        auc_res[r] = 1/(1+exp(-beta-t(X[[t]][i,])%*% X[[t]][j,]))
        r=r+1
      }
    }
  }
  auc_mf <- c(auc_mf,cor(auc_res, auc_mean_mf,method = "pearson"))
  auc_mf_50 <- c(auc_mf_50,cor(auc_res, auc_mean_mf_50,method = "pearson"))
  auc_mf_200 <- c(auc_mf_200,cor(auc_res, auc_mean_mf_200,method = "pearson"))
  
  auc_mix <- c(auc_mix,cor(auc_res, auc_mean_mix,method = "pearson"))
}

cat('Cycles for MF_50:',MF_list_50$iter,'\n')
cat('Running time for adaptive MF_50:',end_time_MF_50-start_time_MF_50,'\n')
cat('Pearson correlation for MF_50:',mean(auc_mf_50),'\n')

cat('Cycles for MF:',MF_list$iter,'\n')
cat('Running time for adaptive MF:',end_time_MF-start_time_MF,'\n')
cat('Pearson correlation for MF:',mean(auc_mf),'\n')

cat('Cycles for MF_200:',MF_list_200$iter,'\n')
cat('Running time for adaptive MF_200:',end_time_MF_200-start_time_MF_200,'\n')
cat('Pearson correlation for MF_200:',mean(auc_mf_200),'\n')

cat('Cycles for SMF:',Mix_list$iter,'\n')
cat('Running time for SMF:',end_time_Mix-start_time_Mix,'\n')

cat('Pearson correlation for SMF:',mean(auc_mix),'\n')

