simu_change_point = function(n = 20,
                              prob =0.99,
                             d=2, k0=2
){

library(changepoints)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)
library(animation)
library(Rcpp)
# file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(paste(file_location,'/function',sep=""))
source('helper.R')
source('change_point_WB_FFS.R')
sourceCpp('MP_gibbs_network_sym_adaptive_cpp.cpp')
#setwd(file_location)



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


positions_to_edges_binary = function(X,beta,sigma,t){
  n = length(X)
  edge_mat = matrix(rep(0,n*n),nrow = n)
  prob_mat = matrix(rep(0,n*n),nrow = n)
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {
        prob_mat[j,i] = 1/(1+exp(-beta-3*t(X[[i]][t,])%*%X[[j]][t,]))
        edge_mat[j,i] = rbinom(n=1, size=1, prob=1/(1+exp(-beta-t(X[[i]][t,])%*%X[[j]][t,])))
      }
    }
  }
  prob_mat = prob_mat + t(prob_mat)
  edge_mat = edge_mat+ t(edge_mat)
  return(list(edge_mat=edge_mat,prob_mat=prob_mat))
}
#--------------------------------Data Generation---------------------



T = 100


WBS_diff =0
FFS_diff =0
RDP_diff =0


sigma = 0.2



beta = 0

  
  t0 =NULL
  while (is.null(t0)) {
    
    
    
    X <- vector("list", n)
    
    for ( i in 1:k0){
      X[[i]] = binary_generate(T,d,prob=prob)
    }
    
    
    
    for ( i in (k0+1):n){
      
      X[[i]] = binary_generate(T,d,prob=1)
      
    }
    
    
    
    Y = array(0,dim=c(n,n,T))
    prob_mat = array(0,dim=c(n,n,T))
    data_mat=matrix(rep(0,n*(n-1)/2*T),ncol = T)
    
    
    
    
    for (t in 1:T){
      obs_list = positions_to_edges_binary(X,beta,sigma,t)
      obs= obs_list$edge_mat
      Y[,,t] = obs
      prob_mat[,,t] = obs_list$prob_mat
      data_mat[,t] = obs[lower.tri(obs)]
      if(t>1 && sum(abs(prob_mat[,,t]-prob_mat[,,t-1]))!=0)
        t0=c(t0,t-1)
    }
  }
  
  
  M = 100 # number of random intervals for WBS
  intervals = WBS.intervals(M = M, lower = 1, upper = ncol(data_mat))
  delta = 2
  
  
  
  data_mat1 = data_mat[,seq(1,ncol(data_mat),2)]
  data_mat2 = data_mat[,seq(2,ncol(data_mat),2)]
  
  
  
  MF_list = MP_binary_weighted_adaptive (Y, gap =1e-3, max_iter=200, d=d)
  
  FFS_result = FFS.WBS(xhat =MF_list$Mean_X,betahat=MF_list$mean_beta, obs_num = T)
  
  ###WBS using changepoints package
  
  WBS_net_result = WBS.network(data_mat1, data_mat2, 1, ncol(data_mat1), 
                               intervals$Alpha, intervals$Beta, delta = delta)
  
  ###RDP using changepoints package
  WBS_RDP_result = WBS.nonpar.RDPG(data_mat, lowerdiag = TRUE, d,
                                   Alpha = intervals$Alpha, Beta = intervals$Beta, delta=delta)
  
  
  K = length(t0)
  
  l = length(WBS_net_result$S)
  
  WBS_diff = Hausdorff.dist(WBS_net_result$S[order(WBS_net_result$Dval,decreasing = TRUE)[1:min(K,l)]]*2,
                                  t0)
  
  l = length(WBS_RDP_result$S)
  
  
  RDP_diff = Hausdorff.dist(WBS_RDP_result$S[order(WBS_RDP_result$Dval,decreasing = TRUE)[1:min(K,l)]],
                                  t0)
  
  FFS_diff = Hausdorff.dist(FFS_result$S[order(FFS_result$Dval,decreasing = TRUE)[1:K]],
                                  t0)


# cat('mean of WBS is ',mean((WBS_diff)), '; sd of WBS is ',sd(WBS_diff),
#     'correct detection ratio of WBS', sum(WBS_diff==0)/iter_num ,'\n')
# cat('mean of RDP is ',mean((RDP_diff)), '; sd of RDP is ',sd(RDP_diff),
#     'correct detection ratio of RDP', sum(RDP_diff==0)/iter_num ,'\n')
# cat('mean of FFS is ',mean((FFS_diff)), '; sd of FFS is ',sd(FFS_diff),
#     'correct detection ratio of FFS', sum(FFS_diff==0)/iter_num ,'\n')
return(c(WBS_diff,RDP_diff,FFS_diff))
}
