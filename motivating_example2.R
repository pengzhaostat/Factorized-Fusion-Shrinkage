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

# trend_generate = function(T,d,sigma_0,p){
#   X = matrix(rep(0,T*d),nrow = T)
#   initial_components = sample(1:2,prob=c(0.5,0.5),size=1, replace=TRUE)
#   initial_mean <- matrix(c(-1,0,1,0), nrow=2)
#   X[1,] = rnorm(n=1, mean = initial_mean[,initial_components], sd = rep(sigma_0,d) )
#   for (t in 2:T) {
#     increment = rep(sample(c(0,1),prob = c(p,1-p),size=1),d)
#     X[t,] <- X[t-1,]+increment*sign(runif(d,-1,1))
#   }
#   return(X)
# }

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
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {
        edge_mat[j,i] = rbinom(n=1, size=1, prob=1/(1+exp(-beta-3*t(X[[i]][t,])%*%X[[j]][t,])))
      }
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}
#--------------------------------Data Generation---------------------

n = 20  

T = 100

rho = 1

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma = 0.2

prob = 0.99

beta = 0


X <- vector("list", n)

for ( i in 1:2){
  X[[i]] = binary_generate(T,d,prob=prob)
}



for ( i in 3:n){
  
  X[[i]] = binary_generate(T,d,prob=1)
  
}

Y = vector("list", T)

p_sparsity = rep(0,T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_binary(X,beta,sigma,t)
  
  p_sparsity[t] = sum(Y[[t]])/n^2
}



MP_list = mix_DN_adaptive_invgamma (Y,mean_beta_prior=0, sigma_beta_prior=sqrt(10), gap =1e-3 )



MF_list = MP_binary_weighted_adaptive (Y, gap =1e-3, max_iter=200,global_prior='Cauthy')



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
  Xm2[[t]] = t(procrustes_r(t(Xm2[[t-1]]),t(Xm2[[t]]))$B.transformed)
}


par(mfrow=c(2,3))
par(mar=rep(3,4))


plot_clus_igraph(Xm[[1]],Y[[1]],c(rep(2,2),rep(1,n-2)),1,1:n,v_shape='square')
mtext('IGLSM', side=2)
plot_clus_igraph(Xm[[75]],Y[[75]],c(rep(2,2),rep(1,n-2)),75,1:n,v_shape='square')
plot_clus_igraph(Xm[[100]],Y[[100]],c(rep(2,2),rep(1,n-2)),100,1:n,v_shape='square')


plot_clus_igraph(Xm2[[1]],Y[[1]],c(rep(2,2),rep(1,n-2)),1,1:n,v_shape='circle')
mtext('FFS', side=2)
plot_clus_igraph(Xm2[[75]],Y[[75]],c(rep(2,2),rep(1,n-2)),75,1:n,v_shape='circle')
plot_clus_igraph(Xm2[[100]],Y[[100]],c(rep(2,2),rep(1,n-2)),100,1:n,v_shape='circle')


Xt = vector("list", T)
for(t in 1:T){
  Xt [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xt [[t]][i,] = X[[i]][t,]
  }
}

heat_matrix = matrix(rep(0,3*n*(T-1)),ncol=T-1)

for (t in 2:T){
  heat_matrix[,t-1]= as.vector(c(sqrt(rowSums((Xm[[t]]-Xm[[t-1]])^2)),sqrt(rowSums((Xm2[[t]]-Xm2[[t-1]])^2)),
                                 sqrt(rowSums((Xt[[t]]-Xt[[t-1]])^2))))
}




library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)

dev.new()

df_plot = heat_matrix%>%
  as_tibble() %>%
  rowid_to_column(var="Index0") %>%
  gather(key="t", value="Value", -1) %>%
  mutate(t=as.numeric(gsub("V","",t)))

df_plot$method = rep(c(rep('IGLSM',n),rep('FSS',n),rep('Truth',n)),T-1)
df_plot$Index = factor(rep(1:n,3*(T-1)))

ggplot(df_plot, aes(Index, t, fill= Value)) + facet_wrap(vars(method))+
  geom_tile()+ scale_fill_gradient(low="white", high="black")+ theme(text = element_text(size = 16))+
  xlab('Row index')+ylab('Transition time')