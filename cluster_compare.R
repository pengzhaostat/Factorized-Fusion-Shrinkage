rm(list=ls())
set.seed(2021)
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

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_sym_adaptive.R')
source('mix_DN_adaptive.R')
source('mix_DN_adaptive_invgamma.R')
source('helper.R')
setwd(file_location)

n_all =NULL

RI_all = NULL

RI2_all = NULL

RI_FSS_all = NULL

prob_all = NULL

for (prob in c(0.85,0.9,0.95,0.8)){
  
    for (n in c(20,40)){

      for (iter in 1:25){
        
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
        edge_mat[j,i] = rbinom(n=1, size=1, prob=1/(1+exp(-beta*t(X[[i]][t,])%*%X[[j]][t,])))
      }
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}


# n = 20  

T = 100

d = 2  

# prob = 0.9

beta = 9


X <- vector("list", n)

for ( i in 1:n){
  X[[i]] = binary_generate(T,d,prob=prob)
}

Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_binary(X,beta,t)
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


MP_list = mix_DN_adaptive_invgamma (Y,mean_beta_prior=0, sigma_beta_prior=sqrt(10), gap = 0.01 )

MP_list2 = mix_DN_adaptive (Y,mean_beta_prior=0, sigma_beta_prior=sqrt(10), gap = 0.01 )


MF_list = MP_binary_weighted_adaptive (Y, gap =0.01, max_iter=50,global_prior='Cauthy')




pred_mean_MF = rep(T*(n-1)*n,0)
pred_mean_MP = rep(T*(n-1)*n,0)
pred_mean_MP2 = rep(T*(n-1)*n,0)

res = rep(T*(n-1)*n,0)
r=1

for (i in 1:n){
  for (j in 1:n){
    if ( j != i){
      for (t in 1:T){
        pred_mean_MF[r] = 1/(1+exp(MF_list$mean_beta-t(MF_list$Mean_X[[i]][t,])%*% MF_list$Mean_X[[j]][t,]))
        pred_mean_MP[r] = 1/(1+exp(MP_list$mean_beta-t(MP_list$Mean_X[[t]][i,])%*% MP_list$Mean_X[[t]][j,]))
        pred_mean_MP2[r] = 1/(1+exp(MP_list2$mean_beta-t(MP_list2$Mean_X[[t]][i,])%*% MP_list2$Mean_X[[t]][j,]))
        res[r] = 1/(1+exp(-t(X[[i]][t,])%*% X[[j]][t,]))
        r=r+1
      }
    }
  }
}

cat('PCC for inverse gamma',cor(res, pred_mean_MP,method = "pearson"),'\n')
cat('PCC for gamma',cor(res, pred_mean_MP2,method = "pearson"),'\n')
cat('PCC for FSS',cor(res, pred_mean_MF,method = "pearson"),'\n')


# cor(res, pred_mean_MP,method = "pearson")
# cor(res, pred_mean_MP2,method = "pearson")
# cor(res, pred_mean_MF,method = "pearson")

Xm = vector("list", T)
for(t in 1:T){
  Xm [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm [[t]][i,] = MP_list$Mean_X[[t]][i,]
  }
  Xm [[t]] = Xm [[t]]/sqrt(apply(Xm[[t]]^2,1,sum))
}
# for ( t in 2:T){
#   Xm[[t]] = t(procrustes_r(t(Xm[[t-1]]),t(Xm[[t]]))$B.transformed)
# }


Xm2 = vector("list", T)
for(t in 1:T){
  Xm2 [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm2 [[t]][i,] = MP_list2$Mean_X[[t]][i,]
  }
  Xm2 [[t]] = Xm2 [[t]]/sqrt(apply(Xm2[[t]]^2,1,sum))
}


Xm_FFS = vector("list", T)
for(t in 1:T){
  Xm_FFS [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm_FFS [[t]][i,] = MF_list$Mean_X[[i]][t,]
  }
  Xm_FFS [[t]] =Xm_FFS[[t]]/sqrt(apply(Xm_FFS[[t]]^2,1,sum))
}
# for ( t in 2:T){
#   Xm2 [[t]] = t(procrustes_r(t(Xm2[[t-1]]),t(Xm2[[t]]))$B.transformed)
# }

Zm = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
km = kmeans(Xm[[t]],centers=2^d,nstart = 100)
Zm[t,] = km$cluster
}

Zm2 = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
  km = kmeans(Xm2[[t]],centers=2^d,nstart = 100)
  Zm2 [t,] = km$cluster
}

Zm_FFS = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
  km = kmeans(Xm_FFS[[t]],centers=2^d,nstart = 100)
  Zm_FFS [t,] = km$cluster
}

RI = rep(0,T)
RI2 = rep(0,T)
RI_FFS = rep(0,T)

for(t in 1:T){
RI[t] = rand.index( Zm[t,], Zt[t,] )
RI2[t]= rand.index( Zm2[t,], Zt[t,] )
RI_FFS[t] =  rand.index( Zm_FFS[t,], Zt[t,] )
}

cat('rand index for inverse gamma',mean(RI),'\n')
cat('rand index for gamma',mean(RI2),'\n')
cat('rand index for FSS',mean(RI_FFS),'\n')

RI_all = c(RI_all,mean(RI))
RI2_all = c(RI2_all,mean(RI2))
RI_FSS_all = c(RI_FSS_all,mean(RI_FFS))
n_all = c(n_all,mean(n))
prob_all = c(prob_all,mean(prob))

      }
    }
}

df_plot = data.frame(list(n=factor(n_all),prob=factor(prob_all), RI = c(RI_FSS_all,RI_all,RI2_all),
                           varaince_prior = c(rep('FFS',200),rep('Inv_Gamma',200),
                                              rep('Gamma',200))))

df_plot = data.frame(list(n=factor(n_all),prob=factor(prob_all), RI = c(RI_FSS_all,RI_all),
                          varaince_prior = c(rep('FFS',200),rep('Inv_Gamma',200) )))

# #
# ggplot(data=df_plot, aes(x=prob, y=RI, group=varaince_prior)) + theme(legend.position="top")+
#    geom_boxplot(size=1.5,aes(color = varaince_prior, shape = varaince_prior)) + ylab("RI")+
#    ggtitle("RI Comparison") + theme(plot.title = element_text(hjust = 0.5))


df_plot$n = as.character(df_plot$n)

df_plot$n[df_plot$n=='20']='n=20'
df_plot$n[df_plot$n=='40']='n=40'
df_plot$n = factor(df_plot$n)

ggplot(df_plot, aes(prob, y=RI, fill=varaince_prior)) +
  geom_boxplot(aes(color = varaince_prior, shape = varaince_prior))+ facet_wrap(vars(factor(n)),scale='free')+ theme(axis.text.x = element_text(colour = 'black', size = 15),
                                                                                               axis.title.x = element_text(size = 15,
                                                                                                                           hjust = 0.5, vjust = 0.2)) +
  theme(axis.text.y = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(size = 15,
                                    hjust = 0.5, vjust = 0.2)) +
  theme(legend.text=element_text(size=15),legend.title=element_text(size=15),
        legend.key.size = unit(1, 'cm'))+
  theme(strip.text = element_text(size=15))



# library(igraph)
# 
# A= matrix(rep(0,25),nrow = 5)
# 
# A[1,3] = 1
# A[2,3] = 1
# A[4,3] =1
# A[5,3] =1
# 
# A= t(A)+A
# 
# 
# G = graph_from_adjacency_matrix(A,mode = "undirected")
# 
# plot(G,layout = layout.random,vertex.size = 20)
#  p2=  ggplot(data=df_change,aes(x=time,y=change))+ geom_line(size=1.5,color='red')
#  
#  grid.arrange(p1, p2, nrow = 2,heights=c(0.6,0.4))
# 
# change =rep(0,T)
# for (t in 2:T){
#   change[t] = sqrt(sum((Xt[[t]]-Xt[[t-1]])^2))
# }
# df_change =  data.frame(list(time=1:100, change = change/max(change))) 
# 
# RI2 =(RI2-min(RI2))/(max(RI2)-min(RI2))
# RI =(RI-min(RI))/(max(RI)-min(RI))
# change =(change-min(change))/(max(change)-min(change))
# 
# df_plot = data.frame(list(time=rep(1:100,3), RI = c(RI2,RI,change),
#                           varaince_prior = c(rep('FSS',100),rep('Inv_Gamma',100),rep('change',100))))
# 
#  ggplot(data=df_plot, aes(x=time, y=RI, group=varaince_prior)) + theme(legend.position="right")+
#    geom_line(size=1.5,aes(color = varaince_prior)) + ylab("RI")+
#    ggtitle("RI Comparison") + theme(plot.title = element_text(hjust = 0.5))
#  
#  
#  ggplot(data=df_plot, aes(x=varaince_prior, y=RI)) + 
#    geom_boxplot()