simu_cluster = function(n = 100,
                             prob =0.99,
                             beta=0.5,
                          k0=2
){
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



# source('MP_gibbs_multi_Sigma_adaptive.R')
# source('MP_gibbs_network_sym_adaptive.R')
# source('mix_DN_adaptive.R')
 
  
sourceCpp('mix_DN_adaptive_inv_gamma_new.cpp')
sourceCpp('MP_gibbs_network_sym_adaptive_cpp.cpp')
source('helper.R')
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


d =2

ind_cluster_num = 0


while (ind_cluster_num !=T) {
  
ind_cluster_num = 0

X <- vector("list", n)

for ( i in 1:k0){
  X[[i]] = binary_generate(T,d,prob=prob)
}

for ( i in (k0+1):n){
  X[[i]] = binary_generate(T,d,prob=1)
}

Xt = vector("list", T)
for(t in 1:T){
  Xt [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xt [[t]][i,] = X[[i]][t,]
  }
  if (dim(unique(Xt[[t]]))[1] ==2^d){
    ind_cluster_num = ind_cluster_num+1
  }
  
}



}

Zt = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
  km = kmeans(Xt[[t]],centers=dim(unique(Xt[[t]]))[1])
  Zt[t,] = km$cluster
}


Ym = vector("list", T)

Y = array(0,dim=c(n,n,T))

for (t in 1:T){
  Y[,,t] = positions_to_edges_binary(X,beta,t)
  Ym[[t]] = Y[,,t]
}



#MP_list2 = mix_DN_adaptive (Y,mean_beta_prior=0, sigma_beta_prior=sqrt(10), gap = 0.01 )


#FFS_list = MP_binary_weighted_adaptive (Y,d=d, gap =1e-3, max_iter=200)

#Y_swapped <- aperm(Y, c(3, 2, 1))

#list.dynsbm <- select.dynsbm(Y_swapped,Qmin=2^d,Qmax=2^d,edge.type="binary")


start_time_MF <- Sys.time()
MP_list = mix_DN_adaptive_invgamma (Ym ,d=d, gap =1e-3, max_iter=200)
end_time_MF <- Sys.time()
time_MP_gibbs <- as.numeric(difftime(end_time_MF, start_time_MF, units = "secs"))

# Record start and end time for MP_binary_weighted_adaptive
start_time_MF <- Sys.time()
FFS_list <- MP_binary_weighted_adaptive(Y, d = d, gap = 1e-3, max_iter = 200)
end_time_MF <- Sys.time()
time_FFS <- as.numeric(difftime(end_time_MF, start_time_MF, units = "secs"))


Y_swapped <- aperm(Y, c(3, 2, 1))

# Record start and end time for select.dynsbm
start_time_MF <- Sys.time()
list.dynsbm <- select.dynsbm(Y_swapped, Qmin = 2^d, Qmax = 2^d, edge.type = "binary")
end_time_MF <- Sys.time()
time_dynsbm <- as.numeric(difftime(end_time_MF, start_time_MF, units = "secs"))



# pred_mean_FFS = rep(T*(n-1)*n,0)
# pred_mean_MP = rep(T*(n-1)*n,0)
# pred_mean_dynsbm = rep(T*(n-1)*n,0)
# 
# res = rep(T*(n-1)*n,0)
# r=1
# 
# for (i in 1:n){
#   for (j in 1:n){
#     if ( j != i){
#       for (t in 1:T){
#         pred_mean_FFS[r] = 1/(1+exp(-FFS_list$mean_beta-t(FFS_list$Mean_X[i,,t])%*% FFS_list$Mean_X[j,,t]))
#         pred_mean_MP[r] = 1/(1+exp(-MP_list$mean_beta-t(MP_list$Mean_X[[t]][i,])%*% MP_list$Mean_X[[t]][j,]))
#         member_vec= list.dynsbm[[1]]$membership[,t]
#         pred_mean_dynsbm[r] = list.dynsbm[[1]]$beta[t,member_vec[i],member_vec[j]]
#         res[r] = 1/(1+exp(-t(X[[i]][t,])%*% X[[j]][t,]))
#         r=r+1
#       }
#     }
#   }
# }
# 
# cat('PCC for inverse gamma',cor(res, pred_mean_MP,method = "pearson"),'\n')
# cat('PCC for dynsbm',cor(res, pred_mean_dynsbm,method = "pearson"),'\n')
# cat('PCC for FSS',cor(res, pred_mean_FFS,method = "pearson"),'\n')


# cor(res, pred_mean_MP,method = "pearson")
# cor(res, pred_mean_MP2,method = "pearson")
# cor(res, pred_mean_MF,method = "pearson")

Zm = matrix(rep(0,n*T),nrow=T)
Xm = matrix(rep(0,n*d),nrow=n)
for (t in 1:T){
  for(i in 1:n){
    Xm[i,] = MP_list$Mean_X[[t]][i,]/(sqrt(sum(MP_list$Mean_X[[t]][i,]^2)))
  }
  km = kmeans(Xm,centers=dim(unique(Xt[[t]]))[1],nstart = 100)
  Zm[t,] = km$cluster
}



Zm_FFS = matrix(rep(0,n*T),nrow=T)
for (t in 1:T){
  Xm_FFS =  FFS_list$Mean_X[,,t]/(sqrt(apply(FFS_list$Mean_X[,,t]^2,1,sum)))
  km = kmeans(Xm_FFS,centers=dim(unique(Xt[[t]]))[1],nstart = 100)
  Zm_FFS [t,] = km$cluster
}




RI = rep(0,T)
RI_dynsbm = rep(0,T)
RI_FFS = rep(0,T)

for(t in 1:T){
RI[t] = rand.index( Zm[t,], Zt[t,] )
RI_dynsbm[t]= rand.index( list.dynsbm[[1]]$membership[,t], Zt[t,] )
RI_FFS[t] =  rand.index( Zm_FFS[t,], Zt[t,] )
}

c(mean(RI),mean(RI_dynsbm),mean(RI_FFS))

# cat('rand index for inverse gamma',mean(RI),'\n')
# cat('rand index for dynsbm',mean(RI_dynsbm),'\n')
# cat('rand index for FSS',mean(RI_FFS),'\n')
# 
# RI_all = c(RI_all,mean(RI))
# RI2_all = c(RI2_all,mean(RI2))
# RI_FSS_all = c(RI_FSS_all,mean(RI_FFS))
# n_all = c(n_all,mean(n))
# prob_all = c(prob_all,mean(prob))
return(c(mean(RI),mean(RI_dynsbm),mean(RI_FFS)))
}

# df_plot = data.frame(list(n=factor(n_all),prob=factor(prob_all), RI = c(RI_FSS_all,RI_all,RI2_all),
#                            varaince_prior = c(rep('FFS',200),rep('dynsbm',200),
#                                               rep('Gamma',200))))
# 
# df_plot = data.frame(list(n=factor(n_all),prob=factor(prob_all), RI = c(RI_FSS_all,RI_all),
#                           varaince_prior = c(rep('FFS',200),rep('dynsbm',200) )))
# 
# # #
# # ggplot(data=df_plot, aes(x=prob, y=RI, group=varaince_prior)) + theme(legend.position="top")+
# #    geom_boxplot(size=1.5,aes(color = varaince_prior, shape = varaince_prior)) + ylab("RI")+
# #    ggtitle("RI Comparison") + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# df_plot$n = as.character(df_plot$n)
# 
# df_plot$n[df_plot$n=='20']='n=20'
# df_plot$n[df_plot$n=='40']='n=40'
# df_plot$n = factor(df_plot$n)
# 
# ggplot(df_plot, aes(prob, y=RI, fill=varaince_prior)) +
#   geom_boxplot(aes(color = varaince_prior, shape = varaince_prior))+ facet_wrap(vars(factor(n)),scale='free')+ theme(axis.text.x = element_text(colour = 'black', size = 15),
#                                                                                                axis.title.x = element_text(size = 15,
#                                                                                                                            hjust = 0.5, vjust = 0.2)) +
#   theme(axis.text.y = element_text(colour = 'black', size = 15),
#         axis.title.y = element_text(size = 15,
#                                     hjust = 0.5, vjust = 0.2)) +
#   theme(legend.text=element_text(size=15),legend.title=element_text(size=15),
#         legend.key.size = unit(1, 'cm'))+
#   theme(strip.text = element_text(size=15))
# 
# 
# 
# # library(igraph)
# # 
# # A= matrix(rep(0,25),nrow = 5)
# # 
# # A[1,3] = 1
# # A[2,3] = 1
# # A[4,3] =1
# # A[5,3] =1
# # 
# # A= t(A)+A
# # 
# # 
# # G = graph_from_adjacency_matrix(A,mode = "undirected")
# # 
# # plot(G,layout = layout.random,vertex.size = 20)
# #  p2=  ggplot(data=df_change,aes(x=time,y=change))+ geom_line(size=1.5,color='red')
# #  
# #  grid.arrange(p1, p2, nrow = 2,heights=c(0.6,0.4))
# # 
# # change =rep(0,T)
# # for (t in 2:T){
# #   change[t] = sqrt(sum((Xt[[t]]-Xt[[t-1]])^2))
# # }
# # df_change =  data.frame(list(time=1:100, change = change/max(change))) 
# # 
# # RI2 =(RI2-min(RI2))/(max(RI2)-min(RI2))
# # RI =(RI-min(RI))/(max(RI)-min(RI))
# # change =(change-min(change))/(max(change)-min(change))
# # 
# # df_plot = data.frame(list(time=rep(1:100,3), RI = c(RI2,RI,change),
# #                           varaince_prior = c(rep('FSS',100),rep('Inv_Gamma',100),rep('change',100))))
# # 
# #  ggplot(data=df_plot, aes(x=time, y=RI, group=varaince_prior)) + theme(legend.position="right")+
# #    geom_line(size=1.5,aes(color = varaince_prior)) + ylab("RI")+
# #    ggtitle("RI Comparison") + theme(plot.title = element_text(hjust = 0.5))
# #  
# #  
# #  ggplot(data=df_plot, aes(x=varaince_prior, y=RI)) + 
# #    geom_boxplot()
# 
# 
# 
# 
