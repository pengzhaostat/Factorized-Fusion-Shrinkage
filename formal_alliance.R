rm(list=ls())
set.seed(2024)
library('rmatio')
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library("ggplot2")
library(reshape2)
library("RColorBrewer")
library(gridExtra)
library(grid)
library(ggpubr)
library(cluster)
library(networkDynamicData)
library(lubridate)
library(ggraph)
library(animation)
library(Bessel)
library(plotly)
library(changepoints)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_sym_adaptive.R')
source('mix_DN_adaptive_invgamma.R')
source('change_point_WB_FFS.R')
source('helper.R')

setwd(paste(file_location,'/data',sep=""))


data = read.csv('alliance_v4.1_by_directed_yearly.csv',sep = ',',header = T)

data_continent = read.csv('COW country codes.csv',header = T)

setwd(file_location)


T = max(data$year)-min(data$year)+1
 
year= data$year-1816+1

n = max(length(unique(data$state_name1)),length(unique(data$state_name2)))


Y = vector("list", T)
for ( t in 1:T){
  Y[[t]] = matrix(rep(0,n^2),nrow=n)
}

data_label = factor(unique(data$state_name1))

data_cont =rep(NA,n)

for ( i in 1:n){
  data_cont[i] = data_continent$Continent[data_continent$StateNme==as.character(data_label[i])]
}

data_cont[is.na(data_cont)] ='NA'


for (k in 1:length(data$version4id)){
  Y[[year[k]]][which(data_label==data$state_name1[k]),which(data_label==data$state_name2[k])] = 1
}

for ( t in 1:T){
  Y[[t]] = t(Y[[t]])+Y[[t]]
  Y[[t]][Y[[t]]>1] = 1
}



d=2


MF_list = MP_binary_weighted_adaptive (Y,
                                       gap =1e-2, d=d ,max_iter=100, global_prior = 'Cauthy')

Xm = vector("list", T)
for(t in 1:T){
  Xm [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm [[t]][i,] = MF_list$Mean_X[[i]][t,]
  }
}


M = 5 # number of random intervals for WBS
delta = 5
intervals = WBS.intervals(M = M, lower = 1, upper = T)

FFS_result = FFS.RDPG(xhat =Xm ,betahat=MF_list$mean_beta, obs_num=T,
                      Alpha = intervals$Alpha, Beta = intervals$Beta, delta=delta)


### Find top 2 change point

t1 = FFS_result$S[order(FFS_result$Dval,decreasing = T)[1]]

t2 = FFS_result$S[order(FFS_result$Dval,decreasing = T)[2]]

for ( t in 2:T){
  Xm[[t]] = t(procrustes_r(t(Xm[[t-1]]),t(Xm[[t]]))$B.transformed)
}







library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(ggExtra)



label.plot = vector("list", T)



for (t in 1:T){
  label.plot[[t]] = rep(NA,n)
  top_ind = which(apply(Y[[t]],1,sum)>0)
  label.plot[[t]][top_ind]=as.character(data_label[top_ind])
}



# dev.new()
# 
# par(mfrow=c(2,4))
# par(mar=rep(2,4))
# ### German 1863 - 1870
# plot_pre = plot_sub_preprocess(Xm,Y,t1=48,t2=55,data_label = label.plot[[50]])
# #
# for (t in 48:55){
#   assign_label = label.plot[[t]]
#   assign_label[c(26,27,28,29,30,31,33)]=NA
#   plot_clus_igraph_sub_8(Xm[[t]],Y[[t]],data_cont,t+1816-1,assign_label,plot_pre$axis_lim)
# }


label.plot = vector("list", T)

for (t in 1:T){
  label.plot[[t]] = rep(NA,n)
  top_ind = which(apply(Y[[t]],1,sum)>2)
  label.plot[[t]][top_ind]=as.character(data_label[top_ind])
}

data_name = c(1,2,16,18,20,22,25,26,27,30,45,48,49,59,67,92,107,115,157)

#data_name = 1:n



### visualize t1
par(mfrow=c(2,2))
par(mar=rep(2,4))
#1935 -1815 =120
###  1935 - 1938
plot_pre = plot_sub_preprocess(Xm,Y,t1=t1,t2=t1+3,data_label = label.plot[[t1+3]])
for (t in t1:(t1+3)){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,label.plot[[t1+3]],plot_pre$axis_lim,data_name)
}


### visualize t2
par(mfrow=c(2,2))
par(mar=rep(2,4))
#1948 -1815 =133
### US 1948 - 1949
plot_pre = plot_sub_preprocess(Xm,Y,t1=t2,t2=t2+3,data_label = label.plot[[t2+3]])
for (t in t2:(t2+3)){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,label.plot[[t2+3]],plot_pre$axis_lim,data_name)
}




