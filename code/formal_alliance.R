rm(list=ls())
set.seed(2021)
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


file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_network_sym_adaptive.R')
source('mix_DN_adaptive_invgamma.R')
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


MF_list = MP_binary_weighted_adaptive (Y, mean_beta_prior=0, sigma_beta_prior=sqrt(100), 
                                       gap =0.1, d=2 ,max_iter=100, global_prior = 'Cauthy')

Xm = vector("list", T)
for(t in 1:T){
  Xm [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm [[t]][i,] = MF_list$Mean_X[[i]][t,]
  }
}


for ( t in 2:T){
  Xm[[t]] = t(procrustes_r(t(Xm[[t-1]]),t(Xm[[t]]))$B.transformed)
}




Xm2 = vector("list", n)
for(i in 1:n){
  Xm2 [[i]] =matrix(rep(0,T*d),nrow=T)
  for (t in 1:T){
    Xm2 [[i]][t,] = Xm [[t]][i,]
  }
}








label.plot = vector("list", T)

for (t in 1:T){
  label.plot[[t]] = rep(NA,n)
top_ind = which(apply(Y[[t]],1,sum)!=0)
label.plot[[t]][top_ind]=as.character(data_label[top_ind])
}





par(mfrow=c(1,1))
ani.record(reset = TRUE)
dev.new(width=12, height=7)
for (t in 1:T){
 plot_clus_igraph(Xm[[t]],Y[[t]],data_cont,t+1816-1,data_label)
   ani.record()
}
oopts = ani.options(interval = 0.5,ani.res=600) 
saveLatex(ani.replay(), img.name = "record_plot")





par(mfrow=c(2,4))
par(mar=rep(2,4))
### German 1865 - 1872
plot_pre = plot_sub_preprocess(Xm,Y,t1=53,t2=60,data_label = label.plot[[60]])
for (t in 50:57){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,plot_pre$label.plot_non_NA,plot_pre$axis_lim)
}



dev.new()
par(mfrow=c(2,2))
par(mar=rep(2,4))
#1948 -1815 =133
### US 1948 - 1949
plot_pre = plot_sub_preprocess(Xm,Y,t1=133,t2=136,data_label = label.plot[[136]])
for (t in 133:136){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,plot_pre$label.plot_non_NA,plot_pre$axis_lim)
}



dev.new()
par(mfrow=c(2,2))
par(mar=rep(2,4))
#1989 -1815 =174
### US 1989 - 1992
plot_pre = plot_sub_preprocess(Xm,Y,t1=174,t2=177,data_label = label.plot[[177]])
for (t in 174:177){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,plot_pre$label.plot_non_NA,plot_pre$axis_lim)
}



