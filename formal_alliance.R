rm(list=ls())
set.seed(2024)
library(Rcpp)
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
library(changepoints)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('mix_DN_adaptive_invgamma.R')
source('change_point_WB_FFS.R')
source('helper.R')
sourceCpp('MP_gibbs_network_sym_adaptive_cpp.cpp')

setwd(paste(file_location,'/data',sep=""))


data = read.csv('alliance_v4.1_by_directed_yearly.csv',sep = ',',header = T)

data_continent = read.csv('COW country codes.csv',header = T)

setwd(file_location)


T = max(data$year)-min(data$year)+1
 
year= data$year-1816+1

n = max(length(unique(data$state_name1)),length(unique(data$state_name2)))

Y = array(0,dim=c(n,n,T))
for ( t in 1:T){
  Y[,,t]= matrix(rep(0,n^2),nrow=n)
}
data_label = factor(unique(data$state_name1))

data_cont =rep(NA,n)

for ( i in 1:n){
  data_cont[i] = data_continent$Continent[data_continent$StateNme==as.character(data_label[i])]
}

data_cont[is.na(data_cont)] ='NA'


for (k in 1:length(data$version4id)){
  Y[which(data_label==data$state_name1[k]),which(data_label==data$state_name2[k]),year[k]] = 1
}

for ( t in 1:T){
  Y[,,t] = t(Y[,,t])+Y[,,t]
  Y[,,t][Y[,,t]>1] = 1
}

index = apply(Y,1,sum)>10

Y = Y[index,index,]

n = sum(index)

data_label = data_label[index]

data_cont = data_cont[index]

d=2


MF_list = MP_binary_weighted_adaptive (Y,  gap =1e-3, d=d ,max_iter=100)

FFS_result = FFS.WBS(xhat =MF_list$Mean_X,betahat=MF_list$mean_beta, obs_num = T)



df_plot <- data.frame(
  Time = 1:T+1815,
  Dval = FFS_result$Dval
)

# Plot the data using ggplot2
font_size <- 16

p <- ggplot(df_plot, aes(x = Time, y = Dval)) +
  geom_line(color = "blue", size = 1.2) +  # Line plot with custom color and size
  geom_point(color = "red", size = 2) +   # Add points for better visibility
  geom_vline(xintercept = 1935, color = "darkgreen", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 1990, color = "purple", linetype = "dashed", size = 1) +
  # Add annotations for the marked points
  annotate("text", x = 1935, y = max(FFS_result$Dval) + 0.5, label = "1935", color = "darkgreen", size = 5, hjust = 0) +
  annotate("text", x = 1990, y = max(FFS_result$Dval) + 0.5, label = "1990", color = "purple", size = 5, hjust = 0) +
  theme_minimal() +                      # Use a clean minimal theme
  theme(
    axis.text.x = element_text(size = font_size, color = "black"),
    axis.text.y = element_text(size = font_size, color = "black"),
    axis.title.x = element_text(size = font_size, hjust = 0.5, vjust = 0.2),
    axis.title.y = element_text(size = font_size, hjust = 0.5, vjust = 0.2),
    legend.text = element_text(size = font_size),
    legend.title = element_text(size = font_size)
  ) +
  labs(
    x = "Year",           # X-axis label
    y = "Changes in estimated probabilities"           # Y-axis label
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = font_size + 4, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = font_size, face = "italic")
  )

print(p)







### Find the top change point

t1 = FFS_result$S[order(FFS_result$Dval,decreasing = T)[1]]



Xm = vector("list", T)
for(t in 1:T){
  Xm [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm [[t]][i,] = MF_list$Mean_X[i,,t]
  }
}

Xm_diff = matrix(rep(0,n*(T-1)),nrow=n)

for ( t in 2:T){
  Xm[[t]] = t(procrustes_r(t(Xm[[t-1]]),t(Xm[[t]]))$B.transformed)
  Xm_diff[,t-1] = apply(Xm[[t]]-Xm[[t-1]],1, function(x) sqrt(sum(x^2)))
}







library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(ggExtra)




label.plot = vector("list", T)

for (t in 1:T){
  label.plot[[t]] = rep(NA,n)
  top_ind = which(apply(Y[,,t],1,sum)>2)
  label.plot[[t]][top_ind]=as.character(data_label[top_ind])
}

#data_name = c(1,6,16,18,19,20,22,25,26,27,30,42,45,48,49,59,60,63,67,88,92)
data_name = c(1,6,16,19,25,26,43,62,90)
data_label[data_name]


#t1 =120
### visualize t1
par(mfrow=c(2,2))
par(mar=rep(2,4))
#1935 -1815 =120
###  1934 - 1937
plot_label = label.plot[[t1+2]]
plot_label[setdiff(1:n,data_name)] =NA
plot_pre = plot_sub_preprocess(Xm,Y,t1=t1-1,t2=t1+2,data_label = plot_label)
for (t in (t1-1):(t1+2)){
  plot_clus_igraph_sub(Xm[[t]],Y[,,t],data_cont,t+1816-1,plot_label,plot_pre$axis_lim,data_label)
}
data_label[order(Xm_diff[,t1],decreasing = T)[1:10]]







