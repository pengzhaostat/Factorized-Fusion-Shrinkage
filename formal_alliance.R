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


MF_list = MP_binary_weighted_adaptive (Y,
                                       gap =0.01, d=2 ,max_iter=100, global_prior = 'Cauthy')

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




# Xm2 = vector("list", n)
# for(i in 1:n){
#   Xm2 [[i]] =matrix(rep(0,T*d),nrow=T)
#   for (t in 1:T){
#     Xm2 [[i]][t,] = Xm [[t]][i,]
#   }
# }


Xm2 = vector("list", T)
for(t in 1:T){
  Xm2 [[t]] =matrix(rep(0,n*d),nrow=n)
  for (i in 1:n){
    Xm2 [[t]][i,] = MF_list$Mean_X[[i]][t,]
  }
  Xm2 [[t]] = Xm2 [[t]]/sqrt(apply(Xm2[[t]]^2,1,sum))
}

 heat_matrix = matrix(rep(0,n*(T-1)),ncol=T-1)
# 
# for (t in 2:T){
#   heat_matrix[,t-1]= sqrt(rowSums((Xm[[t]]-Xm[[t-1]])^2))
# }
# heatmap(heat_matrix,Rowv=NA,Colv=NA)



library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(ggExtra)




df_plot = heat_matrix%>%
  as_tibble() %>%
  rowid_to_column(var="Index0") %>%
  gather(key="t", value="Value", -1) %>%
  mutate(year=as.numeric(gsub("V","",t))) 

df_plot$year = df_plot$year+1815


df_plot$continent = rep(data_cont,T-1)

df_plot = arrange(df_plot, t,continent)

df_plot$country_code = rep(1:n,T-1)

# ggplot(df_plot, aes(country_code, year, fill= Value)) +
#   geom_tile()+ scale_fill_gradient(low="white", high="red")+ theme(text = element_text(size = 16))  

# df_plot$continent = rep(data_cont,T-1)
# ggplot(df_plot, aes(country_code, year, fill= Value)) + facet_wrap(vars(continent))+
#   geom_tile()+ scale_fill_distiller(palette = "RdPu")+ theme(text = element_text(size = 16))  





label.plot = vector("list", T)



for (t in 1:T){
  label.plot[[t]] = rep(NA,n)
  top_ind = which(apply(Y[[t]],1,sum)>0)
  label.plot[[t]][top_ind]=as.character(data_label[top_ind])
}



# par(mfrow=c(1,1))
# ani.record(reset = TRUE)
# dev.new(width=12, height=7)
# for (t in 1:T){
#  plot_clus_igraph(Xm[[t]],Y[[t]],data_cont,t+1816-1,label.plot[[t]])
#    ani.record()
# }
# oopts = ani.options(interval = 0.5,ani.res=600)
# saveLatex(ani.replay(), img.name = "record_plot")




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

#data_name = c(1,2,16,18,20,22,26,45,48,49,59,67,92,107,115,157)

data_name = 1:n

dev.new()

par(mfrow=c(2,2))
par(mar=rep(2,4))
#1948 -1815 =133
### US 1948 - 1949
plot_pre = plot_sub_preprocess(Xm,Y,t1=133,t2=136,data_label = label.plot[[137]])
for (t in 133:136){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,label.plot[[137]],plot_pre$axis_lim,data_name)
}




dev.new()
par(mfrow=c(2,2))
par(mar=rep(2,4))
#1989 -1815 =174
### US 1989 - 1992
#data_name = c(1:5,7:n)
data_name = c(1,2,16,18,20,22,26,45,48,49,59,67,92,107,115,157)
plot_pre = plot_sub_preprocess(Xm,Y,t1=174,t2=177,data_label = label.plot[[174]])
for (t in 174:177){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,label.plot[[t]],plot_pre$axis_lim,data_name)
}




dev.new()
par(mfrow=c(1,2))
par(mar=rep(2,4))
#1989 -1815 =174
### US 1989 - 1992
plot_pre = plot_sub_preprocess(Xm,Y,t1=136,t2=177,data_label = label.plot[[177]])

data_name = !is.na(label.plot[[177]]) | !is.na(label.plot[[136]])
for (t in c(136,177)){
  plot_clus_igraph_sub(Xm[[t]],Y[[t]],data_cont,t+1816-1,label.plot[[177]],plot_pre$axis_lim,data_name=which(data_name) )
}

# km = kmeans(Xm[[136]][data_name,],centers = 3,nstart = 50)
# 
# km$cluster
# 
# plot(Xm[[136]][data_name,],col=km$cluster)
# 
# km2 = kmeans(Xm[[177]][data_name,],centers = 3,nstart = 50)
# 
# plot(Xm[[177]][data_name,],col=km2$cluster)
# 
# rand.index(km$cluster,km2$cluster)
