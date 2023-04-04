
rm(list=ls())
set.seed(2021)
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)
library(extraDistr)
library(ggplot2)
library(LaplacesDemon)


n = 100  

p = 100

T = 100

d = 2  

trend_generate = function(T,d,tau){
  X = matrix(rep(0,T*d),nrow = T)
  initial_mean <- rnorm(d)
  X[1,] = initial_mean
  for (t in 2:T) {
    local_scale = rhcauchy(d, sigma = 1)
    if (tau == 0){
      global_scale = rhcauchy(d, sigma = 1)
    } else{
      global_scale =  rep(tau,d)
    }
    for (k in 1:d){
      X[t,k] <- X[t-1,k]+ rnorm(n=1,mean=0, sd=local_scale[k]*global_scale[k])
    }
  }
  return(X)
}

global_local_prior =function(x, prior){
  n = length(x)
  y = rep(0,n)
  if(prior == 'IG'){
    local_scale = rinvgamma(n, shape = 1)
    global_scale =  rep(1,n)
  } else {
    local_scale = rhcauchy(n, sigma = 1)
    if (prior =='GL'){
      global_scale =  rhcauchy(n, sigma = 1)
    } else {
      global_scale =  rep(1,n)
      local_scale = rgamma(n, shape = 1)
    }
  }
 y=  rnorm(n,mean=0, sd=local_scale*global_scale)
 return(y)
}


xs =seq(-10,10,length.out = 10000)

density_comps <- data.frame(xs)


density_comps['Global_local'] = density(global_local_prior(xs,prior='GL'), n=length(xs), from=-10, to=10)$y

density_comps['Inv_Gamma'] = density(global_local_prior(xs,prior='IG'), n=length(xs), from=-10, to=10)$y

density_comps['halfcauthy'] = density(global_local_prior(xs,prior='HC'), n=length(xs), from=-10, to=10)$y






font_size = 1



density_comps = data.frame(density_comps)

df_plot=NULL

df_plot$x = rep(density_comps$xs,3)
df_plot$value = c(density_comps$`Global_local`,density_comps$`Inv_Gamma`,
                  density_comps$`halfcauthy`)
df_plot$method = c(rep('Global_local',10000),rep('Inv_Gamma',10000),
                   rep('Gamma',10000))

df_plot =data.frame(df_plot)

p <- ggplot(df_plot, aes(x=x,color=method, y=value)) +
  geom_line(aes(linetype =method,shape=method),size=1)+ylim(0,1)+xlim(-2.5,2.5)+
  scale_linetype_manual(values=c("solid","dashed", "dotted","twodash"))+
  theme(axis.text.x = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(size = 15,
                                    hjust = 0.5, vjust = 0.2)) +
  theme(axis.text.y = element_text(colour = 'black', size = 15),
        axis.title.y = element_text(size = 15,
                                    hjust = 0.5, vjust = 0.2)) +
  theme(legend.text=element_text(size=12),legend.title=element_text(size=12))+
  theme(strip.text = element_text(size=15))
p
